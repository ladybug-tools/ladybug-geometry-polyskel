# coding=utf-8
"""Functions for splitting a polygon into sub-polygons based on skeleton topology."""
from __future__ import division

from ladybug_geometry.geometry2d import LineSegment2D, Polygon2D
from ladybug_geometry.geometry3d import Vector3D, Face3D

from .polygraph import PolygonDirectedGraph, \
    skeleton_as_directed_graph


def perimeter_core_subfaces(face, distance, tolerance=1e-5):
    """Compute perimeter and core Face3Ds using the straight skeleton of an input Face3D.

    Args:
        face: A Face3D to split into perimeter and core sub-faces.
        distance: Distance to offset perimeter sub-faces.
        tolerance: Tolerance for point equivalence. (Default: 1e-10).

    Returns:
        A tuple with two items.

        * perimeter_sub_faces -- A list of Face3Ds for perimeter sub-faces.

        * core_sub_faces -- A list of Face3D for core sub-faces.
    """
    # get core and perimeter sub-polygons
    perimeter_sub_polys, core_sub_polys = perimeter_core_subpolygons(
        face.boundary_polygon2d, distance, face.hole_polygon2d, tolerance, False)
    # convert the perimeter polygons into Face3Ds
    perimeter_sub_faces = []
    for poly in perimeter_sub_polys:
        verts_3d = tuple(face.plane.xy_to_xyz(pt) for pt in poly.vertices)
        perimeter_sub_faces.append(Face3D(verts_3d, face.plane))
    # convert the core polygons into Face3Ds
    core_sub_faces = []
    for poly_group in core_sub_polys:
        all_verts_3d = [tuple(face.plane.xy_to_xyz(pt) for pt in poly.vertices)
                        for poly in poly_group]
        if len(all_verts_3d) == 1:  # no holes in the shape
            core_sub_faces.append(Face3D(all_verts_3d[0], face.plane))
        else:
            core_sub_faces.append(Face3D(all_verts_3d[0], face.plane, all_verts_3d[1:]))
    return perimeter_sub_faces, core_sub_faces


def perimeter_core_subpolygons(
        polygon, distance, holes=None, tolerance=1e-5, flat_core=True):
    """Compute perimeter and core sub-polygons using the polygon straight skeleton.

    Args:
        polygon: A Polygon2D to split into perimeter and core sub-polygons.
        distance: Distance to offset perimeter sub-polygons.
        holes: A list of Polygon2D objects representing holes in the
            polygon. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-10).
        flat_core: A boolean to note whether the core_sub_polys should be returned
            as a flat list or as a nested list of lists for boundaries and
            holes within each geometry.

    Returns:
        A tuple with two items.

        * perimeter_sub_polys -- A list of Polygon2Ds for perimeter sub-polygons.

        * core_sub_polys -- When flat_core is True, this will be a list of Polygon2D
            where each Polygon2D represents a loop in the core geometry. When
            flat_core is False, this will be a list of lists where each sub-list
            contains Polygon2Ds and represents one core geometry. The sub-list
            will have has at least one Polygon2D in it and, in the event that
            a core geometry has holes, there will be multiple Polygon2Ds in the
            sub-list. The first item in the list will be the outer boundary of
            the geometry and successive items represent hole polygons.
    """
    # initialize sub polygon lists
    perimeter_sub_polys, core_sub_polys = [], []
    # compute the straight skeleton of the polygon as a directed graph
    dg = skeleton_as_directed_graph(polygon, holes, tolerance)

    # traverse the directed graph to get all of the perimeter polygons
    _perimeter_sub_dg = None
    root_keys = [dg.outer_root_key] + dg.hole_root_keys
    for root_key in root_keys:
        # compute the polygons on the perimeter of the polygon
        _perimeter_sub_polys, _perimeter_sub_dg = _split_perimeter_subpolygons(
            dg, distance, root_key, tolerance, _perimeter_sub_dg)
        perimeter_sub_polys.extend(_perimeter_sub_polys)  # collect perimeter sub-polys

    # compute the polygons on the core of the polygon
    core_sub_polys = _exterior_cycles_as_polygons(_perimeter_sub_dg, tolerance)
    if not flat_core and len(core_sub_polys) != 0:  # remake cores w/ holes
        core_sub_polys = group_boundaries_and_holes(core_sub_polys, tolerance)

    return perimeter_sub_polys, core_sub_polys


def _split_perimeter_subpolygons(dg, distance, root_key, tol, core_graph=None):
    """Split polygons in the PolygonDirectedGraph into perimeter sub-polygons.

    Args:
        dg: A PolygonDirectedGraph object that comes from the skeleton_as_directed_graph
            function and contains separations for how a geometry will be split.
        distance: The distance to offset the perimeter in model units.
        root_key: A string key for a node in the input dg that identifies the
            perimeter loop to be traversed. This will be used as the starting
            point for graph traversal.
        tol: A number representing the tolerance for point equivalence.
        core_graph: An optional PolygonDirectedGraph that is being built to extract
            the core of the shape after all perimeter sub-polygons are extracted.
            This is used when the input dg contains holes and this method for
            _split_perimeter_subpolygons must be run multiple times for the
            boundary and each hole. (Default: None).

    Returns:
        A tuple with two lists:

        * perimeter_subpolygons: A list of perimeter sub-polygons as Polygon2D objects.

        * core_graph: A PolygonDirectedGraph with just the perimeter sub-polygons.
            The edges of this graph that are not bidirectional (exterior edges)
            define the exterior of the polygon and it's core sub-polygons.
    """
    # set up a list lto collect the perimeter polygons and the core graph
    perimeter_subpolygons = []
    core_graph = PolygonDirectedGraph(tol=tol) if core_graph is None else core_graph
    exter_cycle = dg.exterior_cycle(dg.node(root_key))

    # cycle through exterior nodes and get split polygons using the distance
    for i, base_node in enumerate(exter_cycle):
        try:
            next_node = exter_cycle[i + 1]
        except IndexError:  # the last edge of the polygon
            next_node = exter_cycle[0]
        # find the smallest polygon defined by the exterior node
        min_ccw_poly_graph = dg.min_cycle(next_node, base_node)
        # offset edge from specified distance, and cut a perimeter polygon
        split_poly_graph = _split_polygon_graph(
            base_node, next_node, distance, min_ccw_poly_graph, tol)
        # store the perimeter nodes in the core_graph
        pts = [node.pt for node in split_poly_graph]
        for i in range(len(pts) - 1):
            core_graph.add_node(pts[i], [pts[i + 1]])
        # store perimeter points in list as polygon
        perimeter_subpolygons.append(Polygon2D(pts))

    return perimeter_subpolygons, core_graph


def _split_polygon_graph(node1, node2, distance, poly_graph, tol):
    """Split the PolygonDirectedGraph by an edge offset at a distance.

    Args:
        node1: A node defining the edge to offset for splitting. If the graph
            defines a counter-clockwise polygon, this node should represent
            the left hand-most point.
        node2: A node defining the edge to offset for splitting. If the graph
            defines a counter-clockwise polygon, this node should represent
            the right hand-most point.
        distance: The distance to offset the splitting edge, in model units.
        poly_graph: A list of Nodes representing a PolygonDirectedGraph for
            a single polygon.
        tol: Tolerance for point equivalence.

    Returns:
        A list of nodes representing the split polygon adjacent to the exterior
        edge defined by node1 and node2.
    """
    # use the nodes to create a line segment
    ext_seg = LineSegment2D.from_end_points(node1.pt, node2.pt)
    # compute normal facing into the polygon
    ext_arr = ext_seg.v.to_array()
    ext_seg_v = Vector3D(ext_arr[0], ext_arr[1], 0)
    normal = ext_seg_v.cross(Vector3D(0, 0, -1)).normalize() * distance

    # move the segment by normal to get offset segment
    offset_seg = ext_seg.move(normal)
    poly_graph = PolygonDirectedGraph.from_point_array(
        [n.pt for n in poly_graph], tol=tol)
    ordered_nodes = poly_graph.ordered_nodes
    poly_graph.outer_root_key = ordered_nodes[0].key
    root_node = ordered_nodes[0]
    goal_node = ordered_nodes[-1]

    # update graph by intersecting offset segment with other edges
    poly_graph.intersect_graph_with_segment(offset_seg)

    # get the minimum cycle through the graph connecting the nodes
    # since we start at the exterior edge, this will include the perimeter offset
    split_poly_nodes = poly_graph.min_cycle(root_node, goal_node, ccw_only=True)

    return split_poly_nodes


def _exterior_cycles_as_polygons(dg, tol):
    """Convert and sort exterior cycles in a PolygonDirectedGraph to a list of polygons.

    Args:
        dg: A PolygonDirectedGraph.
        tol: The tolerance at which point equivalence is set.

    Returns:
        A list of Polygon2D objects sorted by area.
    """
    ext_cycles = dg.exterior_cycles()
    ext_polygons = []
    for cycle in ext_cycles:
        ext_poly = Polygon2D([node.pt for node in cycle])
        ext_poly = ext_poly.remove_colinear_vertices(tol)
        if ext_poly.is_self_intersecting:
            ext_polygons.extend(ext_poly.split_through_self_intersection(tol))
        else:
            ext_polygons.append(ext_poly)
    return ext_polygons


def group_boundaries_and_holes(polygons, tol):
    """Group polygons by whether they are contained within another.

    Args:
        polygons: A list of Polygon2Ds to be grouped according to boundary and
            holes within those boundaries.

    Returns:
        A list of lists where each sub-list
            contains Polygon2Ds and represents one core geometry. The sub-list
            will have has at least one Polygon2D in it and,, in the event that
            a core geometry has holes, there will be multiple Polygon2Ds in the
            sub-list. The first item in the list will be the outer boundary of
            the geometry and successive items represent hole polygons.
    """
    # first check to be sure that there isn't just one polygon
    if len(polygons) == 1:
        return [polygons]
    # sort the polygons by area and separate base polygon from the remaining
    polygons = sorted(polygons, key=lambda x: x.area, reverse=True)
    base_poly = polygons[0]
    remain_polys = list(polygons[1:])

    # merge the smaller polygons into the larger polygons
    merged_polys = []
    while len(remain_polys) > 0:
        merged_polys.append(_match_holes_to_poly(base_poly, remain_polys, tol))
        if len(remain_polys) > 1:
            base_poly = remain_polys[0]
            del remain_polys[0]
        elif len(remain_polys) == 1:  # lone last Polygon2D
            merged_polys.append([remain_polys[0]])
            del remain_polys[0]
    return merged_polys


def _match_holes_to_poly(base_poly, other_polys, tol):
    """Attempt to merge other polygons into a base polygon as holes.

    Args:
        base_poly: A Polygon2D to serve as the base.
        other_polys: A list of other Polygon2D objects to attempt to merge into
            the base_poly as a hole. This method will remove any Polygon2D
            that are successfully merged into the output from this list.

    Returns:
        A list of Polygon2D where the first item is the base_poly and successive
        items represent holes in this geometry.
    """
    holes = []
    more_to_check = True
    while more_to_check:
        for i, r_poly in enumerate(other_polys):
            if base_poly.polygon_relationship(r_poly, tol) == 1:
                holes.append(r_poly)
                del other_polys[i]
                break
        else:
            more_to_check = False
    return [base_poly] + holes
