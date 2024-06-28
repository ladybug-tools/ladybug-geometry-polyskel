# coding=utf-8
"""Functions for splitting a polygon into sub-polygons based on skeleton topology."""
from __future__ import division

from ladybug_geometry.geometry2d import LineSegment2D, Polygon2D
from ladybug_geometry.geometry3d import Vector3D, LineSegment3D, Face3D

from .polygraph import PolygonDirectedGraph, skeleton_as_directed_graph


def perimeter_core_subfaces(face, distance, tolerance=1e-5):
    """Compute perimeter and core Face3Ds using the straight skeleton of an input Face3D.

    Args:
        face: A Face3D to split into perimeter and core sub-faces.
        distance: Distance to offset perimeter sub-faces.
        tolerance: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        A tuple with two items.

        * perimeter_sub_faces -- A list of Face3Ds for perimeter sub-faces.

        * core_sub_faces -- A list of Face3D for core sub-faces.
    """
    # get core and perimeter sub-polygons
    perimeter, core = perimeter_core_subpolygons(
        face.boundary_polygon2d, distance, face.hole_polygon2d, tolerance, False)
    # convert the polygons into Face3Ds
    return _polygons_to_face3d(face, perimeter, core)


def perimeter_core_subfaces_and_skeleton(face, distance, tolerance=1e-5):
    """Compute perimeter and core Face3Ds using the straight skeleton of an input Face3D.

    Args:
        face: A Face3D to split into perimeter and core sub-faces.
        distance: Distance to offset perimeter sub-faces.
        tolerance: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        A tuple with two items.

        * skeleton -- A list of LineSegment3D for the segments of the straight
            skeleton being used to generate core/perimeter polygons.

        * perimeter_sub_faces -- A list of Face3Ds for perimeter sub-faces.

        * core_sub_faces -- A list of Face3D for core sub-faces.
    """
    # get core and perimeter sub-polygons
    skeleton_2d, perimeter, core = perimeter_core_subpolygons_and_skeleton(
        face.boundary_polygon2d, distance, face.hole_polygon2d, tolerance)
    perimeter_sub_faces, core_sub_faces = _polygons_to_face3d(face, perimeter, core)
    # convert the skeleton segments into LineSegment3Ds
    skeleton = []
    for seg in skeleton_2d:
        verts_3d = tuple(face.plane.xy_to_xyz(pt) for pt in seg.vertices)
        skeleton.append(LineSegment3D.from_end_points(*verts_3d))
    return skeleton, perimeter_sub_faces, core_sub_faces


def perimeter_core_subpolygons(
        polygon, distance, holes=None, tolerance=1e-5, flat_core=True):
    """Compute perimeter and core sub-polygons using the polygon straight skeleton.

    Args:
        polygon: A Polygon2D to split into perimeter and core sub-polygons.
        distance: Distance to offset perimeter sub-polygons.
        holes: A list of Polygon2D objects representing holes in the
            polygon. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-5).
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
        core_sub_polys = Polygon2D.group_boundaries_and_holes(core_sub_polys, tolerance)
    return perimeter_sub_polys, core_sub_polys


def perimeter_core_subpolygons_and_skeleton(
        polygon, distance, holes=None, tolerance=1e-5):
    """Get perimeter and core sub-polygons with line segments for the straight skeleton.

    Args:
        polygon: A Polygon2D to split into perimeter and core sub-polygons.
        distance: Distance to offset perimeter sub-polygons.
        holes: A list of Polygon2D objects representing holes in the
            polygon. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        A tuple with three items.

        * skeleton -- A list of LineSegment2D for the segments of the straight
            skeleton being used to generate core/perimeter polygons.

        * perimeter_sub_polys -- A list of Polygon2Ds for perimeter sub-polygons.

        * core_sub_polys -- A list of lists where each sub-list contains Polygon2Ds
            and represents one core geometry. The sub-list will have has at
            least one Polygon2D and, in the event that a core geometry has holes,
            there will be multiple Polygon2Ds in the sub-list. The first item in
            the list will be the outer boundary of the geometry and successive
            items represent hole polygons.
    """
    # initialize sub polygon lists
    perimeter_sub_polys, core_sub_polys = [], []
    # compute the straight skeleton of the polygon as a directed graph
    dg = skeleton_as_directed_graph(polygon, holes, tolerance)
    skeleton = dg.connection_segments
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
    if len(core_sub_polys) != 0:  # remake cores w/ holes
        core_sub_polys = Polygon2D.group_boundaries_and_holes(core_sub_polys, tolerance)
    return skeleton, perimeter_sub_polys, core_sub_polys


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

        # proceed to extract the polygon and build up the core graph
        if min_ccw_poly_graph is not None:
            # everything is good, offset edge by distance and cut a perimeter polygon
            split_poly_graph = _split_polygon_graph(
                base_node, next_node, distance, min_ccw_poly_graph, tol)
            # store the perimeter nodes in the core_graph
            pts = [node.pt for node in split_poly_graph]
            for i in range(len(pts) - 1):
                core_graph.add_node(pts[i], [pts[i + 1]])
            # store perimeter points in list as polygon
            perimeter_subpolygons.append(Polygon2D(pts))

        else:  # handle the case of unsuccessful minimum cycle polygon extraction
            # if we are just missing a skeleton segment, there may still be hope
            if len(next_node.adj_lst) == 1:  # connect the node to the nearest interior
                rel_nodes, node_dists = [], []
                for node in dg.nodes:
                    if not node.exterior:
                        rel_nodes.append(node)
                        node_dists.append(next_node.pt.distance_to_point(node.pt))
                near_nodes = [n for _, n in sorted(zip(node_dists, rel_nodes),
                                                   key=lambda pair: pair[0])]
                dg.add_adj(next_node, [near_nodes[0].pt])
                dg.add_adj(near_nodes[0], [next_node.pt])  # add node di-directionally
                min_ccw_poly_graph = dg.min_cycle(next_node, base_node)

            if min_ccw_poly_graph is not None:  # connecting the node worked!
                split_poly_graph = _split_polygon_graph(
                    base_node, next_node, distance, min_ccw_poly_graph, tol)
                pts = [node.pt for node in split_poly_graph]
                for i in range(len(pts) - 1):
                    core_graph.add_node(pts[i], [pts[i + 1]])
                perimeter_subpolygons.append(Polygon2D(pts))

            else:  # no hope of extracting a perimeter sub-polygon here
                msg = 'Failed to compute sub-polygons between ' \
                    'Node {} and Node {}.'.format(base_node, next_node)
                print(msg)
                # close the core graph with the exterior segment and move on
                core_graph.add_node(next_node.pt, [base_node.pt])

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


def _polygons_to_face3d(face, perimeter, core):
    """Convert lists of perimeter and core Polygon2D into Face3D."""
    # convert the perimeter polygons into Face3Ds
    perimeter_sub_faces = []
    for poly in perimeter:
        verts_3d = tuple(face.plane.xy_to_xyz(pt) for pt in poly.vertices)
        perimeter_sub_faces.append(Face3D(verts_3d, face.plane))
    # convert the core polygons into Face3Ds
    core_sub_faces = []
    for poly_group in core:
        all_verts_3d = [tuple(face.plane.xy_to_xyz(pt) for pt in poly.vertices)
                        for poly in poly_group]
        if len(all_verts_3d) == 1:  # no holes in the shape
            core_sub_faces.append(Face3D(all_verts_3d[0], face.plane))
        else:
            core_sub_faces.append(Face3D(all_verts_3d[0], face.plane, all_verts_3d[1:]))
    return perimeter_sub_faces, core_sub_faces
