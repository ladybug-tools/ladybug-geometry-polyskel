# coding=utf-8
"""Functions for splitting a polygon into subpolygons based on skeleton topology."""
from __future__ import division

from ladybug_geometry_polyskel import polyskel
from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry.geometry2d.line import LineSegment2D
from ladybug_geometry.geometry3d.pointvector import Vector3D
from ladybug_geometry_polyskel.polygraph import PolygonDirectedGraph, _vector2hash

POLYSKELETON_ERROR_MSG = \
    'The straight skeleton topology calculated for this geometry is incorrect.'


def skeleton_subpolygons(polygon, holes=None, tol=1e-8, recurse_limit=3000,
                         print_recurse=False):
    """Compute the polygon straight skeleton as a list of polygon point arrays.

    Args:
        polygon: list of Polygon2D objects.
        holes: list of Polygon2D objects representing holes in cw order. (Default: None).
        tol: Tolerance for point equivalence. (Default: 1e-10).
        recurse_limit: optional parameter to limit the number of cycles looking for
            polygons in skeleton graph. (Default: 3000).
        print_recurse: optional boolean to print cycle loop cycles, for
            debugging. (Default: False).

    Returns:
        A list of the straight skeleton subpolygons as Polygon2D objects.
    """
    subs1, subs2 = perimeter_core_subpolygons(
        polygon, float('inf'), holes=holes, tol=tol, recurse_limit=recurse_limit,
        print_recurse=print_recurse)

    return subs1 + subs2


def perimeter_core_subpolygons(polygon, distance, holes=None, tol=1e-8,
                               recurse_limit=3000, print_recurse=False):
    """Compute the perimeter and core sub-polygons from the polygon straight skeleton.

    Args:
        polygon: A Polygon2D to split into perimeter and core subpolygons.
        distance: Distance in model units to offset perimeter subpolygon.
        holes: A list of Polygon2D objects representing holes in the
            polygon. (Default: None).
        tol: Tolerance for point equivalence. (Default: 1e-10).
        recurse_limit: optional parameter to limit the number of cycles looking for
            polygons in skeleton graph. (Default: 3000).
        print_recurse: optional boolean to print cycle loop cycles, for
            debugging. (Default: False).

    Returns:
        A tuple with two lists:

        * perimeter_sub_polys -- A list of perimeter subpolygons as Polygon2D objects.

        * core_sub_polys -- A list of core subpolygons as Polygon2D objects. In the
            event of a core sub-polygon with a hole, a list with be returned with
            the first item being a boundary and successive items as hole polygons.
    """
    # Initialize sub polygon lists
    perimeter_sub_polys, core_sub_polys, hole_sub_polys = [], [], []

    # Compute the straight skeleton of the polygon and get exterior edges
    dg = _skeleton_as_directed_graph(polygon, holes, tol)

    if holes is not None:
        holes_exist = all([_hole_exists_in_skeleton(hole, dg) for hole in holes])
        if not holes_exist:
            raise RuntimeError(POLYSKELETON_ERROR_MSG + ' Error calculating the '
                               'holes for the polygon. Try changing the geometry '
                               'of the hole.')

    # Make list of roots for graph traversal of just the exterior cycles
    root_keys = [dg.outer_root_key] + dg.hole_root_keys

    for i, root_key in enumerate(root_keys):
        # Compute the polygons on the perimeter of the polygon
        _perimeter_sub_polys, _perimeter_sub_dg = \
            _split_perimeter_subpolygons(dg, distance, root_key, tol, recurse_limit,
            print_recurse)

        # Add perimeter subpolygons
        perimeter_sub_polys.extend(_perimeter_sub_polys)

        # Compute the polygons on the core of the polygon
        _core_sub_polys = _exterior_cycles_as_polygons(_perimeter_sub_dg)

        # Reverse polygons b/c graph traversals will put interior polys as cw
        _core_sub_polys = [p.reverse() for p in _core_sub_polys]

        # Add to lists
        if i == 0:
            # If polygon is not a hole, get rid of largest area
            core_sub_polys.extend(_core_sub_polys[:-1])
        else:
            # If polygon is a hole, get rid of smallest area
            hole_sub_polys.extend(_core_sub_polys[1:])

    if len(hole_sub_polys) > 0:
        # Remake cores w/ holes
        core_sub_polys = _add_holes_to_polygons(core_sub_polys, hole_sub_polys)

    return perimeter_sub_polys, core_sub_polys


def _skeleton_as_directed_graph(_polygon, holes, tol):
    """Compute the straight skeleton of a polygon as a PolygonDirectedGraph.

    Args:
        polygon: polygon as Polygon2D.
        holes: holes as list of Polygon2Ds.
        tol: Tolerance for point equivalence.

    Returns:
        A PolygonDirectedGraph object.
    """

    if _polygon.is_clockwise:
        # Exterior polygon must be in counter-clockwise order.
        _polygon = _polygon.reverse()
    if holes is not None:
        # Interior polygon must be in counter-clockwise order.
        for i, hole in enumerate(holes):
            if hole.is_clockwise:
                holes[i] = hole.reverse()

    # Make directed graph
    dg = PolygonDirectedGraph(tol=tol)

    # Reverse order to ensure cw order for input
    holes_array = [] if holes is None else [hole.to_array() for hole in holes]
    polygon = _polygon.reverse()  # flip to cw order
    slav = polyskel._SLAV(polygon.to_array(), holes_array, tol)

    # Get the exterior polygon coordinates making sure to flip back to ccw
    for lav in slav._lavs:
        vertices = list(reversed(list(lav)))

        # Get order, account for the fact we flip outer polygons back to ccw order
        is_lav_cw_order = not Polygon2D.from_array(
            [v.point for v in vertices]).is_clockwise

        # Start with last point to be consistent with order of point input, and then
        # add rest of vertices in order.
        for j in range(len(vertices) - 1):
            curr_v = vertices[j]
            next_v = vertices[j + 1]
            k = dg.add_node(curr_v.point, [next_v.point], exterior=True)

            # Add roots
            if j == 0:
                if is_lav_cw_order:  # Outer polygon
                    if dg.outer_root_key is None:
                        dg.outer_root_key = k
                    else:
                        raise RuntimeError('Outer root key is already assigned. '
                                           'Cannot reassign outer root key.')
                else:
                    dg.hole_root_keys.append(k)

        # Close loop
        dg.add_node(vertices[-1].point, [vertices[0].point], exterior=True)

    # Compute the skeleton
    subtree_list = polyskel._skeletonize(slav)

    for subtree in subtree_list:
        event_pt = subtree.source
        for sink_pt in subtree.sinks:
            # Add a bidirectional edge to represent skeleton edges
            dg.add_node(sink_pt, [event_pt])
            dg.add_node(event_pt, [sink_pt], exterior=False)

    return dg


def _split_polygon_graph(node1, node2, distance, poly_graph, tol,
                         recurse_limit=3000, print_recurse=False):
    """Split the PolygonDirectedGraph by an edge offset at a distance.

    Args:
        node1: A node defining the edge to offset for splitting. If the graph
            defines a counter-clockwise polygon, this node should represent
            the left hand-most point.
        node2: A node defining the edge to offset for splitting. If the graph
            defines a counter-clockwise polygon, this node should represent
            the right hand-most point.
        distance: The distance to offset the splitting edge, in model units.
        poly_graph: A PolygonDirectedGraph representing a single polygon.
        tol: Tolerance for point equivalence.
        recurse_limit: optional parameter to limit the number of cycles looking for
            polygons in skeleton graph. (Default: 3000).
        print_recurse: optional boolean to print cycle loop cycles, for
            debugging. (Default: False).

    Returns:
        A list of nodes representing the split polygon adjacent to the exterior
        edge defined by node1 and node2.
    """

    ext_seg = LineSegment2D.from_end_points(node1.pt, node2.pt)

    # Compute normal facing into the polygon
    ext_arr = ext_seg.v.to_array()
    ext_seg_v = Vector3D(ext_arr[0], ext_arr[1], 0)
    normal = ext_seg_v.cross(Vector3D(0, 0, -1)).normalize() * distance

    # Move segment by normal to get offset segment
    offset_seg = ext_seg.move(normal)

    poly_graph = PolygonDirectedGraph.from_point_array(
        [n.pt for n in poly_graph], tol=tol)
    poly_graph.outer_root_key = poly_graph.ordered_nodes[0].key

    # Update graph by intersecting offset segment with other edges
    poly_graph.intersect_graph_with_segment(offset_seg)

    # Get the minimum cycle. Since we start at the exterior edge, this will
    # return the perimeter offset.
    root_node = poly_graph.node(poly_graph.outer_root_key)
    next_node = root_node.adj_lst[0]

    split_poly_nodes = poly_graph.min_ccw_cycle(
        root_node, next_node, recurse_limit=recurse_limit, print_recurse=print_recurse)

    return split_poly_nodes


def _split_perimeter_subpolygons(dg, distance, root_key, tol, recurse_limit=3000,
                                 print_recurse=False):
    """
    Split polygons in the PolygonDirectedGraph into perimeter subpolygons.

    Args:
        dg: A PolygonDirected Graph.
        distance: The distance to offset the perimeter in model units.
        root_key: A key representing any node that identifies the polygon to be split.
            This will be used as the starting point for the graph traversal.
        tol: A number representing the tolerance for point equivalence.
        recurse_limit: optional parameter to limit the number of cycles looking for
            polygons in skeleton graph. (Default: 3000).
        print_recurse: Flag to print loop cycles. (Default: False).

    Returns:
        A tuple with two lists:

        A list of perimeter subpolygons as Polygon2D objects.

        A PolygonDirectedGraph of just the perimer subpolygons. The edges of
            this graph that are not bidirectional (exterior edges) define the
            exterior of the polygon and it's core subpolygons.
    """

    # Start an empty directed graph, and extract exterior cycle associated
    # with the root_key
    perimeter_subpolygons = []
    perimeter_sub_dg = PolygonDirectedGraph(tol=tol)
    exterior = dg.exterior_cycle(dg.node(root_key))

    if exterior is None:
        raise RuntimeError(POLYSKELETON_ERROR_MSG + ' Error traversing the skeleton '
                           'for the following polygon:\n{}'.format(
                               [n.pt.to_array() for n in exterior]))

    # Cycle through exterior nodes and split polygons by distance
    for exterior_node in exterior:

        next_node = exterior_node.adj_lst[0]

        # Find the smallest polygon defined by the exterior node
        min_ccw_poly_graph = PolygonDirectedGraph.min_ccw_cycle(
            exterior_node, next_node, recurse_limit=recurse_limit,
            print_recurse=print_recurse)


        # Offset edge from specified distance, and cut a perimeter polygon
        split_poly_graph = _split_polygon_graph(
            exterior_node, next_node, distance, min_ccw_poly_graph, tol,
            recurse_limit=recurse_limit, print_recurse=print_recurse)

        # Store perimeter nodes in DG
        pts = [node.pt for node in split_poly_graph]
        for i in range(len(pts) - 1):
            perimeter_sub_dg.add_node(pts[i], [pts[i + 1]])
        perimeter_sub_dg.add_node(pts[-1], [pts[0]])

        # Store perimeter points in list as polygon
        perimeter_subpolygons.append(Polygon2D(pts))

    return perimeter_subpolygons, perimeter_sub_dg


def _exterior_cycles_as_polygons(dg):
    """Convert and sort exterior cycles in a PolygonDirectedGraph to a list of polygons.

    Args:
        dg: A PolygonDirectedGraph.

    Returns:
        A list of Polygon2D objects sorted by area.
    """

    # Get all exterior cycles as polygons
    exterior_polys = [Polygon2D([node.pt for node in cycle])
                      for cycle in dg.exterior_cycles]

    # Sort by area and return
    return sorted(exterior_polys, key=lambda p: p.area)


def _add_holes_to_polygons(polys, hole_polys):
    """Add holes to a list of polygons.

    Note that this method will add holes to the polygon only if it is inside or
    partially inside the polygon.

    Args:
        polys: A list of Polygon2Ds to add holes to.
        hole_polys: A list of Polygon2Ds representing the holes.

    Returns:
        A list of the Polygon2D polygons with holes.
    """
    holed_polys = []
    for poly in polys:
        # Check to see if hole is not outside of core
        # Often, the hole is only partially inside the core, so
        # we can't check if it's inside
        inside_holes = [hole for hole in hole_polys
                        if not poly.is_polygon_outside(hole)]

        if len(inside_holes) == 0:
            holed_polys.append(poly)
        else:
            holed_polys.append([poly] + inside_holes)

    return holed_polys


def _hole_exists_in_skeleton(polygon, dg):
    """Check if polygon is in directed graph skeleton.

    This method is based on the PolygonDirectedGraph.polygon_exist method,
    with modifications to detect incorrect straight skeletons.

    Args:
        polygons: A Polygon2D object.
        dg: A PolygonDirectedGraph.

    Return:
        True if exists, else False.
    """
    vertices_loop = list(polygon.vertices)
    vertices_loop = vertices_loop + [vertices_loop[0]]

    for i in range(len(vertices_loop) - 1):
        pt1 = vertices_loop[i]
        pt2 = vertices_loop[i + 1]

        if not dg.pt_exists(pt1):
            return False

        node1 = dg.node(_vector2hash(pt1, dg._tol))
        node2 = dg.node(_vector2hash(pt2, dg._tol))
        key_lst = [n.key for n in node1.adj_lst]

        if node2.key in key_lst:
            return False

        if len(key_lst) == 1:
            # Check to make sure hole is connected to skeleton
            return False

    return True
