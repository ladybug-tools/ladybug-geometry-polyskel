# coding=utf-8
"""Classes for computing straight skeleton for 2D polygons."""
from __future__ import division

from pprint import pprint as pp
from ladybug_geometry_polyskel import polyskel, polysplit
from ladybug_geometry_polyskel.polygraph import \
    PolygonDirectedGraph, _vector2hash

from math import pi

from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry.geometry2d.pointvector import Point2D, Vector2D


TOL = 1e-10


def _cmpstr(item1, item2):
    return '{} vs {}'.format(item1, item2)


def test_dg_noskel():
    """Test the dg with no skeleton"""

    # Points
    pt_array = [[0, 0], [6, 0], [6, 6], [3, 9], [0, 6]]

    # Make the polygon
    polygon = Polygon2D.from_array(pt_array)

    # Make the check cases
    chk_pt_lst = [Point2D.from_array(p) for p in pt_array]

    # Inititalize a dg object
    d = PolygonDirectedGraph(1e-10)
    vertices = polygon.vertices

    # Add edges to dg
    for i in range(len(vertices) - 1):
        k = d.add_node(vertices[i], [vertices[i + 1]])
        if i == 0:
            d.outer_root_key = k
    d.add_node(vertices[-1], [vertices[0]])

    # Test number
    assert len(chk_pt_lst) == d.num_nodes, _cmpstr(len(chk_pt_lst), d.num_nodes)

    # Test root
    assert d.node(d.outer_root_key)._order == 0

    # Test adjacencies are correct
    curr_node = d.node(d.outer_root_key)
    for chk_pt in chk_pt_lst:
        assert chk_pt.is_equivalent(curr_node.pt, TOL), _cmpstr(chk_pt, curr_node.pt)

        # Increment
        curr_node = curr_node.adj_lst[0]

    # Test the adj matrix
    amtx = d.adj_matrix()

    # Adj matrix to test against
    chk_amtx = [
         [0, 1, 0, 0, 0],  # 0
         [0, 0, 1, 0, 0],  # 1
         [0, 0, 0, 1, 0],  # 2
         [0, 0, 0, 0, 1],  # 3
         [1, 0, 0, 0, 0]]  # 4
    #     0, 1, 2, 3, 4

    # Test if the adj matrix is correct
    for i in range(len(chk_amtx)):
        for j in range(len(chk_amtx[0])):
            assert amtx[i][j] == chk_amtx[i][j], _cmpstr(amtx[i][j], chk_amtx[i][j])


def test_dg_skel_rectangle():
    """Test the dg with skeleton"""

    pt_array = [[0, 0], [6, 0], [6, 4], [0, 4]]

    # Make the polygon
    polygon = Polygon2D.from_array(pt_array)

    # Adj matrix to test against
    chk_amtx = [
         [0, 1, 0, 0, 1, 0],  # 0
         [0, 0, 1, 0, 0, 1],  # 1
         [0, 0, 0, 1, 0, 1],  # 2
         [1, 0, 0, 0, 1, 0],  # 3
         [1, 0, 0, 1, 0, 1],  # 4
         [0, 1, 1, 0, 1, 0]]  # 5
    #     0, 1, 2, 3, 4, 5

    chk_lbls = {
        0: '(0.0, 0.0)',
        1: '(6.0, 0.0)',
        2: '(6.0, 4.0)',
        3: '(0.0, 4.0)',
        4: '(2.0, 2.0)',
        5: '(4.0, 2.0)'}

    dg = polysplit._skeleton_as_directed_graph(polygon, [], 1e-10)

    amtx = dg.adj_matrix()
    lbls = dg.adj_matrix_labels()

    # Test the size
    assert len(amtx) == len(chk_amtx), _cmpstr(len(amtx), len(chk_amtx))
    assert len(amtx[0]) == len(chk_amtx[0]), _cmpstr(len(amtx[0]), len(chk_amtx[0]))
    assert len(lbls) == len(chk_lbls)

    # Test the labels
    for key in range(len(lbls)):
        assert lbls[key] == chk_lbls[key], 'key: {} result in '.format(key) + \
            _cmpstr(lbls[key], chk_lbls[key])

    # Test if the adj matrix is correct
    # Flip dict for chk, this ensures if points are deleted or unordered in
    # adj mtx test will still pass
    chk_row_dict = {val: key for key, val in chk_lbls.items()}
    for i in range(len(amtx)):
        key = lbls[i]
        ci = chk_row_dict[key]
        for j in range(len(amtx[0])):
            assert amtx[i][j] == chk_amtx[ci][j], 'at index {},{}: '.format(i, j) + \
                _cmpstr(amtx[i][j], chk_amtx[i][j])


def test_dg_skel_concave():
    """Test the dg with skeleton in concave geom"""

    pt_array = [[0, 0], [6, 0], [6, 6], [3, 4], [0, 6]]

    # Make the polygon
    polygon = Polygon2D.from_array(pt_array)

    # Adj matrix to test against
    chk_amtx = [
         [0, 1, 0, 0, 0, 1, 0, 0],  # 0
         [0, 0, 1, 0, 0, 0, 1, 0],  # 1
         [0, 0, 0, 1, 0, 0, 1, 0],  # 2
         [0, 0, 0, 0, 1, 0, 0, 1],  # 3
         [1, 0, 0, 0, 0, 1, 0, 0],  # 4
         [1, 0, 0, 0, 1, 0, 0, 1],  # 5
         [0, 1, 1, 0, 0, 0, 0, 1],  # 6
         [0, 0, 0, 1, 0, 1, 1, 0]]  # 7
    #     0, 1, 2, 3, 4, 5, 6, 7

    chk_lbls = {
        0: '(0.0, 0.0)',
        1: '(6.0, 0.0)',
        2: '(6.0, 6.0)',
        3: '(3.0, 4.0)',
        4: '(0.0, 6.0)',
        5: '(2.09, 2.09)',
        6: '(3.91, 2.09)',
        7: '(3.0, 1.82)'}

    dg = polysplit._skeleton_as_directed_graph(polygon, [], 1e-2)

    amtx = dg.adj_matrix()
    lbls = dg.adj_matrix_labels()

    # Test the size
    assert len(amtx) == len(chk_amtx), _cmpstr(len(amtx), len(chk_amtx))
    assert len(amtx[0]) == len(chk_amtx[0]), _cmpstr(len(amtx[0]), len(chk_amtx[0]))
    assert len(lbls) == len(chk_lbls)

    # Test the labels
    for key in range(len(lbls)):
        assert lbls[key] == chk_lbls[key], 'key: {} result in '.format(key) + \
            _cmpstr(lbls[key], chk_lbls[key])

    # Test if the adj matrix is correct
    # Flip dict for chk, this ensures if points are deleted or unordered in
    # adj mtx test will still pass
    chk_row_dict = {val: key for key, val in chk_lbls.items()}
    for i in range(len(amtx)):
        key = lbls[i]
        ci = chk_row_dict[key]
        for j in range(len(amtx[0])):
            assert amtx[i][j] == chk_amtx[ci][j], 'at index {},{}: '.format(i, j) + \
                _cmpstr(amtx[i][j], chk_amtx[i][j])


def test_edge_direction():
    """ Tests the bidirection method """

    # Make the polygon
    poly = Polygon2D.from_array([[0, 0], [6, 0], [6, 6], [3, 4], [0, 6]])
    pt_array = poly.vertices

    # Make unidirect graph
    dg = PolygonDirectedGraph(1e-10)

    for i in range(len(pt_array)-1):
        k = dg.add_node(pt_array[i], [pt_array[i+1]])
        if i == 0:
            dg.outer_root_key = k
    dg.add_node(pt_array[-1], [pt_array[0]])

    root = dg.node(dg.outer_root_key)

    # Check
    nodes = dg.ordered_nodes
    for i in range(dg.num_nodes-1):
        assert not dg.is_edge_bidirect(nodes[i], nodes[i+1])

    # Check unidirectionality
    next_node = dg.next_unidirect_node(root)
    assert not dg.is_edge_bidirect(root, next_node)

    # Add bidirectional edge
    dg.add_node(Point2D(0, 0), [Point2D(1, 1)])
    bidir_key = dg.add_node(Point2D(1, 1), [Point2D(0, 0)])

    # Check bidirectionality
    assert dg.is_edge_bidirect(dg.node(bidir_key), root)


def test_exterior_cycle():
    """ Tests the exterior cycles method """

    # Make the polygon
    polygon = Polygon2D.from_array([[0, 0], [6, 0], [6, 4], [0, 4]])
    dg = polysplit._skeleton_as_directed_graph(polygon, [], 1e-10)

    root = dg.node(dg.outer_root_key)
    exterior = dg.exterior_cycle(root)

    for pt, node in zip(polygon.vertices, exterior):
        assert node.pt.is_equivalent(pt, 1e-10)


def test_ccw_angle():
    """ Check edge angle relative to ccw orientation. These tests are to
    ensure the assumed ordering of the DG will provide the correct angle."""

    # Test ccw 1
    theta = Vector2D(1, 0).angle_clockwise(Vector2D(-1, 1)*-1)
    assert abs(theta - pi/4.) < 1e-10, str(theta)

    # Test ccw 2
    theta = Vector2D(-1, 1).angle_clockwise(Vector2D(0, -1)*-1)
    assert abs(theta - pi/4.) < 1e-10, str(theta)

    # Test > 180
    theta = Vector2D(-1, 1).angle_clockwise(Vector2D(0, 1)*-1)
    assert abs(theta - (pi + pi/4.)) < 1e-10, str(theta)

    # Test == 180
    theta = Vector2D(-1, 1).angle_clockwise(Vector2D(-1, 1) * -1)
    assert abs(theta - pi) < 1e-5, str(theta)

    # Test > 270
    theta = Vector2D(-1, 1).angle_clockwise(Vector2D(1, 1)*-1)
    assert abs(theta - (pi + pi / 2.)) < 1e-10, str(theta)


def test_min_ccw_cycle():
    """ Find a closed loop from a PolygonDirectedGraph """

    # Make the polygon
    poly = Polygon2D.from_array([[0, 0], [6, 0], [6, 6], [3, 4], [0, 6]])

    # Make the test cases
    chk_poly = [[0, 0], [6, 0], [3.91, 2.09], [3, 1.82], [2.09, 2.09]]
    chk_poly = Polygon2D.from_array(chk_poly)

    # Skeletonize
    dg = polysplit._skeleton_as_directed_graph(poly, [], 1e-10)

    ref_node = dg.node(dg.outer_root_key)

    next_node = dg.next_unidirect_node(ref_node)
    cycle = dg.min_ccw_cycle(ref_node, next_node)

    cycle_poly = Polygon2D.from_array([n.pt for n in cycle])

    assert cycle_poly.is_equivalent(chk_poly, 1e-2)


def test_vector2hash():
    """Test the vector hash method"""

    # Integer vector
    vec = Vector2D(1, 1)
    hash = _vector2hash(vec, tol=0)
    assert hash == '(1.0, 1.0)', hash

    # Float with 1e-4 points
    vec = Vector2D(1.1111, 1.1111)
    hash = _vector2hash(vec, tol=1e-4)
    assert hash == '(1.1111, 1.1111)', hash

    # Float with 1e-4 points w/ rounding
    vec = Vector2D(1.11116, 1.11116)
    hash = _vector2hash(vec, tol=1e-4)
    assert hash == '(1.1112, 1.1112)', hash

    # Round to the ones
    vec = Vector2D(115.11116, 115.11116)
    hash = _vector2hash(vec, tol=1)
    assert hash == '(115.0, 115.0)', hash

    # Round to the tenths
    vec = Vector2D(116.11116, 116.11116)
    hash = _vector2hash(vec, tol=10)
    assert hash == '(120.0, 120.0)', hash

    # 10 digit tolerance w/o rounding
    vec = Vector2D(1.0123456789, 1.0123456789)
    hash = _vector2hash(vec, tol=1e-10)
    assert hash == '(1.0123456789, 1.0123456789)', hash

    # 10 digit tolerance w/ rounding
    vec = Vector2D(1.01234567888, 1.01234567888)
    hash = _vector2hash(vec, tol=1e-10)
    assert hash == '(1.0123456789, 1.0123456789)', hash

    # 7 digit tolerance w/ rounding
    vec = Vector2D(1.01234567888, 1.01234567888)
    hash = _vector2hash(vec, tol=1e-7)
    assert hash == '(1.0123457, 1.0123457)', hash


def test_sub_polygon_traversal():
    # Graph methods to retrieve subpolygons

    tol = 1e-10
    pts = [[0, 0], [6, 0], [6, 8], [0, 8]]
    poly = Polygon2D.from_array(pts)

    # Test retrieval of exterior cycle form root node
    g = polysplit._skeleton_as_directed_graph(poly, None, tol)
    chk_nodes = [[0, 0], [6, 0], [6, 8], [0, 8], [3, 5], [3, 3]]
    assert len(chk_nodes) == len(g.ordered_nodes)
    for i, n in enumerate(g.ordered_nodes):
        assert n.pt.is_equivalent(Point2D.from_array(chk_nodes[i]), tol)

    # Check root
    assert g.outer_root_key == _vector2hash(Vector2D(0, 0), 1e-5)

    # Get exterior cycle
    exterior = g.exterior_cycle(g.node(g.outer_root_key))
    chk_nodes = [[0, 0], [6, 0], [6, 8], [0, 8]]
    assert len(chk_nodes) == len(exterior)
    for i, n in enumerate(exterior):
        assert n.pt.is_equivalent(Point2D.from_array(chk_nodes[i]), tol)


