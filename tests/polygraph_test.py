# coding=utf-8
"""Test the class for handling pathways through the straight skeleton network."""
from __future__ import division
import pytest

from ladybug_geometry_polyskel.polygraph import PolygonDirectedGraph, \
    skeleton_as_directed_graph, skeleton_as_cycle_polygons, _vector2hash
from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry.geometry2d.pointvector import Point2D, Vector2D


def _cmpstr(item1, item2):
    return '{} vs {}'.format(item1, item2)


def test_vector2hash():
    """Test the vector hash method"""
    # Integer vector
    vec = Vector2D(1, 1)
    hash_val = _vector2hash(vec, tol=0)
    assert hash_val == '(1.0, 1.0)'

    # Float with 1e-4 points
    vec = Vector2D(1.1111, 1.1111)
    hash_val = _vector2hash(vec, tol=1e-4)
    assert hash_val == '(1.1111, 1.1111)'

    # Float with 1e-4 points w/ rounding
    vec = Vector2D(1.11116, 1.11116)
    hash_val = _vector2hash(vec, tol=1e-4)
    assert hash_val == '(1.1112, 1.1112)'

    # Round to the ones
    vec = Vector2D(115.11116, 115.11116)
    hash_val = _vector2hash(vec, tol=1)
    assert hash_val == '(115.0, 115.0)'

    # Round to the tenths
    vec = Vector2D(116.11116, 116.11116)
    hash_val = _vector2hash(vec, tol=10)
    assert hash_val == '(120.0, 120.0)'

    # 10 digit tolerance w/o rounding
    vec = Vector2D(1.0123456789, 1.0123456789)
    hash_val = _vector2hash(vec, tol=1e-10)
    assert hash_val == '(1.0123456789, 1.0123456789)'

    # 10 digit tolerance w/ rounding
    vec = Vector2D(1.01234567888, 1.01234567888)
    hash_val = _vector2hash(vec, tol=1e-10)
    assert hash_val == '(1.0123456789, 1.0123456789)'

    # 7 digit tolerance w/ rounding
    vec = Vector2D(1.01234567888, 1.01234567888)
    hash_val = _vector2hash(vec, tol=1e-7)
    assert hash_val == '(1.0123457, 1.0123457)'


def test_vector2hash_non_unity():
    """Test the vector hash method with tolerance that don't have a base of 1."""
    vec = Vector2D(1.115, 1.115)
    hash_val = _vector2hash(vec, tol=0.001)
    assert hash_val == '(1.115, 1.115)'

    vec = Vector2D(1.115, 1.115)
    hash_val = _vector2hash(vec, tol=0.002)
    assert hash_val == '(1.114, 1.114)'

    hash_val = _vector2hash(vec, tol=0.003)
    assert hash_val == '(1.116, 1.116)'

    hash_val = _vector2hash(vec, tol=0.004)
    assert hash_val == '(1.116, 1.116)'

    hash_val = _vector2hash(vec, tol=0.005)
    assert hash_val == '(1.115, 1.115)'

    hash_val = _vector2hash(vec, tol=0.006)
    assert hash_val == '(1.116, 1.116)'

    vec = Vector2D(1.1156, 1.1157)
    hash_val = _vector2hash(vec, tol=0.002)
    assert hash_val == '(1.116, 1.116)'


def test_dg_noskel():
    """Test the dg with no skeleton"""
    # Points
    pt_array = [[0, 0], [6, 0], [6, 6], [3, 9], [0, 6]]

    # Make the polygon
    polygon = Polygon2D.from_array(pt_array)

    # Make the check cases
    chk_pt_lst = [Point2D.from_array(p) for p in pt_array]

    # Inititalize a dg object
    d = PolygonDirectedGraph(1e-5)
    vertices = polygon.vertices

    # Add edges to dg
    for i in range(len(vertices) - 1):
        k = d.add_node(vertices[i], [vertices[i + 1]])
        if i == 0:
            d.outer_root_key = k
    d.add_node(vertices[-1], [vertices[0]])

    # Test number
    assert len(chk_pt_lst) == d.node_count, _cmpstr(len(chk_pt_lst), d.node_count)

    # Test root
    assert d.node(d.outer_root_key)._order == 0

    # Test adjacencies are correct
    curr_node = d.node(d.outer_root_key)
    for chk_pt in chk_pt_lst:
        assert chk_pt.is_equivalent(curr_node.pt, 1e-5), _cmpstr(chk_pt, curr_node.pt)

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


def test_skeleton_as_directed_graph_triangle():
    """Test the skeleton_as_directed_graph function."""
    polygon_verts = [
        [0., 0.],
        [7., 0.],
        [4., 4.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    cycle_polys = skeleton_as_cycle_polygons(polygon, tolerance=1e-5)
    for poly in cycle_polys:
        assert isinstance(poly, Polygon2D)
        assert len(poly.vertices) == 3


def test_skeleton_as_directed_graph_square():
    """Test the skeleton_as_directed_graph function."""
    polygon_verts = [
        [0.,  0.],
        [10., 0.],
        [10., 10.],
        [0., 10.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    cycle_polys = skeleton_as_cycle_polygons(polygon, tolerance=1e-5)
    for poly in cycle_polys:
        assert isinstance(poly, Polygon2D)
        assert len(poly.vertices) == 3


def test_skeleton_as_directed_graph_one_hole():
    """Test skeleton_as_directed_graph using a shape with one hole."""
    polygon_verts = [
        [0., 0.],
        [3., 0.],
        [3., 3.],
        [0., 3.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    hole_verts = [
        [1., 1.],
        [2., 1.],
        [2., 2.],
        [1., 2.]
    ]
    hole = Polygon2D.from_array(hole_verts)

    dg = skeleton_as_directed_graph(polygon, [hole], tolerance=1e-5)
    ext_cycle = dg.exterior_cycle(dg.outer_root_node)
    min_cycle = dg.min_cycle(ext_cycle[1], ext_cycle[0])
    assert len(min_cycle) == 4
    assert min_cycle[-1].pt.x == pytest.approx(0., rel=1e-3)
    assert min_cycle[-1].pt.y == pytest.approx(0., rel=1e-3)
    assert min_cycle[0].pt.x == pytest.approx(3., rel=1e-3)
    assert min_cycle[0].pt.y == pytest.approx(0., rel=1e-3)
    assert not dg.is_intersect_topology

    cycle_polys = skeleton_as_cycle_polygons(polygon, [hole], tolerance=1e-5)
    for poly in cycle_polys:
        assert isinstance(poly, Polygon2D)
        assert len(poly.vertices) == 4 or len(poly.vertices) == 5


def test_skeleton_as_directed_graph_concave_two_holes():
    """Test concave with two holes."""
    polygon_verts = [
        [0.7, 0.2],
        [2, 0],
        [2, 2],
        [1, 1],
        [0, 2],
        [0, 0]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    hole1_verts = [
        [0.6, 1.2],
        [1, 0.8],
        [1.5, 0.6],
        [0.6, 0.6]
    ]
    hole1 = Polygon2D.from_array(hole1_verts)

    hole2_verts = [
        [1.3, 0.5],
        [1.5, 0.25],
        [1.1, 0.25]
    ]
    hole2 = Polygon2D.from_array(hole2_verts)

    cycle_polys = skeleton_as_cycle_polygons(polygon, [hole1, hole2], tolerance=1e-5)
    for poly in cycle_polys:
        assert isinstance(poly, Polygon2D)
        assert 4 <= len(poly.vertices) <= 7


def test_skeleton_as_directed_graph_bad_topology():
    """Test a known case where polyskel returns bad topology"""
    polygon_verts = [
        [-3.10661836926, 21.3923815366],
        [-3.10661836926, 19.7597225480],
        [-7.59797976254, 19.7597225480],
        [-7.59797976254, 21.3923815366],
        [-6.59662650790, 21.3923815366],
        [-6.59662650790, 21.8240139474],
        [-4.11444320173, 21.8240139474],
        [-4.11444320173, 21.3923815366]
    ]
    polygon = Polygon2D.from_array(polygon_verts)
    cycle_polys = skeleton_as_cycle_polygons(polygon, tolerance=0.01)
    assert len(cycle_polys) == 8
    for poly in cycle_polys:
        assert isinstance(poly, Polygon2D)
        assert 3 <= len(poly.vertices) <= 7


def test_skeleton_as_directed_graph_bad_topology_v2():
    """Test another known case where polyskel returns bad topology"""
    polygon_verts = [
        [-3.69996219119, 17.0824478496],
        [-3.69996219119, 15.9029877055],
        [-8.19132358447, 15.9029877055],
        [-8.19132358447, 17.5356466941],
        [-7.18997032982, 17.5356466941],
        [-7.18997032982, 17.9672791049],
        [-4.70778702366, 17.9672791049],
        [-4.70778702366, 17.0824478496]
    ]
    polygon = Polygon2D.from_array(polygon_verts)
    cycle_polys = skeleton_as_cycle_polygons(polygon, tolerance=0.01)
    assert len(cycle_polys) == 8
    for poly in cycle_polys:
        assert isinstance(poly, Polygon2D)
        assert 3 <= len(poly.vertices) <= 8
