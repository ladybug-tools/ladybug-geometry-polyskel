# coding=utf-8
"""Tests for polyskel classes."""
from __future__ import division

from ladybug_geometry.geometry2d import LineSegment2D, Polygon2D
from ladybug_geometry_polyskel import polyskel


def helper_check_lavertex(v1, v2):
    """ Checking equality of different LAVertex properties
    """
    tol = 1e10
    assert (v1.point.x - v2.point.x) < tol and \
           (v1.point.y - v2.point.y) < tol

    assert v1.is_reflex == v2.is_reflex


def helper_assert_polygon_equality(polygon, chk_edges, holes=None, tolerance=0.01):
    """
    Consumes polygons and holes as a list of list of vertices, and the
    corresponding list of skeleton edges for checking. This function compares
    the passed chk_edges with the equivalent produced in the ladybug-geometry
    library and returns a boolean if both are equal.

    Args:
         polygon: list of list of polygon vertices as floats in ccw order.
         chk_edges: list of list of line segments as floats.
         holes: list of list of polygon hole vertices as floats in cw order.

    Returns:
        Boolean of equality.
    """

    tst_edges = polyskel.skeleton_as_edge_list(polygon, holes, tolerance)
    print(len(tst_edges))
    # Check number of edges
    assert len(chk_edges) == len(tst_edges)

    # Check lines
    # Assumes order is constant
    is_equal = True
    for i, (chk_edge, tst_edge) in enumerate(zip(chk_edges, tst_edges)):
        for chk_pt, tst_pt in zip(chk_edge.vertices, tst_edge.vertices):
            for chk_coord, tst_coord in zip(chk_pt, tst_pt):
                is_equal = abs(chk_coord - tst_coord) < tolerance
                if not is_equal:
                    break

    if not is_equal:
        print('Edges not equal btwn chk_edgs:')
        print(chk_edges)
        print('and actual edges:')
        print(tst_edges)
        print('\nSpecifically at ln {}:'.format(i))
        print(chk_edge)
        print('!=')
        print(tst_edge)
        print('\n')
        raise Exception('Equality error.')

    return is_equal


def test_polyskel_triangle():
    """Test simplest geometry, a triangle."""
    polygon_verts = [
        [0., 0.],
        [7., 0.],
        [4., 4.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    # Make actual geom that we already solved
    chk_edges = [
        [(3.8284271247461903, 1.5857864376269049), (4.0, 4.0)],
        [(3.8284271247461903, 1.5857864376269049), (7.0, 0.0)],
        [(3.8284271247461903, 1.5857864376269049), (0.0, 0.0)]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    assert helper_assert_polygon_equality(polygon, chk_edges)


def test_polyskel_square():
    """Test square."""
    polygon_verts = [
        [0.,  0.],
        [10., 0.],
        [10., 10.],
        [0., 10.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    # Make actual geom that we already solved
    chk_edges = [
        [(5.0, 5.0), (0.,   0.)],
        [(5.0, 5.0), (0.,  10.)],
        [(5.0, 5.0), (5.0, 5.0)],
        [(5.0, 5.0), (10., 10.)],
        [(5.0, 5.0), (10.,  0.)]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    assert helper_assert_polygon_equality(polygon, chk_edges)


def test_polyskel_pentagon():
    """Test polygon."""
    polygon_verts = [
        [0.,  0.],
        [10., 0.],
        [10., 10.],
        [5., 15.],
        [0., 10.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    # Make actual geom that we already solved
    chk_edges = [
        [(5.0, 7.9289321881345245), (0.,  10.)],
        [(5.0, 7.9289321881345245), (5.0, 15.)],
        [(5.0, 7.9289321881345245), (5., 7.9289321881345245)],
        [(5.0, 7.9289321881345245), (10., 10.)],
        [(5.0, 5.0),                (5., 7.9289321881345245)],
        [(5.0, 5.0),                (10.,  0.)],
        [(5.0, 5.0),                (0.,   0.)]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    assert helper_assert_polygon_equality(polygon, chk_edges)


def test_polyskel_complex_convex():
    """Test complex convex with many edges."""
    polygon_verts = [
        [0.,  0.],
        [2., -1.],
        [4., -1.5],
        [10., 0.],
        [11., 5.],
        [11., 7.],
        [10., 10.],
        [5., 15.],
        [2., 15.],
        [0., 10.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    # Make actual geom that we already solved
    chk_edges = [
        [(3.8613, 12.2509), (2.0, 15.0)],
        [(3.8613, 12.2509), (5.0, 15.0)],
        [(3.0723, 1.8988), (2.0, -1.0)],
        [(3.0723, 1.8988), (0.0, 0.0)],
        [(4.0, 2.6231), (4.0, -1.5)],
        [(4.0, 2.6231), (3.0723, 1.8988)],
        [(4.5012, 9.1331), (0.0, 10.0)],
        [(4.5012, 9.1331), (3.8613, 12.2509)],
        [(5.3361, 7.1175), (4.5012, 9.1331)],
        [(5.3361, 7.1175), (10.0, 10.0)],
        [(5.3866, 4.399), (10.0, 0.0)],
        [(5.3866, 4.399), (4.0, 2.6231)],
        [(5.5, 6.1075), (5.3361, 7.1175)],
        [(5.5, 6.1075), (11.0, 7.0)],
        [(5.5, 5.5446), (5.5, 6.1075)],
        [(5.5, 5.5446), (11.0, 5.0)],
        [(5.5, 5.5446), (5.3866, 4.399)]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    assert helper_assert_polygon_equality(polygon, chk_edges)


def test_polyskel_simple_concave():
    """Test simplest possible concave polygon: A triangle with poked in side."""
    polygon_verts = [
        [0.,  0.],
        [1.,  0.5],
        [2.,  0.],
        [1.,  1.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    # Make actual geom that we already solved
    chk_edges = [
        [(1.0, 0.7207592200561265), (1.0, 1.0)],
        [(1.0, 0.7207592200561265), (2.0, 0.0)],
        [(1.0, 0.7207592200561265), (1.0, 0.7207592200561265)],
        [(1.0, 0.7207592200561265), (1.0, 0.5)],
        [(1.0, 0.7207592200561265), (0.0, 0.0)]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    assert helper_assert_polygon_equality(polygon, chk_edges)


def test_polyskel_concave():
    """Testa concave shape: a rectangle with side poked in."""
    polygon_verts = [
        [0., 0.],
        [2., 0.],
        [2., 2.],
        [1., 1.],
        [0., 2.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    # make actual geom that we already solved
    chk_edges = [
        [(0.585786437626905, 0.585786437626905), (0.0, 0.0)],
        [(0.585786437626905, 0.585786437626905), (0.0, 2.0)],
        [(1.414213562373095, 0.585786437626905), (2.0, 2.0)],
        [(1.414213562373095, 0.585786437626905), (2.0, 0.0)],
        [(1.0, 0.41421356237309503), (0.585786437626905, 0.585786437626905)],
        [(1.0, 0.41421356237309503), (1.0, 1.0)],
        [(1.0, 0.41421356237309503), (1.414213562373095, 0.585786437626905), ]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    assert helper_assert_polygon_equality(polygon, chk_edges)


def test_polyskel_concave_two_holes():
    """ Test concave with two holes"""
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

    # Make actual geoms we already solved
    chk_edges = [
        [(1.3, 0.5616), (1.3, 0.5)],
        [(1.0, 0.9), (1.0, 1.0)],
        [(1.9193, 0.5193), (1.5, 0.6)],
        [(0.5293, 1.3707), (1.0, 0.9)],
        [(0.5293, 1.3707), (0.6, 1.2)],
        [(0.3, 1.2757), (0.0, 2.0)],
        [(0.3, 1.2757), (0.5293, 1.3707)],
        [(0.3, 0.3977), (0.0, 0.0)],
        [(0.3, 0.3977), (0.3, 1.2757)],
        [(0.3557, 0.3557), (0.3, 0.3977)],
        [(0.3557, 0.3557), (0.6, 0.6)],
        [(0.6873, 0.4021), (0.7, 0.2)],
        [(0.6873, 0.4021), (0.3557, 0.3557)],
        [(0.9298, 0.3836), (0.6873, 0.4021)],
        [(0.9298, 0.3836), (1.3, 0.5616)],
        [(1.0005, 0.2022), (0.9298, 0.3836)],
        [(1.0005, 0.2022), (1.1, 0.25)],
        [(1.7129, 0.1477), (1.0005, 0.2022)],
        [(1.7129, 0.1477), (1.5, 0.25)],
        [(1.7775, 0.2594), (2.0, 0.0)],
        [(1.7775, 0.2594), (1.7129, 0.1477)],
        [(1.7468, 0.3468), (1.9193, 0.5193)],
        [(1.7468, 0.3468), (1.7775, 0.2594)],
        [(1.7468, 0.3468), (1.3, 0.5616)],
        [(1.5888, 1.0073), (2.0, 2.0)],
        [(1.5888, 1.0073), (1.9193, 0.5193)],
        [(1.0659, 0.9), (1.5888, 1.0073)],
        [(1.0659, 0.9), (1.0, 0.8)],
        [(1.0659, 0.9), (1.0, 0.9)]
    ]
    chk_edges = [LineSegment2D.from_array(seg) for seg in chk_edges]

    holes = [hole1, hole2]
    assert helper_assert_polygon_equality(polygon, chk_edges, holes)
