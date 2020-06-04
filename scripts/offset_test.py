# coding=utf-8
"""Classes for computing straight skeleton for 2D concave polygons."""
from __future__ import division

from pprint import pprint as pp
import math
from ladybug_geometry_polyskel import polyskel

from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry.geometry2d.pointvector import Point2D, Vector2D
from ladybug_geometry.geometry3d.pointvector import Vector3D
from ladybug_geometry_polyskel.polygon_directed_graph import \
    PolygonDirectedGraph, _vector2hash


def test_polygon_offset_inward():
    """Test the offset method"""

    # Construct a simple rectangle
    poly = [[0, 0], [4, 0], [4, 6], [0, 6]]
    poly = Polygon2D.from_array(poly)

    # Make solution polygon (list of polygons)
    chk_off = Polygon2D.from_array([[1, 1], [3, 1], [3, 5], [1, 5]])

    # Run method
    offset = polyskel.offset(poly, 1, [])

    assert offset[0].is_equivalent(chk_off, 1e-2)


def test_polygon_offset_outward():
    """Test the offset method"""

    # Construct a simple rectangle
    poly = [[1, 1], [3, 1], [3, 5], [1, 5]]
    poly = Polygon2D.from_array(poly)

    # Make solution polygon (list of polygons)
    chk_off = Polygon2D.from_array([[0, 0], [4, 0], [4, 6], [0, 6]])

    # Run method
    offset = polyskel.offset(poly, -1, [])

    #assert offset[0].is_equivalent(chk_off, 1e-2)


def test_convex_angles():
    """ Test simple rectangle angles"""
    poly = Polygon2D.from_array([[0, 0], [4, 0], [4, 6], [0, 6]])
    chk_angles = [90, 90, 90, 90]

    angles = polyskel.interior_angles(poly, radian=False)
    angles = list(angles)

    for i in range(len(angles)):
        assert abs(chk_angles[i] - angles[i]) < 1e-10

    # Test pentagon
    poly = Polygon2D.from_array(
        [[0, 0], [4, 0], [4, 6], [2, 8], [0, 6]])

    chk_angles = [90, 135, 90, 135, 90]

    angles = polyskel.interior_angles(poly, radian=False)
    angles = list(angles)

    pp(angles)
    for i in range(len(angles)):
        assert abs(chk_angles[i] - angles[i]) < 1e-10


def test_concave_angles():
    """ Test simple concave angles"""
    poly = Polygon2D.from_array(
        [[0, 0], [4, 0], [4, 6], [2, 4], [0, 6]])

    conv_theta = math.acos(2 / math.sqrt(8)) / math.pi * 180
    conc_theta = 180. + (2 * conv_theta)
    chk_angles = [90, conv_theta, conc_theta, conv_theta, 90]

    angles = polyskel.interior_angles(poly, radian=False)
    angles = list(angles)

    for i in range(len(angles)):
        assert abs(chk_angles[i] - angles[i]) < 1e-10


if __name__ == "__main__":

    # Inward offset
    #test_polygon_offset_inward()
    test_polygon_offset_outward()
    #test_convex_angles()
    #test_concave_angles()

