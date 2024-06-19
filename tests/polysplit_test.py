# coding=utf-8
"""Test the methods for splitting polygons with a straight skeleton."""
from __future__ import division
import pytest

from ladybug_geometry_polyskel.polysplit import perimeter_core_subpolygons
from ladybug_geometry.geometry2d import Polygon2D


def test_perimeter_core_subpolygons_triangle():
    """Test the perimeter_core_subpolygons function with a triangle."""
    polygon_verts = [
        [0., 0.],
        [8., 0.],
        [4., 4.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 1.0, tolerance=1e-5)
    assert len(perim_polys) == 3
    for poly in perim_polys:
        assert len(poly) == 4
    assert perim_polys[0].area == pytest.approx(5.58578643, rel=1e-3)
    assert perim_polys[1].area == pytest.approx(3.949747468, rel=1e-3)
    assert perim_polys[2].area == pytest.approx(3.949747468, rel=1e-3)
    assert len(core_polys) == 1
    assert len(core_polys[0]) == 3
    assert core_polys[0].area == pytest.approx(2.51471862, rel=1e-3)

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 6.0, tolerance=1e-5)
    assert len(perim_polys) == 3
    for poly in perim_polys:
        assert len(poly) == 3
    assert len(core_polys) == 0


def test_perimeter_core_subpolygons_square():
    """Test the perimeter_core_subpolygons function with a square."""
    polygon_verts = [
        [0.,  0.],
        [10., 0.],
        [10., 10.],
        [0., 10.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 2.0, tolerance=1e-5)
    assert len(perim_polys) == 4
    for poly in perim_polys:
        assert len(poly) == 4
        assert poly.area == pytest.approx(16., rel=1e-3)
    assert len(core_polys) == 1
    assert len(core_polys[0]) == 4
    assert core_polys[0].area == pytest.approx(36., rel=1e-3)

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 6.0, tolerance=1e-5)
    assert len(perim_polys) == 4
    for poly in perim_polys:
        assert len(poly) == 3
        assert poly.area == pytest.approx(25., rel=1e-3)
    assert len(core_polys) == 0

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 5.0, tolerance=1e-5)
    assert len(perim_polys) == 4
    for poly in perim_polys:
        assert len(poly) == 3
        assert poly.area == pytest.approx(25., rel=1e-3)
    assert len(core_polys) == 0


def test_perimeter_core_subpolygons_one_hole():
    """Test perimeter_core_subpolygons using a shape with one hole."""
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

    perim_polys, core_polys = perimeter_core_subpolygons(
        polygon, 0.2, [hole], tolerance=1e-5, flat_core=False)
    assert len(perim_polys) == 8
    for poly in perim_polys:
        assert len(poly) == 4
        assert poly.area == pytest.approx(0.56, rel=1e-3) or \
            poly.area == pytest.approx(0.24, rel=1e-3)
    assert len(core_polys) == 1
    assert len(core_polys[0]) == 2
    assert core_polys[0][0].area == pytest.approx(6.76, rel=1e-3)
    assert core_polys[0][1].area == pytest.approx(1.96, rel=1e-3)

    perim_polys, core_polys = perimeter_core_subpolygons(
        polygon, 1.0, [hole], tolerance=1e-5, flat_core=False)
    assert len(perim_polys) == 8
    for poly in perim_polys:
        assert len(poly) == 4 or len(poly) == 5
        assert poly.area == pytest.approx(1.25, rel=1e-3) or \
            poly.area == pytest.approx(0.75, rel=1e-3)
    assert len(core_polys) == 0


def test_perimeter_core_subpolygons_two_holes():
    """Test perimeter_core_subpolygons with a concave geometry that has two holes."""
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

    perim_polys, core_polys = perimeter_core_subpolygons(
        polygon, 0.1, [hole1, hole2], tolerance=1e-5)
    print(perim_polys)
    print(core_polys)
