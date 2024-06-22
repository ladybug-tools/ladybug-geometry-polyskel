# coding=utf-8
"""Test the methods for splitting polygons with a straight skeleton."""
from __future__ import division
import pytest

from ladybug_geometry.geometry2d import Polygon2D

from ladybug_geometry_polyskel.offset import offset_polygon, offset_perimeter_polygons


def test_offset_square():
    """Test the offset_polygon function with a square."""
    polygon_verts = [
        [0.,  0.],
        [10., 0.],
        [10., 10.],
        [0., 10.]
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    core_polys = offset_polygon(polygon, 2.0, tolerance=1e-5)
    assert len(core_polys) == 1
    assert len(core_polys[0]) == 4
    assert core_polys[0].area == pytest.approx(36., rel=1e-3)

    core_polys = offset_polygon(polygon, -2.0, tolerance=1e-5)
    assert len(core_polys) == 1
    assert len(core_polys[0]) == 4
    assert core_polys[0].area == pytest.approx(196., rel=1e-3)

    perim_polys = offset_perimeter_polygons(polygon, 2.0, tolerance=1e-5)
    assert len(perim_polys) == 4
    for poly in perim_polys:
        assert len(poly) == 4
        assert poly.area == pytest.approx(16., rel=1e-3)

    perim_polys = offset_perimeter_polygons(polygon, -2.0, tolerance=1e-5)
    assert len(perim_polys) == 4
    for poly in perim_polys:
        assert len(poly) == 4
        assert poly.area == pytest.approx(24., rel=1e-3)
