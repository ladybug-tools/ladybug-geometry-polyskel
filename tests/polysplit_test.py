# coding=utf-8
"""Test the methods for splitting polygons with a straight skeleton."""
from __future__ import division
import pytest

from ladybug_geometry.geometry2d import Polygon2D
from ladybug_geometry.geometry3d import LineSegment3D, Face3D

from ladybug_geometry_polyskel.polysplit import perimeter_core_subpolygons, \
    perimeter_core_subfaces, perimeter_core_subfaces_and_skeleton


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
    assert len(perim_polys) == 13
    assert len(core_polys) == 2
    assert sum(poly.area for poly in perim_polys + core_polys) == \
        pytest.approx(polygon.area - hole1.area - hole2.area, rel=1e-3)


def test_perimeter_core_subpolygons_signed_zero():
    """Test a case with a signed zero that was failing."""
    polygon_verts = [
        (0.0, -0.0),
        (-66.365025351950564, -64.111893276245013),
        (30.121656801816368, -64.111893276245013),
        (83.372463222004896, 3.5354591123066541),
        (22.759770196994474, -23.559633773719245)
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 5.0, tolerance=1e-5)
    assert len(perim_polys) == 5
    for poly in perim_polys:
        assert len(poly) == 4
    assert len(core_polys) == 1
    assert len(core_polys[0]) == 5
    assert sum(poly.area for poly in perim_polys + core_polys) == \
        pytest.approx(polygon.area, rel=1e-3)


def test_perimeter_core_concave_intersect():
    """Test a case where an offset depth creates 4 intersections with a concave shape."""
    polygon_verts = [
        (342.55419132538884, -10.017859160602118),
        (315.77976924081884, 15.204009322638644),
        (304.52663639962361, 49.624331671429445),
        (323.99843057415865, 59.584601981416057),
        (330.48830165015289, 35.040359098751637),
        (345.50412057008964, 17.152954528454099),
        (377.18511138328233, 73.378466870092268),
        (396.15898612775987, 53.204057854261976),
        (351.58292353598546, -18.52304727135888)
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    perim_polys, core_polys = perimeter_core_subpolygons(polygon, 11.0, tolerance=1e-5)
    assert len(perim_polys) == 9
    assert len(core_polys) == 2
    assert sum(poly.area for poly in perim_polys + core_polys) == \
        pytest.approx(polygon.area, rel=1e-3)


def test_perimeter_core_self_intersecting_core():
    """Test a case where the raw extracted core polygon is self-intersecting."""
    polygon_verts = [
        (-37.735447736091537, -14.307192311918065),
        (-37.735447736091537, -7.0708061795016199),
        (-33.572580444862723, -1.9353063435931741),
        (-22.017705814068719, -1.3906321185725812),
        (-16.143005244203756, -9.2884083813711751),
        (-29.837671473292943, -7.2653326884375442),
        (-25.752614785638496, -14.65734002800273),
        (-24.624361033810125, -19.325976242464954),
        (-33.494769841288345, -18.470059603146879)
    ]
    polygon = Polygon2D.from_array(polygon_verts)
    hole_verts = [
        (-35.012076610988572, -13.37346506902562),
        (-29.254091946485165, -14.696245329789917),
        (-33.261338030565241, -9.7552720028173976)
    ]
    hole = Polygon2D.from_array(hole_verts)

    perim_polys, core_polys = perimeter_core_subpolygons(
        polygon, 1.4, [hole], tolerance=1e-5)
    assert len(perim_polys) == 12
    assert len(core_polys) == 2
    assert sum(poly.area for poly in perim_polys + core_polys) == \
        pytest.approx(polygon.area - hole.area, rel=1e-3)
    for core_poly in core_polys:
        assert not core_poly.is_self_intersecting


def test_perimeter_core_subpolygons_infinite_loop_in_core():
    """Test perimeter_core_subpolygons with a core polygraph having a 3-way intersection.

    This can create an infinite loop if the method searching for the core cycle
    has no way to mark already-explored nodes.
    """
    polygon_verts = [
        (192.26116569601979, -80.969790470169656),
        (206.26116569601979, -76.969790470169656),
        (232.26116569601979, -80.969790470169656),
        (232.26116569601979, -40.969790470169656),
        (212.26116569601979, -60.969790470169656),
        (192.26116569601979, -40.969790470169656)
    ]
    polygon = Polygon2D.from_array(polygon_verts)

    hole1_verts = [
        (204.26116569601979, -68.969790470169656),
        (204.26116569601979, -56.969790470169649),
        (212.26116569601979, -64.969790470169656),
        (222.26116569601979, -68.969790470169656)
    ]
    hole1 = Polygon2D.from_array(hole1_verts)

    hole2_verts = [
        (214.26116569601979, -75.969790470169656),
        (218.26116569601979, -70.969790470169656),
        (222.26116569601979, -75.969790470169656)
    ]
    hole2 = Polygon2D.from_array(hole2_verts)

    perim_polys, core_polys = perimeter_core_subpolygons(
        polygon, 1.0, [hole1, hole2], tolerance=1e-5)
    assert len(perim_polys) == 13
    assert len(core_polys) == 2
    assert core_polys[0].polygon_relationship(core_polys[1], 0.01) == 1


def test_perimeter_core_subfaces():
    """Test the perimeter_core_subfaces function."""
    polygon_verts = [[
        [0.,  0., 0.],
        [10., 0., 0.],
        [10., 10., 0.],
        [0., 10., 0.]
    ]]
    face = Face3D.from_array(polygon_verts)

    perim_faces, core_faces = perimeter_core_subfaces(face, 2.0, tolerance=1e-5)
    assert len(perim_faces) == 4
    for f in perim_faces:
        assert len(face) == 4
        assert f.area == pytest.approx(16., rel=1e-3)
    assert len(core_faces) == 1
    assert len(core_faces[0]) == 4
    assert core_faces[0].area == pytest.approx(36., rel=1e-3)
    assert sum(f.area for f in perim_faces + core_faces) == \
        pytest.approx(face.area, rel=1e-3)


def test_perimeter_core_subfaces_and_skeleton():
    """Test the perimeter_core_subfaces_and_skeleton function"""
    polygon_verts = [[
        [0.,  0., 0.],
        [10., 0., 0.],
        [10., 10., 0.],
        [0., 10., 0.]
    ]]
    face = Face3D.from_array(polygon_verts)

    skeleton, perim_faces, core_faces = perimeter_core_subfaces_and_skeleton(
        face, 2.0, tolerance=1e-5)
    assert len(perim_faces) == 4
    for f in perim_faces:
        assert len(face) == 4
        assert f.area == pytest.approx(16., rel=1e-3)
    assert len(core_faces) == 1
    assert len(core_faces[0]) == 4
    assert core_faces[0].area == pytest.approx(36., rel=1e-3)
    assert sum(f.area for f in perim_faces + core_faces) == \
        pytest.approx(face.area, rel=1e-3)
    assert len(skeleton) == 8
    for seg in skeleton:
        assert isinstance(seg, LineSegment3D)
        assert seg.length == pytest.approx(10.0, rel=1e-3) or \
            seg.length == pytest.approx(7.0710678, rel=1e-3)
