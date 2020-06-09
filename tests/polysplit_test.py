"""Tests for polysplit functions."""
from __future__ import division
from pprint import pprint as pp
import pytest

from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry_polyskel import polysplit, polygraph


def test_skeleton_subpolygons():
    """Test splitting polygon into skeleton polygons"""

    # Make polygon
    poly = Polygon2D.from_array([[2, 0], [2, 2], [1, 1], [0, 2], [0, 0]])

    # Tested Polygons
    test_subpolys = [
        Polygon2D.from_array(
            [[2.        , 0.        ],
             [2.        , 2.        ],
             [1.41421356, 0.58578644]]),
        Polygon2D.from_array(
            [[2.        , 2.        ],
             [1.        , 1.        ],
             [1.        , 0.41421356],
             [1.41421356, 0.58578644]]),
        Polygon2D.from_array(
            [[1.        , 1.        ],
             [0.        , 2.        ],
             [0.58578644, 0.58578644],
             [1.        , 0.41421356]]),
        Polygon2D.from_array(
            [[0.        , 2.        ],
             [0.        , 0.        ],
             [0.58578644, 0.58578644]]),
        Polygon2D.from_array(
            [[0.        , 0.        ],
             [2.        , 0.        ],
             [1.41421356, 0.58578644],
             [1.        , 0.41421356],
             [0.58578644, 0.58578644]])]

    # Split
    subpolys = polysplit.skeleton_subpolygons(poly)

    # Assert
    for subpoly, test_subpoly in zip(subpolys, test_subpolys):
        assert subpoly.is_equivalent(test_subpoly, 1e-7)


def test_perimeter_core_subpolygons_hole_error():
    """Test throwing an exception if hole doesn't get computed in straight skeleton."""
    # Construct a simple rectangle
    poly = Polygon2D.from_array([[0, 0], [6, 0], [6, 8], [0, 8]])
    holes = [Polygon2D.from_array([[2, 2], [4, 2], [4, 6], [2, 6]])]

    # Run method
    with pytest.raises(RuntimeError):
        _ = polysplit.perimeter_core_subpolygons(poly, 1, holes)


def test_perimeter_core_subpolygons():
    """Test splitting perimeter and core subpolygons from polygon."""
    # Construct a simple rectangle
    poly = Polygon2D.from_array([[0, 0], [6, 0], [6, 8], [0, 8]])

    # Make solution polygons (list of polygons)
    test_perims = [
        Polygon2D.from_array(((0.0, 0.0), (6.0, 0.0), (5.0, 1.0), (1.0, 1.0))),
        Polygon2D.from_array(((6.0, 0.0), (6.0, 8.0), (5.0, 7.0), (5.0, 1.0))),
        Polygon2D.from_array(((6.0, 8.0), (0.0, 8.0), (1.0, 7.0), (5.0, 7.0))),
        Polygon2D.from_array(((0.0, 8.0), (0.0, 0.0), (1.0, 1.0), (1.0, 7.0)))]

    test_core = Polygon2D.from_array([[1, 1], [5, 1], [5, 7], [1, 7]])

    # Run method
    perims, cores = polysplit.perimeter_core_subpolygons(poly, 1)

    # Check equality
    for perim, test_perim in zip(perims, test_perims):
        assert perim.is_equivalent(test_perim, 1e-10)

    assert test_core.is_equivalent(cores[0], 1e-10)


def test_complex_perimeter_core_subpolygons():
    """Test splitting perimeter and core subpolygons from polygon."""

    poly = Polygon2D.from_array(
        [[0.7, 0.2], [2, 0], [2, 2], [1, 1], [0, 2], [0, 0]])
    holes = [Polygon2D.from_array([[0.6, 0.6], [1.5, 0.6], [1, 0.8], [0.6, 1.2]]),
             Polygon2D.from_array([[1.1, 0.25], [1.5, 0.25], [1.3, 0.5]])]

    p, c = polysplit.perimeter_core_subpolygons(poly, .1, holes=holes, tol=1e-10)

    # Some very simple tests
    assert len(p) == 13
    assert len(c) == 1


def test_polysplit_hashing():
    """Test that tolerance sets hashing correctly from parent functions."""

    poly = Polygon2D.from_array(
        [[1.0123456789, 3.0123456789],
         [3.0123456789, 5.0123456789],
         [1.0123456789, 5.0123456789]])

    # 7 digit tolerance w/ rounding
    tol = 1e-7
    g = polysplit._skeleton_as_directed_graph(poly, holes=None, tol=tol)
    k = g.ordered_nodes[0].key
    assert k == '(1.0123457, 3.0123457)', k

    # 10 digit tolerance w/o rounding
    tol = 1e-10
    g = polysplit._skeleton_as_directed_graph(poly, holes=None, tol=tol)
    k = g.ordered_nodes[0].key
    assert k == '(1.0123456789, 3.0123456789)', k

    # Test with number x 100
    poly = Polygon2D.from_array(
        [[100.0123456789, 300.0123456789],
         [300.0123456789, 500.0123456789],
         [100.0123456789, 500.0123456789]])

    # 10 digit tolerance
    tol = 1e-7
    g = polysplit._skeleton_as_directed_graph(poly, holes=None, tol=tol)
    k = g.ordered_nodes[0].key
    assert k == '(100.0123457, 300.0123457)', k

    # 10 digit tolerance w/o rounding
    tol = 1e-10
    g = polysplit._skeleton_as_directed_graph(poly, holes=None, tol=tol)
    k = g.ordered_nodes[0].key
    assert k == '(100.0123456789, 300.0123456789)', k

if __name__ == "__main__":
    test_polysplit_hashing()