# coding=utf-8
"""Classes for computing straight skeleton for 2D polygons."""
from __future__ import division

from pprint import pprint as pp
from ladybug_geometry_polyskel import polyskel
from ladybug_geometry_polyskel.polygon_dag import PolygonDAG

# FIXME: temp while prototyping. Do not PR
import sys
lbgeom_path = "/app/ladybug-geometry/"
if lbgeom_path not in sys.path:
    sys.path.insert(0, lbgeom_path)
# FIXME: temp import
import numpy as np

from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry.geometry2d.pointvector import Point2D, Vector2D
from ladybug_geometry.geometry2d.line import LineSegment2D


def helper_skeletonize(polygon, holes, tol=1e-10):
    # A short skeletonize fx for testing.
    # Skeletonize fx
    dag = None
    polygon = polygon[::-1]
    slav = polyskel._SLAV(polygon, [], 1e-10)
    prioque = polyskel._EventQueue()

    for lav in slav:
        for vertex in lav:
            v = vertex.next_event()
            prioque.put(v)

    while not (prioque.empty() or slav.empty()):
        i = prioque.get()  # vertex a, b is self or next vertex
        # Handle edge or split events.
        # arc: subtree(event.intersection_point, event.distance, sinks)
        # events: updated events with new vertex
        if isinstance(i, polyskel._EdgeEvent):
            if not i.vertex_a.is_valid or not i.vertex_b.is_valid:
                # Discarded outdated edge event
                continue
            (arc, events) = slav.handle_edge_event(i)
        elif isinstance(i, polyskel._SplitEvent):
            if not i.vertex.is_valid:
                # Discarded outdated split event
                continue
            (arc, events) = slav.handle_split_event(i)

        prioque.put_all(events)

    return dag


def test_polygon_dict():
    """Test the dictionary"""
    # Make the polygon
    polygon = Polygon2D.from_array(
        [[0, 0], [6, 0], [6, 6], [3, 9], [0, 6]])

    chk_key_lst = [
        [[0, 0], [6, 0]],
        [[6, 0], [6, 6]],
        [[6, 6], [3, 9]],
        [[3, 9], [0, 6]],
        [[0, 6], [0, 0]]
    ]

    # Init dictionary
    d = PolygonDAG()
    segments = polygon.segments

    # Add edges to dictionary
    [d.add_segment_adj(segments[i], [segments[i + 1]])
        for i in range(len(segments) - 1)]
    d.add_segment_adj(segments[-1], [segments[0]])

    # Test the key index
    for i, chk_key in enumerate(chk_key_lst):
        line = LineSegment2D.from_array(chk_key)
        adjs = d.get_segment_adj(line)
        print(adjs)
        #print(d.keys[i])
        #assert LineSegment2D.from_array(chk_key) is d.keys[i]

    # Test the dictionary
    pass

def test_dag_dict():
    """Test the DAG."""
    # Make the polygon
    polygon = Polygon2D.from_array(
        [[0, 0], [6, 0], [6, 6], [3, 9], [0, 6]])

    chk_edges = [
        [[3., 4.75735931], [3., 9.]],
        [[3., 4.75735931], [6., 6.]],
        [[3., 3.], [6., 0.]],
        [[3., 3.], [0., 0.]],
        [[3., 4.75735931], [0., 6.]],
        [[3., 4.75735931], [3., 4.75735931]],
        [[3., 4.75735931], [3., 3.]]]

    chk_amtx = [
        [[0., 0.]]
    ]

    skel = polyskel.skeletonize(polygon, [])
    #dag = helper_skeletonize(polygon, [])


    # Test if hash results in all the same edges in dag.dict
    assert True

    # Test if the edges are correct in the dag.dict
    assert True


def test_dag_amtx():
    # Test if the adj matrix is correct
    assert True


def test_dag_ccw():
    #  Skeletonize and retrieve DAG from polyskeleton
    #  dag = polyskel.skeletonize(p.to_array(), [], 1e-10)
    pass


if __name__ == "__main__":

    # Test if hash works
    #test_dag_dict()
    test_polygon_dict()
