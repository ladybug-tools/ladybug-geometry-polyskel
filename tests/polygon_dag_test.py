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


TOL = 1e-10


def _cmpstr(item1, item2):
    return '{} vs {}'.format(item1, item2)


def test_dag_noskel():
    """Test the dag with no skeleton"""

    # Points
    pt_array = [[0, 0], [6, 0], [6, 6], [3, 9], [0, 6]]

    # Make the polygon
    polygon = Polygon2D.from_array(pt_array)

    # Make the check cases
    chk_pt_lst = [Point2D.from_array(p) for p in pt_array]

    # Inititalize a DAG object
    d = PolygonDAG()
    vertices = polygon.vertices

    # Add edges to DAG
    [d.add_node(vertices[i], [vertices[i + 1]], True)
        for i in range(len(vertices) - 1)]
    d.add_node(vertices[-1], [vertices[0]], True)

    # Test number
    assert len(chk_pt_lst) == d.num_nodes, _cmpstr(len(chk_pt_lst), d.num_nodes)

    # Test root
    assert d.root._order == 0

    # Test adjacencies are correct
    curr_node = d.root
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

    dag = polyskel._skeleton_as_dag(polygon)
    print('test')

def test_dag_skel():
    """Test the DAG with skeleton"""

    pt_array = [[0, 0], [6, 0], [6, 6], [0, 6]]

    # Make the polygon
    polygon = Polygon2D.from_array(pt_array)

    # Adj matrix to test against
    chk_amtx = [
         [0, 1, 0, 1, 1, 0],  # 0
         [1, 0, 1, 0, 0, 1],  # 1
         [0, 1, 0, 1, 0, 1],  # 2
         [1, 0, 1, 0, 1, 0],  # 3
         [1, 0, 0, 1, 0, 1],  # 4
         [0, 1, 1, 0, 1, 0]]  # 5
    #     0, 1, 2, 3, 4, 5

    dag = polyskel._skeleton_as_dag(polygon)

    assert dag is not None

    # amtx, lbls = dag.adj_matrix, dag.adj_labels

    # # Check the labels are integers in order
    # vertices = polygon.vertices
    # for li in lbls.keys():
    #     assert lbls[li].pt.is_equivalent(vertices[li]), \
    #         _cmpstr(lbls[li].pt, vertices[li])

    # # Check the size
    # assert len(amtx) == len(chk_amtx), _cmpstr(len(amtx), len(chk_amtx))
    # assert len(amtx[0]) == len(chk_amtx[0]), _cmpstr(len(amtx[0]), len(chk_amtx[0]))

    # # Test the labels
    # # Test if the adj matrix is correct
    # for i in range(len(chk_amtx)):
    #     for j in range(len(chk_amtx[0])):
    #         assert amtx[i][j] == chk_amtx[i][j], _cmpstr(amtx[i][j], chk_amtx[i][j])


def test_dag_ccw():
    #  Skeletonize and retrieve DAG from polyskeleton
    #  dag = polyskel.skeletonize(p.to_array(), [], 1e-10)
    pass


if __name__ == "__main__":

    test_dag_noskel()
    #test_dag_skel()

