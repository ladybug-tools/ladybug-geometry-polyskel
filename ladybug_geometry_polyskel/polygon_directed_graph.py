# coding=utf-8
"""
Implementation of a dg.
"""

from __future__ import division

# FIXME: temp while prototyping. Do not PR
import sys
lbgeom_path = "/app/ladybug-geometry/"
if lbgeom_path not in sys.path:
    sys.path.insert(0, lbgeom_path)

# Geometry classes
from ladybug_geometry.geometry2d.pointvector import Point2D, Vector2D
from ladybug_geometry.geometry2d.line import LineSegment2D
from ladybug_geometry.geometry2d.ray import Ray2D
from ladybug_geometry import intersection2d

import math


def _vector2hash(vector, tol=4):
    #  FIXME: find a better way of hashing spatial coordinates
    # myhash = "("
    # for i in range(len(vector)):
    #     coordinate = vector[i]
    #     myhash += str(round(coordinate, tol))
    #     if i < len(vector)-1:
    #         myhash += ","
    # myhash += ")"
    myhash = vector.__repr__()
    return myhash


class _Node(object):
    """Private Node class for dg containment"""
    def __init__(self, key, val, order, adj_lst):
        """ Private class to handle nodes in Polygondg.

        Args:
            val: Any python object
            key: Hash of passed object
            order: integer counting order of Node (based on dg propagation)
            adj_lst: list of keys: ['key1', 'key2' ....  'keyn']
            exterior: Boolean value representing if point on exterior edge
        """
        self.key = key
        self.pt = val
        self._order = order
        self.adj_lst = adj_lst

    def num_children(self):
        return len(self.adj_lst)

    def __repr__(self):
        return '{}: {}'.format(self._order, self.key)

    def adj_num(self):
        """ Number of adjacent nodes"""
        return len(self.adj_lst)


class PolygonDirectedGraph(object):
    """An directed graph that represents polygon adjacency relationships
    for points and edges. This class assumes that exterior edges are naked
    (unidirectional) and interior edges are bidirectional.
    """

    def __init__(self):
        self._directed_graph = {}
        self._root = None
        self.num_nodes = 0

    def __repr__(self):
        return '\n'.join((n.__repr__() for n in self.ordered_nodes()))

    @property
    def root(self):
        """ Root node, used for traversal of dg """
        if self._root is None:
            self._root = self.ordered_nodes()[0]
        return self._root

    def get_node(self, key):
        """ Retrieves the pt node based on passed value.

        Args:
            val: A Linept2D object

        Returns:
            List of adjacent Linept2D objects, as keys
        """
        if key in self._directed_graph:
            return self._directed_graph[key]

    def add_node(self, val, adj_lst):
        """ Consumes a polygon point, and computes its key value, and adds it in the
        dg if it doesn't exist. If it does exist it appends adj_lst to existing pt.

        Args:
            val: A Point2D object
            adj_lst: A list of Point2D objects

        Returns:
            The hashed key from the existing or new node.
        """

        # Get key
        key = _vector2hash(val)

        # Get node if it exists
        if key not in self._directed_graph:
            # If key doesn't exist, add to dg
            self.num_nodes += 1
            self._directed_graph[key] = _Node(key, val, self.num_nodes-1, [])

        # Add the adj_lst to dg
        adj_lst = (self.get_node(self.add_node(adj, [])) for adj in adj_lst)

        # Add nodes to node in dg
        self._directed_graph[key].adj_lst += adj_lst

        return key

    def nodes(self):
        """Returns an iterable of pt nodes"""
        return self._directed_graph.values()

    def ordered_nodes(self):
        """Returns an iterable of pt nodes in order of addition"""
        nodes = list(self.nodes())
        nodes.sort(key=lambda v: v._order)
        return nodes

    def adj_matrix(self):
        """ Returns an adjacency matrix of dg.
        1 = adjacency from row node to col node.
        0 = no adjacency.

        Returns:
            N x N square matrix where N is number of nodes.
        """
        nodes = self.ordered_nodes()

        # Initialize amtx with no adjacencies
        amtx = [[0 for i in range(self.num_nodes)]
                for j in range(self.num_nodes)]

        for i in range(self.num_nodes):
            adj_indices = [adj._order for adj in nodes[i].adj_lst]
            for adj_idx in adj_indices:
                amtx[i][adj_idx] = 1

        return amtx

    def adj_matrix_labels(self):
        """ Returns a dictionary where label key corresponds to index in adj_matrix
        and value is node key"""
        return {i: node.key for i, node in enumerate(self.ordered_nodes())}

    @classmethod
    def is_edge_bidirect(cls, node1, node2):
        """ Are two nodes bidirectional.

        Args:
            node1: _Node object
            node2: _Node object

        Returns:
            True if node1 and node2 are in each other's adjacency list,
            else False.
        """
        return node1.key in (n.key for n in node2.adj_lst) and \
            node2.key in (n.key for n in node1.adj_lst)

    @classmethod
    def get_next_unidirect_node(cls, node):
        """ Retrieves the first unidirectional point adjacent
        to consumed point. They define an exterior or naked edge.

        Args:
            node: _Node

        Returns:
            Next node that defines unidirectional edge, or None if all
            adjacencies are bidirectional.
        """

        # Check bidirectionality
        next_node = None
        for _next_node in node.adj_lst:
            if not cls.is_edge_bidirect(node, _next_node):
                next_node = _next_node
                break

        return next_node

    def get_exterior_cycle(self, cycle_root):
        """ Retreives exterior boundary. If there is a no continuous outer boundary
        cycle from the node it will return None.

        This method assumes that exterior edges are naked (unidirectional) and
        interior edges are bidirectional.

        Args:
            cycle_root: Starting _Node in exterior cycle.
        Returns:
            List of nodes on exterior if a cycle exists, else None.
        """

        # Get the first exterior edge
        curr_node = cycle_root
        next_node = self.get_next_unidirect_node(curr_node)

        if not next_node:
            return None

        ext_cycle = [curr_node]
        while next_node.key != cycle_root.key:
            ext_cycle.append(next_node)
            next_node = self.get_next_unidirect_node(next_node)
            if not next_node:
                return None

        return ext_cycle

    def get_smallest_closed_cycles(self):
        """ Traverses a directed graph of connected straight skeleton edges and
        produces a list of the smallest individual polygons defined by the edges. This
        is achieved by looping through the exterior edges, and identifying the closed
        loop with the smallest counter-clockwise angle of rotation between edges. Since
        the exterior edges of a polygon split by a straight skeleton will always result
        in either a split or edge event, of the interior skeleton, this will identify the
        smallest polygon nested in the directed graph.

        # TODO: Holes?

        Returns:
            A list of polygon point arrays.
        """

        polygon_lst = []

        # TODO: Better way to get the exterior root node? Define in DG class?

        # Get continous exterior nodes list.
        for node in self.ordered_nodes():
            ext_nodes = self.get_exterior_cycle(node)
            if ext_nodes is not None:
                break

        # Add first node to ensure complete cycle
        ext_nodes += [ext_nodes[0]]
        for i, ext_node in enumerate(ext_nodes[:-1]):
            cycle = [ext_node]
            next_node = ext_nodes[i + 1]
            cycle = self.min_ccw_cycle(ext_node, next_node, next_node.adj_lst, cycle)
            polygon_lst.append(cycle)

        return polygon_lst

    def min_ccw_cycle(self, ref_node, next_node, adj_lst, cycle):
        """
        Recursively identifes most counter-clockwise adjacent node and returns closed
        loop.

        Args:
            ref_node: The first node, for first edge.
            next_node: The node next to ref_node that constitues a polygon edge.
            adj_lst: List of adjacent nodes to next_node.
            cycle: Current list of nodes that will form a polygon.

        Returns:
            A list of nodes that form a polygon.
        """

        # Base case: loop is completed
        if next_node.key == cycle[0].key:
            return cycle

        cycle.append(next_node)

        # Get current edge direction vector
        edge_dir = Vector2D(*(next_node.pt - ref_node.pt).to_array())

        # Initialize values for comparison
        min_theta = float("inf")
        min_node = None

        # Identify the node with the smallest ccw angle
        for adj_node in adj_lst:

            # Make sure this node isn't backtracking traversal
            if adj_node.key == cycle[-2].key:
                continue

            # Get next edge direction vector
            next_edge_dir = Vector2D(*(adj_node.pt - next_node.pt).to_array())

            # Flip the next_edge and calculate theta
            theta = edge_dir.angle_clockwise(next_edge_dir * -1)

            if theta < min_theta:
                min_theta = theta
                min_node = adj_node

        return self.min_ccw_cycle(next_node, min_node, min_node.adj_lst, cycle)
