# coding=utf-8
"""
Implementation of a Directed Graph.
"""

from __future__ import division

# Geometry classes
from ladybug_geometry.geometry2d.line import LineSegment2D
from ladybug_geometry import intersection2d
from math import log10


def _vector2hash(vector, tol):
    """Hashes spatial coordinates for use in dictionary.

    Args:
        vector: A Vector2D object.
        tol: floating point precision tolerance.

    Returns:
        Hash of vector as a string of rounded coordinates.
    """
    try:
        rtol = int(log10(tol)) * -1
    except ValueError:
        rtol = 0

    return str((round(vector.x, rtol), round(vector.y, rtol)))


class _Node(object):
    """Private class to handle nodes in PolygonDirectedGraph.

        Args:
            val: A Point2D object.
            key: Hash of Point2D object.
            order: integer counting order of Node (based on dg propagation)
            adj_lst: list of keys adjacent to this node.
            exterior: Node boundary condition. None if not set by user, else True
                or False according to user.

        Properties:
        * pt: A Point2D object.
        * key: Hash of Point2D object.
        * adj_lst: A list of keys adjacent to this node.
        * exterior: Node boundary condition. None if not set by user, else True or
            False according to user.
        * adj_count: Number of adjacent nodes to this node.
    """
    __slots__ = ('key', 'pt', '_order', 'adj_lst', 'exterior')

    def __init__(self, key, val, order, adj_lst, exterior):
        """Initialize _Node"""

        self.key = key
        self.pt = val
        self._order = order
        self.adj_lst = adj_lst
        # Potentially change exterior to data (similar to networkX)
        # and pass conditional function to get_exterior
        # this resolves redundancy between unidirect and exterior
        # node/edge properties.
        self.exterior = exterior

    def __repr__(self):
        return '{}: {}'.format(self._order, self.key)

    @property
    def adj_count(self):
        """Number of adjacent nodes"""
        return len(self.adj_lst)


class PolygonDirectedGraph(object):
    """A directed graph for point and edge adjacency relationships.

    This class assumes that exterior edges are naked (unidirectional), oriented
    counter-clockwise, and interior edges are bidirectional.

    Args:
        tol: floating point precision used for hashing points.

    Properties:
        * outer_root_key: Root key for outside exterior boundary (i.e not holes).
        * hole_root_keys: List of root keys for inside exterior boundary (holes).
        * num_nodes: Number of nodes in graph.
        * nodes: An iterable of nodes in graph.
        * ordered_nodes: An interable of nodes in graph in order they were added.
        * exterior_cycles: A list of unidirectional edge arrays.
    """

    def __init__(self, tol):
        """Initialize a PolygonDirectedGraph."""
        self._directed_graph = {}
        self._tol = tol
        self.outer_root_key = None
        self.hole_root_keys = []

    def __repr__(self):
        s = ''
        for n in self.ordered_nodes:
            s += '{}, [{}]\n'.format(
                n.pt.to_array(),
                ', '.join([str(_n.pt.to_array()) for _n in n.adj_lst]))
        return s

    @classmethod
    def from_polygon(cls, polygon, tol):
        """Generate a directed graph from a polygon.

        Args:
            polygon: A Polygon2D object.
        """
        return cls.from_point_array(polygon.vertices, tol, loop=True)

    @classmethod
    def from_point_array(cls, point_array, tol, loop=True):
        """Generate a directed graph from a 1-dimensional array of points.

        Args:
            point_array: Array of Point2D objects.
            loop: Optional parameter to connect 1d array
            tol: Tolerance for point equivalence.
        """

        dg = cls(tol)
        for i in range(len(point_array) - 1):
            dg.add_node(point_array[i], [point_array[i+1]], exterior=True)

        if loop:
            dg.add_node(point_array[-1], [point_array[0]], exterior=True)

        return dg

    @property
    def num_nodes(self):
        return len(self.nodes)

    @property
    def nodes(self):
        """Get an iterable of pt nodes"""
        return self._directed_graph.values()

    @property
    def ordered_nodes(self):
        """Get an iterable of pt nodes in order of addition"""
        nodes = list(self.nodes)
        nodes.sort(key=lambda v: v._order)
        return nodes

    @property
    def exterior_cycles(self):
        """Computes all exterior boundaries.

        Returns:
            List of boundaries as list of nodes. The first polygon will
            be the outer exterior edge (in counter-clockwise order), and
            subsequent edges will be the edges of the holes in the graph
            (in clockwise order).
        """

        exterior_poly_lst = []
        exterior_check = {}

        for root_node in self.ordered_nodes:

            # Store node in check
            exterior_check[root_node.key] = None

            # Get next exterior adjacent node
            next_node = self.next_exterior_node(root_node)

            is_valid = (next_node is not None) and \
                (next_node.key not in exterior_check)

            if not is_valid:
                continue

            # Create list of exterior points
            exterior_poly = [root_node]
            # Add to dict to prevent repetition
            exterior_check[next_node.key] = None

            while next_node.key != root_node.key:
                exterior_poly.append(next_node)
                exterior_check[next_node.key] = None
                next_node = self.next_exterior_node(next_node)

            exterior_poly_lst.append(exterior_poly)

        return exterior_poly_lst

    def node(self, key):
        """Retrieves the node based on passed value.

        Args:
            val: The key for a node in the directed graph.

        Returns:
            The node for the passed key.
        """

        try:
            return self._directed_graph[key]
        except KeyError:
            return None

    def add_adj(self, node, adj_val_lst):
        """Adds nodes to node.adj_lst.

        This method will ensure no repetitions will occur in adj_lst.

        Args:
            node: _Node to add adjacencies to.
            adj_val_lst: List of Point2D objects to add as adjacent nodes.
        """
        adj_keys = {n.key: None for n in node.adj_lst}
        adj_keys[node.key] = None
        for adj_val in adj_val_lst:
            adj_key = _vector2hash(adj_val, self._tol)
            if adj_key in adj_keys:
                continue

            self._add_node(adj_key, adj_val, exterior=None)
            adj_keys[adj_key] = None
            node.adj_lst.append(self.node(adj_key))

    def remove_adj(self, node, adj_key_lst):
        """Removes nodes in node.adj_lst.

        Args:
            node: _Node to remove adjacencies to.
            adj_val_lst: List of adjacency keys to remove as adjacent nodes.
        """
        node.adj_lst = [n for n in node.adj_lst if n.key not in set(adj_key_lst)]

    def _add_node(self, key, val, exterior=None):
        """If key doesn't exist, add to dg.

        Helper function for add_node.
        """
        if key not in self._directed_graph:
            self._directed_graph[key] = _Node(key, val, self.num_nodes, [], exterior)
        return self._directed_graph[key]

    def add_node(self, val, adj_lst, exterior=None):
        """Consumes a polygon point, and computes its key value, and adds it in the
        graph if it doesn't exist. If it does exist it appends adj_lst to existing pt.

        Args:
            val: A Point2D object
            adj_lst: A list of Point2D objects

        Returns:
            The hashed key from the existing or new node.
        """

        # Get key
        key = _vector2hash(val, self._tol)

        # Get node if it exists
        self._add_node(key, val, exterior)

        node = self._directed_graph[key]

        # Add the adj_lst to dg, and leave exterior None
        self.add_adj(node, adj_lst)

        # If pass exterior boolean, change node attribute
        if exterior is not None:
            node.exterior = exterior

        return node.key

    def insert_node(self, node, new_val, next_node, exterior=None):
        """Insert node in the middle of an edge defined by node and next_node.

        Args:
            node: _Node to left.
            new_val: Value for middle node.
            next_node: _Node to right.
            exterior: Optional boolean for exterior attribute.

        Returns:
            key of new_val node.
        """
        # Add new_val as a node, with next_node as an adjacency
        new_key = self.add_node(new_val, [next_node.pt], exterior=exterior)

        # Update parent by adding new adjacency, and removing old adjacency
        self.add_adj(node, [self.node(new_key).pt])

        # Edge case where the new point is coincident to parent or next_point.
        # This occurs when intersection passes through a corner.
        if (new_key == next_node.key) or (new_key == node.key):
            return new_key

        self.remove_adj(node, [next_node.key])

        return new_key

    def node_exists(self, key):
        """True if node in directed graph else False"""
        return key in self._directed_graph

    def pt_exists(self, pt):
        """True if a point (as Point2D) in directed graph exists as node else False.
        """
        return self.node_exists(_vector2hash(pt, self._tol))

    def polygon_exists(self, polygon):
        """Check if polygon is in directed graph.

        Args:
            polygons: A Polygon2D object.
            dg: A PolygonDirectedGraph.

        Return:
            True if exists, else False.
        """
        vertices_loop = list(polygon.vertices)
        vertices_loop = vertices_loop + [vertices_loop[0]]

        for i in range(len(vertices_loop) - 1):
            pt1 = vertices_loop[i]
            pt2 = vertices_loop[i + 1]

            if not self.pt_exists(pt1):
                return False

            node1 = self.node(_vector2hash(pt1, self._tol))
            node2 = self.node(_vector2hash(pt2, self._tol))
            if node2.key in [n.key for n in node1.adj_lst]:
                return False

        return True

    def adj_matrix(self):
        """Gets an adjacency matrix of the directed graph where:

        * 1 = adjacency from row node to col node.
        * 0 = no adjacency.

        Returns:
            N x N square matrix where N is number of nodes.
        """
        nodes = self.ordered_nodes

        # Initialize amtx with no adjacencies
        amtx = [[0 for i in range(self.num_nodes)]
                for j in range(self.num_nodes)]

        for i in range(self.num_nodes):
            adj_indices = [adj._order for adj in nodes[i].adj_lst]
            for adj_idx in adj_indices:
                amtx[i][adj_idx] = 1

        return amtx

    def adj_matrix_labels(self):
        """Returns a dictionary where label key corresponds to index in adj_matrix
        and value is node key"""
        return {i: node.key for i, node in enumerate(self.ordered_nodes)}

    def intersect_graph_with_segment(self, segment):
        """Update graph with intersection of partial segment that crosses through polygon.

        Args:
            segment: LineSegment2D to intersect. Does not need to be contained within
            polygon.
        """
        int_key_lst = []

        for node in self.ordered_nodes:

            # Convert graph edge to trimming segment
            next_node = node.adj_lst[0]
            trim_seg = LineSegment2D.from_end_points(node.pt, next_node.pt)
            int_pt = intersection2d.intersect_line2d_infinite(trim_seg, segment)

            # Add intersection point as new node in graph
            if int_pt:
                int_key = self.insert_node(
                    node, int_pt, next_node, exterior=False)
                int_key_lst.append(int_key)

        # Add intersection edges
        if len(int_key_lst) == 2:
            # Typical case with convex cases
            # Make edge between intersection nodes
            n1, n2 = self.node(int_key_lst[0]), self.node(int_key_lst[1])
            self.add_node(n1.pt, [n2.pt], exterior=False)
            self.add_node(n2.pt, [n1.pt], exterior=False)

        elif len(int_key_lst) > 2:
            # Edge case with concave geometry creates multiple intersections
            # Sort distance and add adjacency
            n = self.node(int_key_lst[0])
            distances = [(0, 0.0)]

            for i, k in enumerate(int_key_lst[1:]):
                distance = LineSegment2D.from_end_points(n.pt, self.node(k).pt).length
                distances.append((i + 1, distance))

            distances = sorted(distances, key=lambda t: t[1])

            for i in range(len(distances)-1):
                k1, k2 = distances[i][0], distances[i+1][0]
                n1, n2 = self.node(int_key_lst[k1]), self.node(int_key_lst[k2])

                # Add bidirection so the min cycle works
                self.add_node(n1.pt, [n2.pt], exterior=False)
                self.add_node(n2.pt, [n1.pt], exterior=False)

    @staticmethod
    def is_edge_bidirect(node1, node2):
        """Are two nodes bidirectional.

        Args:
            node1: _Node object
            node2: _Node object

        Returns:
            True if node1 and node2 are in each other's adjacency list,
            else False.
        """
        return node1.key in (n.key for n in node2.adj_lst) and \
            node2.key in (n.key for n in node1.adj_lst)

    @staticmethod
    def next_unidirect_node(node):
        """Retrieves the first unidirectional point adjacent
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
            if not PolygonDirectedGraph.is_edge_bidirect(node, _next_node):
                next_node = _next_node
                break

        return next_node

    @staticmethod
    def next_exterior_node(node):
        """Retrieves the first exterior node adjacent
        to consumed node. They define an exterior or naked edge.

        Args:
            node: _Node

        Returns:
            Next node that defines exterior edge, or None if all
            adjacencies are bidirectional.
        """

        # Check bidirectionality
        next_node = None
        for _next_node in node.adj_lst:

            if _next_node.exterior is None:
                # If user-assigned attribute isn't defined, check bidirectionality.
                if not PolygonDirectedGraph.is_edge_bidirect(node, _next_node):
                    next_node = _next_node
                    break
            elif _next_node.exterior is True:
                next_node = _next_node
                break

        return next_node

    @staticmethod
    def exterior_cycle(cycle_root):
        """Computes exterior boundary from a given node.

        This method assumes that exterior edges are naked (unidirectional) and
        interior edges are bidirectional.

        Args:
            cycle_root: Starting _Node in exterior cycle.
        Returns:
            List of nodes on exterior if a cycle exists, else None.
        """

        # Get the first exterior edge
        curr_node = cycle_root
        next_node = PolygonDirectedGraph.next_exterior_node(curr_node)
        if not next_node:
            return None

        ext_cycle = [curr_node]
        while next_node.key != cycle_root.key:
            ext_cycle.append(next_node)
            next_node = PolygonDirectedGraph.next_exterior_node(next_node)
            if not next_node:
                return None

        return ext_cycle

    @staticmethod
    def _min_ccw_cycle(curr_node, next_node, cycle):
        """Identify the counter-clockwise adjacent node with the minimum angle.

        Args:
            curr_node: A node that defines first point of a counter-clockwise polygon edge.
            next_node: A node connected to the curr_node that defines the second point of a
                counter-clockwise polygon edge.
            cycle: Current list of nodes that will form a polygon.
        Returns:
            The next connected node that contains the minimum counter-clockwise angle
            by the edge defined by the curr_node and next_node.
        """
        # Get current edge direction vector
        # N.B point subtraction or addition results in Vector2D
        edge_dir = next_node.pt - curr_node.pt

        # Initialize values for comparison
        min_theta = float("inf")
        min_node = None

        # Identify the node with the smallest ccw angle
        for adj_node in next_node.adj_lst:
            # Make sure this node isn't backtracking by checking
            # new node isn't parent of next_node
            if adj_node.key == cycle[-2].key:
                continue

            # Get next edge direction vector
            next_edge_dir = adj_node.pt - next_node.pt
            theta = edge_dir.angle_clockwise(next_edge_dir * -1)

            if theta < min_theta:
                min_theta = theta
                min_node = adj_node

        return min_node

    @staticmethod
    def min_ccw_cycle(curr_node, next_node, recurse_limit=3000, print_recurse=False):
        """Recursively identify most counter-clockwise adjacent node and get closed loop.

        Args:
            curr_node: The first node, for first edge.
            next_node: The node next to ref_node that constitutes a polygon edge.
            recurse_limit: optional parameter to limit the number of while loop cycles.
                Default: 3000.
            print_recurse: optional boolean to print cycle loop cycles, for debugging.
                Default: False.

        Returns:
            A list of nodes that form a polygon if the cycle exists, else None.
        """
        # Set parameters
        count = 0
        cycle = [curr_node, next_node]
        while next_node.key != cycle[0].key:
            if print_recurse:
                print('"recursion" cycle: {}'.format(count))
            # Checks to ensure not trapped in while loop by degenerate skeleton.
            if recurse_limit and count >= recurse_limit:
                # Base case 1: recursion limit is hit
                raise RuntimeError('Error finding the minimum counterclockwise cycle '
                                   'in this polygon. Recursion limit of {} is '
                                   'exceeded.'.format(recurse_limit))
            elif next_node is None:
                # Base case 2: No node exists in adjacency list
                raise RuntimeError('Error finding the minimum counterclockwise cycle '
                                   'in this polygon')

            # Calculate minimum counter-clockwise cycle
            min_ccw_node = PolygonDirectedGraph._min_ccw_cycle(
                curr_node, next_node, cycle)

            # Update parameter names
            curr_node = next_node
            next_node = min_ccw_node
            cycle.append(next_node)
            count += 1

        # Remove final node which is duplicate of first node
        cycle.pop(0)

        return cycle
