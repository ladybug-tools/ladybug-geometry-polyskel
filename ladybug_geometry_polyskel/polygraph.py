# coding=utf-8
"""Implementation of a Directed Graph data structure."""
from __future__ import division
import math

from ladybug_geometry.geometry2d import LineSegment2D, Polygon2D
from ladybug_geometry.intersection2d import intersect_line2d_infinite

from .polyskel import skeleton_as_edge_list, _intersect_skeleton_segments, \
    _remove_segments_outside_boundary


def _vector2hash(vector, tol):
    """Hashes spatial coordinates for use in dictionary.

    Args:
        vector: A Vector2D object.
        tol: floating point precision tolerance.

    Returns:
        Hash of vector as a string of rounded coordinates.
    """
    try:
        rtol = (int(math.log10(tol)) * -1)
    except ValueError:
        rtol = 0
    # avoid cases of signed zeros messing with keys
    x_val = 0.0 if vector.x == 0 else vector.x
    y_val = 0.0 if vector.y == 0 else vector.y
    return str((round(x_val, rtol), round(y_val, rtol)))


class _Node(object):
    """Private class to handle nodes in PolygonDirectedGraph.

    Args:
        val: A Point2D object.
        key: Hash of Point2D object.
        order: integer counting order of the Node (based on dg propagation)
        adj_lst: list of keys adjacent to this node.
        exterior: Node boundary condition. None if not set by user, else True
            or False according to user.

    Properties:
        * pt
        * key
        * adj_lst
        * exterior
        * adj_count
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

    @property
    def adj_count(self):
        """Number of adjacent nodes"""
        return len(self.adj_lst)

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return isinstance(other, _Node) and \
            self.key == other.key

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return '{}: {}'.format(self._order, self.key)


class PolygonDirectedGraph(object):
    """A directed graph data structure for point relationships.

    The PolygonDirectedGraph effectively represents the network of a straight
    skeleton and assists with finding the shortest pathways through it. It also
    helps differentiate interior from exterior parts of the graph. Typically,
    interior pathways are bi-directional in the graph while exterior pathways
    are uni-directional.

    Args:
        tolerance: Tolerance for point equivalence. (Default: 1e-5). This is used
            for hashing points within the network.

    Properties:
        * node_count: Integer for the number of nodes in graph.
        * nodes: An iterable of nodes in graph.
        * ordered_nodes: An iterable of nodes in graph in order they were added.
        * connection_segments: List of LineSegment2D for the node connections
        * outer_root_node: A node for the outer root key
        * hole_root_nodes: A list of nodes for the hole root keys
        * is_intersect_topology: A boolean for whether the graph self-intersects
    """

    def __init__(self, tol):
        """Initialize a PolygonDirectedGraph."""
        self._directed_graph = {}
        self._tol = tol
        self.outer_root_key = None  # will be set during skeleton creation
        self.hole_root_keys = []  # will be set during skeleton creation
        self.is_intersect_topology = None  # will be set during skeleton creation

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
            dg.add_node(point_array[i], [point_array[i + 1]], exterior=None)
        if loop:
            dg.add_node(point_array[-1], [point_array[0]], exterior=None)
        return dg

    @property
    def node_count(self):
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
    def outer_root_node(self):
        """Get the node of the outer boundary root."""
        return self.node(self.outer_root_key)

    @property
    def hole_root_nodes(self):
        """Get a list of nodes for the roots of the holes."""
        return [self.node(hole_key) for hole_key in self.hole_root_keys]

    @property
    def connection_segments(self):
        """Get a list of LineSegment2D for the node connections in the graph."""
        traversed = set()
        connections = []
        for node in self.nodes:
            for conn_node in node.adj_lst:
                if (conn_node.key, node.key) not in traversed:
                    conn_seg = LineSegment2D.from_end_points(node.pt, conn_node.pt)
                    connections.append(conn_seg)
                    traversed.add((node.key, conn_node.key))
        return connections

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
            return None  # broken connection

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

    def add_node(self, val, adj_lst, exterior=None):
        """Add a node into the PolygonDirectedGraph.

        This method consumes a Point2D, computes its key value, and adds it in the
        graph if it doesn't exist. If it does exist it appends adj_lst to existing pt.

        Args:
            val: A Point2D object.
            adj_lst: A list of Point2D objects adjacent to the node.
            exterior: Optional boolean for exterior attribute.

        Returns:
            The hashed key from the existing or new node.
        """
        key = _vector2hash(val, self._tol)  # get key
        self._add_node(key, val, exterior)  # get node if it exists
        node = self._directed_graph[key]
        self.add_adj(node, adj_lst)  # add the adj_lst to dg
        # if the exterior boolean was passed, change the node attribute
        if exterior is not None:
            node.exterior = exterior
        return node.key

    def _add_node(self, key, val, exterior=None):
        """Helper function for add_node. If key doesn't exist, add to dg."""
        if key not in self._directed_graph:
            self._directed_graph[key] = _Node(key, val, self.node_count, [], exterior)
        return self._directed_graph[key]

    def insert_node(self, base_node, new_val, next_node, exterior=None):
        """Insert node in the middle of an edge defined by node and next_node.

        Args:
            base_node: _Node object to the left.
            new_val:  A Point2D object for the new node in the middle.
            next_node: _Node object to the right.
            exterior: Optional boolean for exterior attribute.

        Returns:
            key of new_val node.
        """
        # add new_val as a node, with next_node as an adjacency
        new_key = self.add_node(new_val, [next_node.pt], exterior=exterior)
        # update parent by adding new adjacency, and removing old adjacency
        self.add_adj(base_node, [self.node(new_key).pt])

        # catch the edge case where the new point is coincident to parent or next_point.
        # this occurs when intersection passes through a corner.
        if (new_key == next_node.key) or (new_key == base_node.key):
            return new_key
        self.remove_adj(base_node, [next_node.key])
        return new_key

    def node_exists(self, key):
        """Check if a node is in the graph. True if node in directed graph else False."""
        return key in self._directed_graph

    def pt_exists(self, pt):
        """True if a point (as Point2D) in directed graph exists as node else False.
        """
        return self.node_exists(_vector2hash(pt, self._tol))

    def polygon_exists(self, polygon):
        """Check if polygon is in the directed graph.

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
        # initialize a mtx with no adjacencies
        amtx = [[0 for i in range(self.node_count)]
                for j in range(self.node_count)]

        for i in range(self.node_count):
            adj_indices = [adj._order for adj in nodes[i].adj_lst]
            for adj_idx in adj_indices:
                amtx[i][adj_idx] = 1

        return amtx

    def adj_matrix_labels(self):
        """Returns a dictionary where label key corresponds to index in adj_matrix
        and value is node key"""
        return {i: node.key for i, node in enumerate(self.ordered_nodes)}

    def intersect_graph_with_segment(self, segment):
        """Update graph with intersection of partial segment crossing through polygon.

        Args:
            segment: LineSegment2D to intersect with the graph. The LineSegment2D
                does not need to be contained within polygon and will be intersected
                infinitely.
        """
        # loop through all nodes in the graph and find intersection points
        int_key_lst = []
        for node in self.ordered_nodes:
            # convert graph edge to trimming segment
            next_node = node.adj_lst[0]
            trim_seg = LineSegment2D.from_end_points(node.pt, next_node.pt)
            int_pt = intersect_line2d_infinite(trim_seg, segment)

            # add intersection point as new node in graph
            if int_pt:
                int_key = self.insert_node(node, int_pt, next_node, exterior=False)
                int_key_lst.append(int_key)

        # add intersection edges between the newly-found nodes
        if len(int_key_lst) == 2:
            # typical case with convex cases
            # make edge between intersection nodes
            n1, n2 = self.node(int_key_lst[0]), self.node(int_key_lst[1])
            self.add_node(n1.pt, [n2.pt], exterior=False)
            self.add_node(n2.pt, [n1.pt], exterior=False)

        elif len(int_key_lst) > 2:
            # edge case with concave geometry creates multiple intersections
            # sort distance and add adjacency
            n = self.node(int_key_lst[0])
            distances = [(0, 0.0)]

            for i, k in enumerate(int_key_lst[1:]):
                distance = LineSegment2D.from_end_points(n.pt, self.node(k).pt).length
                distances.append((i + 1, distance))

            distances = sorted(distances, key=lambda t: t[1])

            for i in range(len(distances) - 1):
                k1, k2 = distances[i][0], distances[i + 1][0]
                n1, n2 = self.node(int_key_lst[k1]), self.node(int_key_lst[k2])
                # add bi-direction so the min cycle works
                self.add_node(n1.pt, [n2.pt], exterior=False)
                self.add_node(n2.pt, [n1.pt], exterior=False)

    def min_cycle(self, base_node, goal_node, ccw_only=False):
        """Identify the shortest interior cycle between two exterior nodes.

        Args:
            base_node: The first exterior node of the edge.
            goal_node: The end exterior node of the cycle that, together with
                the base_node, constitutes an exterior edge.
            ccw_only: A boolean to note whether the search should be limited
                to the counter-clockwise direction only. (Default: False).

        Returns:
            A list of nodes that form a polygon if the cycle exists, else None.
        """
        # set up a queue for exploring the graph
        explored = []
        queue = [[base_node]]
        orig_dir = base_node.pt - goal_node.pt  # yields a vector
        # loop to traverse the graph  with the help of the queue
        while queue:
            path = queue.pop(0)
            node = path[-1]
            # make sure that the current node has not been visited
            if node not in explored:
                prev_dir = node.pt - path[-2].pt if len(path) > 1 else orig_dir
                # iterate over the neighbors to determine relevant nodes
                rel_neighbors, rel_angles = [], []
                for neighbor in node.adj_lst:
                    if neighbor == goal_node:  # the shortest path was found!
                        path.append(goal_node)
                        return path
                    elif neighbor.exterior:
                        continue  # don't traverse the graph exterior
                    edge_dir = neighbor.pt - node.pt
                    cw_angle = prev_dir.angle_clockwise(edge_dir * -1)
                    if not (1e-5 < cw_angle < (2 * math.pi) - 1e-5):
                        continue  # prevent back-tracking along the search
                    rel_neighbors.append(neighbor)
                    rel_angles.append(cw_angle)
                # sort the neighbors by clockwise angle
                if len(rel_neighbors) > 1:
                    rel_neighbors = [n for _, n in sorted(zip(rel_angles, rel_neighbors),
                                                          key=lambda pair: pair[0])]
                # add the relevant neighbors to the path and the queue
                if ccw_only:
                    new_path = list(path)
                    new_path.append(rel_neighbors[0])
                    queue.append(new_path)
                else:  # add all neighbors to the search
                    for neighbor in rel_neighbors:
                        new_path = list(path)
                        new_path.append(neighbor)
                        queue.append(new_path)
                explored.append(node)
        # if we reached the end of the queue, then no path was found
        return None

    def exterior_cycle(self, cycle_root):
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

        # loop through the cycle until we get it all or run out of points
        max_iter = self.node_count + 1  # maximum length a cycle can be
        ext_cycle = [curr_node]
        iter_count = 0
        while next_node.key != cycle_root.key:
            ext_cycle.append(next_node)
            next_node = PolygonDirectedGraph.next_exterior_node(next_node)
            if not next_node:
                return None  # we have hit a dead end in the cycle
            iter_count += 1
            if iter_count > max_iter:
                break  # we have gotten stuck in a loop

        return ext_cycle

    def exterior_cycles(self):
        """Get a list of lists where each sub-list is an exterior cycle of Nodes."""
        exterior_poly_lst = []  # list to store cycles
        explored_nodes = set()  # set to note explored exterior nodes
        max_iter = self.node_count + 1  # maximum length a cycle can be

        # loop through all of the nodes of the graph and find cycles
        for root_node in self.ordered_nodes:
            # make a note that the current node has been explored
            explored_nodes.add(root_node.key)  # mark the node as explored
            # get next exterior adjacent node and check that it's valid
            next_node = self.next_exterior_node(root_node)  # mark the node as explored
            is_valid = (next_node is not None) and (next_node.key not in explored_nodes)
            if not is_valid:
                continue
            # make a note that the next node has been explored
            explored_nodes.add(next_node.key)

            # traverse the loop of points until we get back to start or hit a dead end
            exterior_poly = [root_node]
            prev_node = root_node
            iter_count = 0
            while next_node.key != root_node.key:
                exterior_poly.append(next_node)
                explored_nodes.add(next_node.key)  # mark the node as explored
                follow_node = self.next_exterior_node_no_backtrack(
                    next_node, prev_node, explored_nodes)
                prev_node = next_node  # set as the previous node for the next step
                next_node = follow_node
                if next_node is None:
                    break  # we have hit a dead end in the cycle
                iter_count += 1
                if iter_count > max_iter:
                    print('Extraction of core polygons hit an endless loop.')
                    break  # we have gotten stuck in a loop
            exterior_poly_lst.append(exterior_poly)

        # return all of the exterior loops that were found
        return exterior_poly_lst

    @staticmethod
    def next_exterior_node_no_backtrack(node, previous_node, explored_nodes):
        """Get the next exterior node adjacent to the input node.

        This method is similar to the next_exterior_node method but it includes
        extra checks to handle intersections with 3 or more segments in the
        graph exterior cycles. In these cases a set of previously explored_nodes
        is used to ensure that no back-tracking happens over the search of the
        network, which can lead to infinite looping through the graph. Furthermore,
        the previous_node is used to select the pathway with the smallest angle
        difference with the previous direction. This leads the result towards
        minimal polygons with fewer self-intersecting loops.

        Args:
            node: A _Node object for which the next node will be returned.
            previous_node: A _Node object for the node that came before
                the current one in the loop. This will be used in the event that
                multiple exterior nodes are found connecting to the input node.
                In this case, the exterior node with the smallest angle difference
                with the previous direction will be returned. This leads the
                result towards minimal polygons and away from self-intersecting
                exterior loops like a bowtie.

        Returns:
            Next node that defines exterior edge, or None if all adjacencies are
            bidirectional.
        """
        # loop through the all adjacent nodes and determine if they are exterior
        next_nodes = []
        for _next_node in node.adj_lst:
            if _next_node.exterior:  # user has labeled it as exterior; we're done!
                return _next_node
            elif _next_node.exterior is None:  # don't know if it's interior or exterior
                # if user-assigned attribute isn't defined, check bi-directionality
                if not PolygonDirectedGraph.is_edge_bidirect(node, _next_node):
                    next_nodes.append(_next_node)

        # evaluate whether there is one obvious choice for the next node
        if len(next_nodes) <= 1:
            return next_nodes[0] if len(next_nodes) == 1 else None
        next_nodes = [nn for nn in next_nodes if nn.key not in explored_nodes]
        if len(next_nodes) <= 1:
            return next_nodes[0] if len(next_nodes) == 1 else None

        # if we have multiple exterior nodes, use the previous node to find the best one
        prev_dir = previous_node.pt - node.pt  # yields a vector
        next_angles = []
        for next_node in next_nodes:
            edge_dir = next_node.pt - node.pt  # yields a vector
            next_angles.append(prev_dir.angle(edge_dir * -1))
        sorted_nodes = [n for _, n in sorted(zip(next_angles, next_nodes),
                                             key=lambda pair: pair[0])]
        return sorted_nodes[0]  # return the node making the smallest angle

    @staticmethod
    def next_exterior_node(node):
        """Get the next exterior node adjacent to consumed node.

        If there are adjacent nodes that are labeled as exterior, with True or
        False defining the _Node.exterior property, the first of such nodes in
        the adjacency list will be returned as the next one. Otherwise, the
        bi-directionality will be used to determine whether the next node is
        exterior.

        Args:
            node: A _Node object for which the next node will be returned.

        Returns:
            Next node that defines exterior edge, or None if all adjacencies are
            bidirectional.
        """
        # loop through the adjacency and find an exterior node
        for _next_node in node.adj_lst:
            if _next_node.exterior:  # user has labeled it as exterior; we're done!
                return _next_node
            elif _next_node.exterior is None:  # don't know if it's interior or exterior
                # if user-assigned attribute isn't defined, check bi-directionality
                if not PolygonDirectedGraph.is_edge_bidirect(node, _next_node):
                    return _next_node
        return None

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

    def __repr__(self):
        """Represent PolygonDirectedGraph."""
        s = ''
        for n in self.ordered_nodes:
            s += '{}, [{}]\n'.format(
                n.pt.to_array(),
                ', '.join([str(_n.pt.to_array()) for _n in n.adj_lst]))
        return s


def skeleton_as_directed_graph(boundary, holes=None, tolerance=1e-5):
    """Compute the straight skeleton of a shape as a PolygonDirectedGraph.

    Args:
        boundary: A ladybug-geometry Polygon2D for the boundary around the shape
            for which the straight skeleton will be computed.
        holes: An optional list of ladybug-geometry Polygon2D for the holes within
            the shape for which a straight skeleton will be computed. If None,
            it will be assumed that no holes exist in the shape. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        A PolygonDirectedGraph object that represents the network and boundary
        of the straight skeleton. All interior connections in the graph are
        assumed to be bi-directional. The edges (including the boundary
        and holes) are uni-directional with the outer boundary being counterclockwise
        and the holes being clockwise. In other words, the fill of the shape is
        always to the left of each exterior edge. The nodes at the boundary
        and the holes have the exterior property set to True.
    """
    # get the segments representing the straight skeleton
    skeleton = skeleton_as_edge_list(boundary, holes, tolerance)
    skeleton, is_intersect_topology = _intersect_skeleton_segments(skeleton, tolerance)
    skeleton = _remove_segments_outside_boundary(skeleton, boundary, tolerance)

    # ensure the boundary and holes are oriented correctly for the graph
    if boundary.is_clockwise:
        boundary = boundary.reverse()
    loops = [boundary.vertices]
    if holes is not None:
        for hole in holes:
            if hole.is_clockwise:
                loops.append(hole.vertices)
            else:
                loops.append(hole.reverse().vertices)

    # make the directed graph and add the nodes for the boundary + holes
    dg = PolygonDirectedGraph(tol=tolerance)
    for loop_count, vertices in enumerate(loops):
        for j in range(len(vertices) - 1):
            curr_v = vertices[j]
            next_v = vertices[j + 1]
            k = dg.add_node(curr_v, [next_v], exterior=True)
            if j == 0:
                if loop_count == 0:
                    dg.outer_root_key = k
                else:
                    dg.hole_root_keys.append(k)
        dg.add_node(vertices[-1], [vertices[0]], exterior=True)  # close loop

    # add the straight skelton to the graph
    for seg in skeleton:
        # add a bidirectional edge to represent skeleton edges
        dg.add_node(seg.p2, [seg.p1])
        dg.add_node(seg.p1, [seg.p2], exterior=False)

    # set the property to track whether the graph is self-intersecting
    dg.is_intersect_topology = is_intersect_topology
    return dg


def skeleton_as_cycle_polygons(boundary, holes=None, tolerance=1e-5):
    """Compute the straight skeleton of a shape as a list of cycle polygons.

    Args:
        boundary: A ladybug-geometry Polygon2D for the boundary around the shape
            for which the straight skeleton will be computed.
        holes: An optional list of ladybug-geometry Polygon2D for the holes within
            the shape for which a straight skeleton will be computed. If None,
            it will be assumed that no holes exist in the shape. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        A list of Polygon2D that represent the straight skeleton of the input shape.
        There will be one Polygon2D for each edge of the shape (including both
        the boundary and the holes). Together, the Polygon2Ds will fill the
        entire input shape. There may be some overlap in between the output
        polygons if the straight skeleton did not have the correct topology.
    """
    # get the straight skeleton as a PolygonDirectedGraph
    dg = skeleton_as_directed_graph(boundary, holes, tolerance)

    # function to add cycles to the list of polygons to be returned
    cycle_polys = []

    def _add_cycle_polygon(min_cycle):
        """Add a cycle of Nodes to the list of polygons."""
        if min_cycle is not None:
            cycle_poly = Polygon2D([node.pt for node in min_cycle])
            cycle_polys.append(cycle_poly)

    # convert the edge cycles to Polygon2D
    exter_cycle = dg.exterior_cycle(dg.outer_root_node)
    for i, base_node in enumerate(exter_cycle[:-1]):
        next_node = exter_cycle[i + 1]
        _add_cycle_polygon(dg.min_cycle(next_node, base_node))
    _add_cycle_polygon(dg.min_cycle(exter_cycle[0], exter_cycle[-1]))
    if holes is not None:
        for hole_root in dg.hole_root_nodes:
            exter_cycle = dg.exterior_cycle(hole_root)
            for i, base_node in enumerate(exter_cycle[:-1]):
                next_node = exter_cycle[i + 1]
                _add_cycle_polygon(dg.min_cycle(next_node, base_node))
            _add_cycle_polygon(dg.min_cycle(exter_cycle[0], exter_cycle[-1]))
    return cycle_polys
