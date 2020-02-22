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
from ladybug_geometry.geometry2d.pointvector import Point2D
from ladybug_geometry.geometry2d.line import LineSegment2D
from ladybug_geometry.geometry2d.ray import Ray2D
from ladybug_geometry import intersection2d


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
        return len(adj_lst)

class PolygonDirectedGraph(object):
    """Private dictionary that coordinates polygon pt adjacencies."""
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
        """ Retreives exterior boundary, by identifying unidirectional (naked) edges.

        Args:
            cycle_root: Starting _Node in exterior cycle.
        Returns:
            List of nodes on exterior.
        """

        gap_error = 'There is a no continuous outer boundary cycle from this node.'

        # Get the first exterior edge
        curr_node = cycle_root
        next_node = self.get_next_unidirect_node(curr_node)

        if not next_node:
            raise Exception(gap_error)

        ext_cycle = [curr_node]
        while next_node.key != cycle_root.key:
            ext_cycle.append(next_node)
            next_node = self.get_next_unidirect_node(next_node)
            if not next_node:
                raise Exception(gap_error)
                break

        return ext_cycle

    def find_most_ccw_cycle(self):
        """ Traverses a directed graph of connected edges and produces a list of
        individual polygons they define.

        Returns:
            A list of polygon point arrays.
        """

        polygon_lst = []

        nodes = self.ordered_nodes()
        for node in nodes:

            # Confirm if node on exterior boundary
            if not node._exterior:
                continue

            # Identify the next adjacent node on exterior.
            for i in range(node.adj_num):
                _next_node = node.adj_lst[i]
                if _next_node.exterior:
                    next_node = neighbor
                    break

            #Now we recursively check most ccw
            next_node_adj_lst = next_node.adj_lst
            cycle = [node]
            try:
                cycle = self.recurse_ccw(node, next_node, n_adj_lst, cycle, 0)
            except:
                pass
            # print '-------\n-----FINISHED CYCLE\n', cycle, '---\---\n'
            polygon_lst.append(cycle)

        return polygon_lst

    # def recurse_ccw(self,refn,nextn,lok,cycle,count):
    #     def get_ccw_angle(prev_dir,next_dir):
    #         #Input prev_dir vector and next_dir vector in CCW ordering
    #         #Output CCW angle between them in radians

    #         #Reverse prev_dir order for angle checking
    #         #We create a new vector b/c must be ccw order for reflex check
    #         reverse_prev_dir = prev_dir * -1.0
    #         #Use the dot product to find the angle
    #         dotprod = rc.Geometry.Vector3d.Multiply(reverse_prev_dir,next_dir)
    #         try:
    #             cos_angle = dotprod/(prev_dir.Length * next_dir.Length)
    #         except ZeroDivisionError:
    #             print 'ZeroDivisionError'
    #             cos_angle = 0.0

    #         # Get angle from dot product
    #         # This will be between 0 and pi
    #         # b/c -1 < cos theta < 1
    #         dotrad = math.acos(cos_angle)

    #         #Use 2d cross product (axby - bxay) to see if next_vector is right/left
    #         #This requires ccw ordering of vectors
    #         #If cross is positive (for ccw ordering) then next_vector is to left (inner)
    #         #If cross is negative (for ccw ordering) then next_vector is to right (outer)
    #         #If cross is equal then zero vector, then vectors are colinear. Assume inner.

    #         cross_z_sign = prev_dir[0] * next_dir[1] - prev_dir[1] * next_dir[0]
    #         #print 'deg: ', round(math.degrees(dotrad),2)
    #         #If reflex we must subtract 2pi from it to get reflex angle
    #         if cross_z_sign < 0.0:
    #             dotrad = 2*math.pi - dotrad
    #         return dotrad

    #     #Input list of keys
    #     #output key with min ccw
    #     cycle.append(nextn)
    #     cycle_id_lst = map(lambda n: n.id, cycle)

    #     # Check the base case where loop is completed
    #     if nextn.id == cycle_id_lst[0] or count > 20:
    #         return cycle

    #     #print 'cycle', cycle_id_lst

    #     #reference direction vector
    #     ref_edge_dir =  nextn.value - refn.value

    #     min_rad = float("Inf")
    #     min_node = None

    #     # Identify the node with the smallest ccw angle
    #     for i in xrange(len(lok)):
    #         k = lok[i]
    #         n2chk = self.adj_graph[k]
    #         #Make sure we don't backtrack
    #         if n2chk.id == cycle_id_lst[-2]:
    #             continue
    #         chk_edge_dir = n2chk.value - nextn.value
    #         #print 'chkccw', refn.id, '--', nextn.id, '--->', n2chk.id
    #         rad = get_ccw_angle(ref_edge_dir,chk_edge_dir)
    #         if rad < min_rad:
    #             min_rad = rad
    #             min_node = n2chk
    #         #print '---'
    #     #print 'min is', n2chk.id,':', round(math.degrees(rad),2)
    #     #print '---'
    #     alok = min_node.adj_lst

    #     return self.recurse_ccw(nextn, min_node, alok, cycle, count+1)

# class AdjGraph(object):
#     """
#     # Graph as adjacency list
#     #Good ref for adj graphs: http://interactivepython.org/runestone/static/pythonds/Graphs/Implementation.html
#     #adj_graph is a dict like this:
#     #{key1: _AdjGraphNode.adj_lst = [2,4],
#     # key2: _AdjGraphNode.adj_lst = [3,4],
#     # key3: _AdjGraphNode.adj_lst = [4],
#     # key4: _AdjGraphNode.adj_lst = [1]}

#     Graph
#     1     4
#        /  |
#     2     3

#     Adj Matrix:
#     [[x 1 2 3 4],
#      [1 - - - -],
#      [2 - - - x],
#      [3 - - - x],
#      [4 - - - -]]
#     """

#     def __init__(self,adj_graph=None):
#         self.adj_graph = adj_graph if adj_graph != None else {}
#         self.num_node = len(self.adj_graph.keys())
#     def vector2hash(self,vector,tol=4):
#         #Tolerance set to
#         myhash = "("
#         for i in xrange(len(vector)):
#             coordinate = vector[i]
#             myhash += str(round(coordinate,tol))
#             if i < len(vector)-1:
#                 myhash += ","
#         myhash += ")"
#         return myhash
#     def add_node(self,key,value,is_out_edge=False):
#         #_AdjGraphNode is a private class
#         #Instantiate _AdjGraphNode, we creates key = num_node
#         if key in self.adj_graph:
#             n = self.adj_graph[key]
#             print n.id, ' key already exists in adj_graph!'
#             return self.adj_graph[key]
#         id = len(self.adj_graph.keys())
#         adj_graph_node = _AdjGraphNode(key,value,id,is_out_edge)
#         #Now add it to the adj_graph
#         self.adj_graph[key] = adj_graph_node
#         self.num_node += 1
#         return adj_graph_node
#     def __getitem__(self,k):
#         if k in self.adj_graph:
#             return self.adj_graph[k]
#         else:
#             return None
#     def keylst_2_nodelst(self,keylst):
#         return map(lambda k: self.adj_graph[k],keylst)
#     def add_directed_edge(self,key,key2add):
#         #This function will add existing node key to adjacency list of
#         #another node indicating a directed edge
#         if key in self.adj_graph and key2add in self.adj_graph:
#             node = self.adj_graph[key]
#             if key2add in node.adj_lst or key2add == key:
#                 print 'key2add already in adj list or self-intersection'
#                 return None
#             node.adj_lst.append(key2add)
#         else:
#             print 'key not in adj graph'
#     def recurse_ccw(self,refn,nextn,lok,cycle,count):
#         def get_ccw_angle(prev_dir,next_dir):
#             #Input prev_dir vector and next_dir vector in CCW ordering
#             #Output CCW angle between them in radians

#             #Reverse prev_dir order for angle checking
#             #We create a new vector b/c must be ccw order for reflex check
#             reverse_prev_dir = prev_dir * -1.0
#             #Use the dot product to find the angle
#             dotprod = rc.Geometry.Vector3d.Multiply(reverse_prev_dir,next_dir)
#             try:
#                 cos_angle = dotprod/(prev_dir.Length * next_dir.Length)
#             except ZeroDivisionError:
#                 print 'ZeroDivisionError'
#                 cos_angle = 0.0

#             # Get angle from dot product
#             # This will be between 0 and pi
#             # b/c -1 < cos theta < 1
#             dotrad = math.acos(cos_angle)

#             #Use 2d cross product (axby - bxay) to see if next_vector is right/left
#             #This requires ccw ordering of vectors
#             #If cross is positive (for ccw ordering) then next_vector is to left (inner)
#             #If cross is negative (for ccw ordering) then next_vector is to right (outer)
#             #If cross is equal then zero vector, then vectors are colinear. Assume inner.

#             cross_z_sign = prev_dir[0] * next_dir[1] - prev_dir[1] * next_dir[0]
#             #print 'deg: ', round(math.degrees(dotrad),2)
#             #If reflex we must subtract 2pi from it to get reflex angle
#             if cross_z_sign < 0.0:
#                 dotrad = 2*math.pi - dotrad
#             return dotrad

#         if True: pass #weird code folding glitch neccessitatest this
#         #Input list of keys
#         #output key with most ccw
#         #Base case
#         #print 'startid, chkid:', cycle[0].id, nextn.id
#         cycle.append(nextn)
#         cycle_id_lst = map(lambda n: n.id, cycle)
#         if nextn.id == cycle_id_lst[0] or count > 20:
#             return cycle

#         #print 'cycle', cycle_id_lst

#         #reference direction vector
#         ref_edge_dir =  nextn.value - refn.value

#         min_rad = float("Inf")
#         min_node = None

#         for i in xrange(len(lok)):
#             k = lok[i]
#             n2chk = self.adj_graph[k]
#             #Make sure we don't backtrack
#             if n2chk.id == cycle_id_lst[-2]:
#                 continue
#             chk_edge_dir = n2chk.value - nextn.value
#             #print 'chkccw', refn.id, '--', nextn.id, '--->', n2chk.id
#             rad = get_ccw_angle(ref_edge_dir,chk_edge_dir)
#             if rad < min_rad:
#                 min_rad = rad
#                 min_node = n2chk
#             #print '---'
#         #print 'min is', n2chk.id,':', round(math.degrees(rad),2)
#         #print '---'
#         alok = min_node.adj_lst

#         return self.recurse_ccw(nextn,min_node,alok,cycle,count+1)
#     def find_most_ccw_cycle(self):
#         #def helper_most_ccw(lok):

#         #Input adjacency graph
#         #Output loc: listof (listof (listof pts in closed cycle))
#         LOC = []
#         keylst = self.get_sorted_keylst()
#         for i in xrange(len(keylst)):
#             key = keylst[i]
#             root_node = self.adj_graph[key]
#             if not root_node.is_out_edge:
#                 continue

#             #Identify the next node on outer edge
#             #b/c outer edge vertexes are placed first in adj graph
#             #worst complexity <= O(n)
#             for i in xrange(root_node.num_neighbor()):
#                 adj_key = root_node.adj_lst[i]
#                 neighbor = self.adj_graph[adj_key]
#                 if neighbor.is_out_edge:
#                     next_node = neighbor
#                     break

#             #Now we recursively check most ccw
#             n_adj_lst = next_node.adj_lst
#             cycle = [root_node]
#             try:
#                 cycle = self.recurse_ccw(root_node,next_node,n_adj_lst,cycle,0)
#             except:
#                 pass
#             #print '-------\n-----FINISHED CYCLE\n', cycle, '---\---\n'
#             LOC.append(cycle)
#         #print '-'
#         return LOC
#     def get_sorted_keylst(self):
#         valuelst = self.adj_graph.values()
#         valuelst.sort(key=lambda v: v.id)
#         keylst = map(lambda v: v.key,valuelst)
#         return keylst
#     def is_near_zero(self,num,eps=1E-10):
#         return abs(float(num)) < eps
#     def __repr__(self):
#         keylst = self.get_sorted_keylst()
#         strgraph = ""
#         for i in xrange(len(keylst)):
#             key = keylst[i]
#             node = self.adj_graph[key]
#             strgraph += str(node.id) + ': '
#             strgraph += str(map(lambda k: self.adj_graph[k].id,node.adj_lst))
#             strgraph += '\n'
#         return strgraph
#     def __contains__(self,key):
#         return self.adj_graph.has_key(key)