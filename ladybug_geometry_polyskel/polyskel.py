# coding=utf-8
"""Implementation of the straight skeleton algorithm by Felkel and Obdrzalek[1].

The functions and classes here here are derived directly from the polyskel python
library by Armin Scipiades (@Bottfy), which is available at:
https://github.com/Botffy/polyskel

[1] Felkel, Petr and Stepan Obdrzalek. 1998. "Straight Skeleton Implementation." In
Proceedings of Spring Conference on Computer Graphics, Budmerice, Slovakia. 210 - 218.
"""
from __future__ import division

import heapq
from itertools import tee, islice, cycle, chain
from collections import namedtuple
import operator

# Geometry classes
from ladybug_geometry.geometry2d import Point2D, Ray2D, LineSegment2D, Polygon2D
from ladybug_geometry.intersection2d import intersect_line_segment2d, \
    intersect_line2d

# Polygon sorting classes
_OriginalEdge = namedtuple('_OriginalEdge', 'edge bisector_left, bisector_right')
Subtree = namedtuple('Subtree', 'source, height, sinks')
_SplitEventSubClass = namedtuple(
    '_SplitEvent', 'distance, intersection_point, vertex, opposite_edge')
_EdgeEventSubClass = namedtuple(
    '_EdgeEvent', 'distance intersection_point vertex_a vertex_b')


def _window(lst):
    """Window operator for lists.

    Consumes a list of items, and returns a zipped list of the previous,
    current and next items in the list, accessible by the same index.

    Args:
        lst: list

    Returns:
        Zipped list of previous, current and next items in list.
    """
    prevs, items, nexts = tee(lst, 3)
    prevs = islice(cycle(prevs), len(lst) - 1, None)
    nexts = islice(cycle(nexts), 1, None)
    return zip(prevs, items, nexts)


def _cross(a, b):
    """Get the determinant between two Vector2Ds and/or Point2Ds."""
    return a.x * b.y - b.x * a.y


def _approximately_equals(a, b):
    """Determine whether two Point2Ds or Vector2Ds are equal within a relative tolerance.
    """
    return a == b or (abs(a - b) <= max(abs(a), abs(b)) * 0.001)


def _normalize_contour(contour, tol):
    """Consumes list of x,y coordinate tuples and returns list of Point2Ds.

    Args:
        contour: list of x,y tuples from contour.
        tol: Number for point equivalence tolerance.

    Return:
         list of Point2Ds of contour.
    """
    contour = [Point2D(float(x), float(y)) for (x, y) in contour]
    normed_contour = []
    for prev, point, next in _window(contour):
        normed_prev = (point - prev).normalize()
        normed_next = (next - point).normalize()

        if not point.is_equivalent(next, tol) or \
                normed_prev.is_equivalent(normed_next, tol):
            normed_contour.append(point)

    return normed_contour


class _SplitEvent(_SplitEventSubClass):
    """A SplitEvent is a reflex vertex that splits an Edge Event.

    They therefore split the entire polygon and create new adjacencies between
    the split edge and each of the two edges incident to the reflex
    vertex (Felkel and Obdrzalek 1998, 1).
    """
    __slots__ = ()

    def __lt__(self, other):
        return self.distance < other.distance

    def __str__(self):
        return "{} Split event @ {} from {} to {}".format(
            self.distance, self.intersection_point, self.vertex, self.opposite_edge)


class _EdgeEvent(_EdgeEventSubClass):
    """An EdgeEvent is an edge extended from a perimeter edge, that shrinks to zero.

    This will make its neighboring edges adjacent (Felkel and Obdrzalek 1998, 2).
    """
    __slots__ = ()

    def __lt__(self, other):
        return self.distance < other.distance

    def __str__(self):
        return "{} Edge event @ {} between {} and {}".format(
            self.distance, self.intersection_point, self.vertex_a, self.vertex_b)


class _LAVertex:

    def __init__(self, point, edge_left, edge_right, direction_vectors=None, tol=1e-5):
        """A vertex in a double connected circular list of active vertices."""
        self.tol = tol
        self.point = point
        self.edge_left = edge_left
        self.edge_right = edge_right
        self.prev = None
        self.next = None
        self.lav = None
        # TODO this might be handled better. Maybe membership in lav implies validity?
        self._valid = True

        creator_vectors = (edge_left.v.normalize() * -1, edge_right.v.normalize())
        if direction_vectors is None:
            direction_vectors = creator_vectors

        self._is_reflex = (_cross(*direction_vectors)) < 0
        self._bisector = Ray2D(
            self.point, operator.add(*creator_vectors) * (-1 if self.is_reflex else 1))

    @property
    def bisector(self):
        return self._bisector

    @property
    def is_reflex(self):
        return self._is_reflex

    @property
    def original_edges(self):
        return self.lav._slav._original_edges

    def next_event(self):
        events = []
        if self.is_reflex:
            # A reflex vertex may generate a split event.
            # Split events happen when a vertex hits an opposite edge,
            # thereby splitting the polygon in two.
            for edge in self.original_edges:
                if edge.edge == self.edge_left or edge.edge == self.edge_right:
                    continue

                # A potential b is at the intersection between our own bisector and
                # the bisector of the angle between the tested edge and any one of our
                # own edges.
                # We choose 'less parallel' edge (to exclude a potentially parallel edge)

                # Make normalized copies of vectors
                norm_edge_left_v = self.edge_left.v.normalize()
                norm_edge_right_v = self.edge_right.v.normalize()
                norm_edge_v = edge.edge.v.normalize()

                # Compute dot
                leftdot = abs(norm_edge_left_v.dot(norm_edge_v))
                rightdot = abs(norm_edge_right_v.dot(norm_edge_v))
                selfedge = self.edge_left if leftdot < rightdot else self.edge_right

                # Make copies of edges and compute intersection
                self_copy = LineSegment2D(selfedge.p, selfedge.v)
                edge_copy = LineSegment2D(edge.edge.p, edge.edge.v)
                i = intersect_line_segment2d(edge_copy, self_copy)

                if i is not None and not _approximately_equals(i, self.point):
                    # locate candidate b
                    linvec = (self.point - i).normalize()
                    edvec = edge.edge.v.normalize()
                    if linvec.dot(edvec) < 0:
                        edvec = -edvec

                    bisecvec = edvec + linvec
                    if abs(bisecvec) == 0:
                        continue
                    bisector = LineSegment2D(i, bisecvec)
                    b = intersect_line2d(self.bisector, bisector)

                    if b is None:
                        continue

                    # Check eligibility of b.
                    # A valid b should lie within the area limited by the edge
                    # and the bisectors of its two vertices.
                    _left_bisector_norm = edge.bisector_left.v.normalize()
                    _left_b_norm = (b - edge.bisector_left.p).normalize()
                    xleft = _left_bisector_norm.determinant(_left_b_norm) > -self.tol

                    _right_bisector_norm = edge.bisector_right.v.normalize()
                    _right_b_norm = (b - edge.bisector_right.p).normalize()
                    xright = _right_bisector_norm.determinant(_right_b_norm) < self.tol

                    _edge_edge_norm = edge.edge.v.normalize()
                    _b_to_edge_norm = (b - edge.edge.p).normalize()
                    xedge = _edge_edge_norm.determinant(_b_to_edge_norm) < self.tol

                    if not (xleft and xright and xedge):
                        continue  # discarded candidate
                    # found valid candidate
                    _d2b = LineSegment2D(edge.edge.p, edge.edge.v).distance_to_point(b)
                    _new_split_event = _SplitEvent(_d2b, b, self, edge.edge)
                    events.append(_new_split_event)

        i_prev = self.bisector.intersect_line_ray(self.prev.bisector)
        i_next = self.bisector.intersect_line_ray(self.next.bisector)

        if i_prev is not None:
            left_seg = LineSegment2D(self.edge_left.p, self.edge_left.v)
            dist_to_i_prev = left_seg.distance_to_point(i_prev)
            events.append(_EdgeEvent(dist_to_i_prev, i_prev, self.prev, self))
        if i_next is not None:
            right_seg = LineSegment2D(self.edge_right.p, self.edge_right.v)
            dist_to_i_next = right_seg.distance_to_point(i_next)
            events.append(_EdgeEvent(dist_to_i_next, i_next, self, self.next))

        if not events:
            return None

        ev = min(events, key=lambda event:
                 self.point.distance_to_point(event.intersection_point))
        return ev

    def invalidate(self):
        if self.lav is not None:
            self.lav.invalidate(self)
        else:
            self._valid = False

    @property
    def is_valid(self):
        return self._valid

    def __str__(self):
        return "Vertex ({:.2f};{:.2f})".format(self.point.x, self.point.y)

    def __repr__(self):
        return 'Vertex ({}) ({:.2f};{:.2f}), bisector {}, edges {} {}'.format(
            'reflex' if self.is_reflex else 'convex',
            self.point.x, self.point.y, self.bisector,
            self.edge_left, self.edge_right)


class _SLAV:

    def __init__(self, polygon, tol=1e-5):
        """ A set of circular lists of active vertices.

        It stores a loop of vertices for the polygon created during the straight
        skeleton computation (Felkel and Obdrzalek 1998).
        """
        self.tol = tol
        contours = [_normalize_contour(polygon, tol)]
        self._lavs = [_LAV.from_polygon(contour, self) for contour in contours]

        # store original polygon edges for calculating split events
        self._original_edges = [
            _OriginalEdge(
                LineSegment2D(vertex.prev.point, vertex.point),
                vertex.prev.bisector, vertex.bisector
            ) for vertex in chain.from_iterable(self._lavs)
        ]

    def __iter__(self):
        for lav in self._lavs:
            yield lav

    def __len__(self):
        return len(self._lavs)

    def empty(self):
        return len(self._lavs) == 0

    def handle_edge_event(self, event):
        """Resolves adjacency of new edge event.

        This function resolves the edge event with previous edges, LAV, and then stores
        edge information in a Subtree.

        Args:
            event: EdgeEvent

        Returns:
            Subtree namedTuple
        """
        sinks = []
        events = []

        lav = event.vertex_a.lav
        # triangle; one sink point
        if event.vertex_a.prev == event.vertex_b.next:
            self._lavs.remove(lav)
            for vertex in list(lav):
                sinks.append(vertex.point)
                vertex.invalidate()
        else:
            new_vertex = lav.unify(event.vertex_a, event.vertex_b,
                                   event.intersection_point)
            if lav.head in (event.vertex_a, event.vertex_b):
                lav.head = new_vertex
            sinks.extend((event.vertex_a.point, event.vertex_b.point))
            next_event = new_vertex.next_event()
            if next_event is not None:
                events.append(next_event)

        return (Subtree(event.intersection_point, event.distance, sinks), events)

    def handle_split_event(self, event):
        """Consumes a split event.

        This function resolves the adjacency of new split event with previous edges, LAV,
        and then stores edge information in a Subtree.

        Args:
            event: EdgeEvent

        Returns:
            Subtree namedTuple
        """
        lav = event.vertex.lav
        sinks = [event.vertex.point]
        vertices = []
        x = None  # right vertex
        y = None  # left vertex
        norm = event.opposite_edge.v.normalize()
        for v in chain.from_iterable(self._lavs):
            if norm == v.edge_left.v.normalize() and \
                    event.opposite_edge.p == v.edge_left.p:
                x = v
                y = x.prev
            elif norm == v.edge_right.v.normalize() and \
                    event.opposite_edge.p == v.edge_right.p:
                y = v
                x = y.next

            if x:
                xleft = y.bisector.v.normalize().determinant(
                    (event.intersection_point - y.point).normalize()) >= -self.tol

                xright = x.bisector.v.normalize().determinant(
                    (event.intersection_point - x.point).normalize()) <= self.tol

                if xleft and xright:
                    break
                else:
                    x = None
                    y = None

        if x is None:
            return (None, [])

        v1 = _LAVertex(event.intersection_point, event.vertex.edge_left,
                       event.opposite_edge, tol=self.tol)
        v2 = _LAVertex(event.intersection_point, event.opposite_edge,
                       event.vertex.edge_right, tol=self.tol)

        v1.prev = event.vertex.prev
        v1.next = x
        event.vertex.prev.next = v1
        x.prev = v1

        v2.prev = y
        v2.next = event.vertex.next
        event.vertex.next.prev = v2
        y.next = v2

        new_lavs = None
        self._lavs.remove(lav)
        if lav != x.lav:
            # the split event actually merges two lavs
            self._lavs.remove(x.lav)
            new_lavs = [_LAV.from_chain(v1, self)]
        else:
            new_lavs = [_LAV.from_chain(v1, self), _LAV.from_chain(v2, self)]

        for lv in new_lavs:
            if len(lv) > 2:
                self._lavs.append(lv)
                vertices.append(lv.head)
            else:
                sinks.append(lv.head.next.point)
                for v in list(lv):
                    v.invalidate()

        events = []
        for vertex in vertices:
            next_event = vertex.next_event()
            if next_event is not None:
                events.append(next_event)

        event.vertex.invalidate()
        return (Subtree(event.intersection_point, event.distance, sinks), events)


class _LAV:

    def __init__(self, slav):
        """A single circular list of active vertices.

        Stored in a SLAV (Felkel and Obdrzalek 1998, 2).
        """
        self.head = None
        self._slav = slav
        self._len = 0
        self.tol = slav.tol

    @classmethod
    def from_polygon(cls, polygon, slav):
        """Construct a LAV from a list of point coordinates representing a polygon.

        Args:
            polygon: list of points (tuple of x,y coordinates).
            slav: SLAV (a set of circular lists of active vertices).

        Returns:
            LAV (single circular list of active vertices).
        """
        lav = cls(slav)
        for prev, point, next in _window(polygon):
            lav._len += 1
            vertex = _LAVertex(
                point,
                LineSegment2D.from_end_points(prev, point),
                LineSegment2D.from_end_points(point, next),
                tol=slav.tol
            )
            vertex.lav = lav
            if lav.head is None:
                lav.head = vertex
                vertex.prev = vertex.next = vertex
            else:
                vertex.next = lav.head
                vertex.prev = lav.head.prev
                vertex.prev.next = vertex
                lav.head.prev = vertex
        return lav

    @classmethod
    def from_chain(cls, head, slav):
        """Creates new LAV from consumed _LAVertex, and reference _SLAV.

        Args:
            head: Head _LAVertex that creates new _LAV loop.
            slav: Reference SLAV (a set of circular lists of active vertices).

        Returns:
            LAV (a single circular list of active vertices).
        """
        lav = cls(slav)
        lav.head = head
        for vertex in lav:
            lav._len += 1
            vertex.lav = lav
        return lav

    def invalidate(self, vertex):
        """Sets vertex validity to False, and handles head relationship.

        Args:
            vertex: _LAVertex to be invalidated.
        """
        assert vertex.lav is self, 'Tried to invalidate a vertex that is not mine'
        vertex._valid = False
        if self.head == vertex:
            self.head = self.head.next
        vertex.lav = None

    def unify(self, vertex_a, vertex_b, point):
        """Generate new _LAVertex from input Point2D and resolve adjacency to old edges.

        Args:
            vertex_a: _LAVertex
            vertex_b: _LAVertex
            point: Intersection point from angle bisectors from vertex_a and
                vertex_b as a Point2D.

        Returns:
            _LAVertex of intersection point.
        """
        replacement = _LAVertex(
            point, vertex_a.edge_left, vertex_b.edge_right,
            (vertex_b.bisector.v.normalize(), vertex_a.bisector.v.normalize()),
            tol=vertex_a.tol
        )
        replacement.lav = self

        if self.head in [vertex_a, vertex_b]:
            self.head = replacement

        vertex_a.prev.next = replacement
        vertex_b.next.prev = replacement
        replacement.prev = vertex_a.prev
        replacement.next = vertex_b.next

        vertex_a.invalidate()
        vertex_b.invalidate()

        self._len -= 1
        return replacement

    def __str__(self):
        return "LAV {}".format(id(self))

    def __repr__(self):
        return "{} = {}".format(str(self), [vertex for vertex in self])

    def __len__(self):
        return self._len

    def __iter__(self):
        cur = self.head
        while True:
            yield cur
            cur = cur.next
            if cur == self.head:
                return

    def _show(self):
        cur = self.head
        while True:
            print(cur.__repr__())
            cur = cur.next
            if cur == self.head:
                break


class _EventQueue:

    def __init__(self):
        """A priority queue that stores vertices of the polygon."""
        self.__data = []

    def put(self, item):
        if item is not None:
            heapq.heappush(self.__data, item)

    def put_all(self, iterable):
        for item in iterable:
            heapq.heappush(self.__data, item)

    def get(self):
        return heapq.heappop(self.__data)

    def empty(self):
        return len(self.__data) == 0

    def peek(self):
        return self.__data[0]

    def show(self):
        for item in self.__data:
            print(item)


def _merge_sources(skeleton):
    """Merge the sources of a straight skeleton.

    In highly symmetrical shapes with reflex vertices, multiple sources may share
    the same  location. This function merges those sources.
    """
    sources = {}
    to_remove = []
    for i, p in enumerate(skeleton):
        source = tuple(i for i in p.source)
        if source in sources:
            source_index = sources[source]
            # source exists, merge sinks
            for sink in p.sinks:
                if sink not in skeleton[source_index].sinks:
                    skeleton[source_index].sinks.append(sink)
            to_remove.append(i)
        else:
            sources[source] = i
    for i in reversed(to_remove):
        skeleton.pop(i)


def _clean_sinks(skeleton, tol=1e-5):
    """Remove cases where the source and sink point are the same.

    This can occur as an unintended effect of merging sources.
    """
    clean_subtrees = []
    for s_tree in skeleton:
        source, clean_sinks, is_invalid = s_tree.source, [], False
        for sink in s_tree.sinks:
            if source.is_equivalent(sink, tol):
                is_invalid = True
            else:
                clean_sinks.append(sink)
        if is_invalid:  # create a new Subtree
            clean_subtrees.append(Subtree(source, s_tree.height, clean_sinks))
        else:
            clean_subtrees.append(s_tree)
    return clean_subtrees


def _skeletonize(polygon, tol=1e-5):
    """Compute the straight skeleton of a Polygon.

    Args:
        polygon: A list of 2D vertices in clockwise order (when viewed from
            above the XY plane). For example, a square can be represented with
            the following: [[0,1], [1,1], [0,0], [1,0]]
        tol: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        The straight skeleton as a list of "subtrees", which are in the form of
        (source, height, sinks). Source is the highest points, height is its
        height, and sinks are the point connected to the source.
    """
    slav = _SLAV(polygon, tol=tol)
    output = []
    prioque = _EventQueue()

    for lav in slav:
        for vertex in lav:
            v = vertex.next_event()
            prioque.put(v)

    # While the priority queue or SLAV is not empty, compute the intersection of
    # bisectors at edge. The 'handle_edge_event' method computes the next event
    # from the new intersection via the next_event method. Thus, this while loop
    # iteratively adds new events without recursion.
    while not (prioque.empty() or slav.empty()):
        i = prioque.get()  # vertex a, b is self or next vertex
        # Handle edge or split events.
        # arc: subtree(event.intersection_point, event.distance, sinks)
        # events: updated events with new vertex
        if isinstance(i, _EdgeEvent):
            if not i.vertex_a.is_valid or not i.vertex_b.is_valid:
                continue
            (arc, events) = slav.handle_edge_event(i)
        elif isinstance(i, _SplitEvent):
            if not i.vertex.is_valid:
                continue
            (arc, events) = slav.handle_split_event(i)
        prioque.put_all(events)

        # As we traverse priorque, output list of "subtrees", which are in the form
        # of (source, height, sinks) where source is the highest points, height is
        # its distance to an edge, and sinks are the point connected to the source.
        if arc is not None:
            output.append(arc)
    # merge the sources and clean the result
    _merge_sources(output)
    output = _clean_sinks(output, tol=1e-5)
    return output


def _intersect_skeleton_segments(skeleton, tolerance=1e-5):
    """Intersect all of the LineSegment2D of the skeleton and split them.

    Args:
        skeleton: A list of ladybug-geometry LineSegment2D for the segments
            of the straight skeleton.
        tolerance: The tolerance at which the intersection will be computed.

    Returns:
        A tuple with two items.

        * split_skeleton -- A list of LineSegment2D objects for the skeleton
            split through self-intersection.

        * is_intersect_topology -- A boolean value for whether the topology of the
            input skeleton is self-intersecting. This value is True whenever
            there are segments that were split as part of this operation.
    """
    # extend skeleton segments a little to ensure intersections happen
    under_tol = tolerance * 0.99
    ext_skeleton = []
    for seg in skeleton:
        m_v = seg.v.normalize() * under_tol
        ext_seg = LineSegment2D.from_end_points(seg.p1.move(-m_v), seg.p2.move(m_v))
        ext_skeleton.append(ext_seg)

    # compute all of the intersection points across the skeleton
    intersect_pts = [[] for _ in skeleton]
    for i, seg in enumerate(skeleton):
        try:
            for other_seg in ext_skeleton[:i] + ext_skeleton[i + 1:]:
                int_pt = intersect_line_segment2d(seg, other_seg)
                if int_pt is None or int_pt.is_equivalent(seg.p1, tolerance) or \
                        int_pt.is_equivalent(seg.p2, tolerance):
                    continue
                # we have found an intersection point where segments should be split
                intersect_pts[i].append(int_pt)
        except IndexError:
            pass  # we have reached the end of the list

    # loop through the segments and split them at the intersection points
    split_skeleton, is_intersect_topology = [], False
    for seg, split_pts in zip(skeleton, intersect_pts):
        if len(split_pts) == 0:
            split_skeleton.append(seg)
        elif len(split_pts) == 1:  # split the segment in two
            is_intersect_topology = True
            int_pt = split_pts[0]
            split_skeleton.append(LineSegment2D.from_end_points(seg.p1, int_pt))
            split_skeleton.append(LineSegment2D.from_end_points(int_pt, seg.p2))
        else:  # sort the points along the segment to split it
            is_intersect_topology = True
            pt_dists = [seg.p1.distance_to_point(ipt) for ipt in split_pts]
            sort_obj = sorted(zip(pt_dists, split_pts), key=lambda pair: pair[0])
            sort_pts = [x for _, x in sort_obj]
            sort_pts.append(seg.p2)
            pr_pt = seg.p1
            for s_pt in sort_pts:
                if not pr_pt.is_equivalent(s_pt, tolerance):
                    split_skeleton.append(LineSegment2D.from_end_points(pr_pt, s_pt))
                pr_pt = s_pt

    return split_skeleton, is_intersect_topology


def _remove_segments_outside_boundary(skeleton, boundary, tolerance=1e-5):
    """Remove LineSegment2D that are outside the boundary of the parent shape.

    This can be used to clean up the result after intersection of segments
    since the skeleton routines can sometimes fail to produce a result that
    is completely within the parent shape.

    Args:
        skeleton: A list of ladybug-geometry LineSegment2D for the segments
            of the straight skeleton.
        boundary: A Polygon2D for the boundary of the shape. Segments that lie
            outside of this boundary beyond the tolerance will be removed from
            the result.
        tolerance: The tolerance for distinguishing whether skeleton points lie
            outside the boundary.

    Returns:
        A tuple with two items.

        * clean_skeleton -- A list of LineSegment2D objects with segments removed
            that outside of the boundary.
    """
    clean_skeleton = []
    for seg in skeleton:
        p1, p2 = seg.p1, seg.p2
        if boundary.point_relationship(p2, tolerance) >= 0 and \
                boundary.point_relationship(p1, tolerance) >= 0:
            clean_skeleton.append(seg)
    return clean_skeleton


def skeleton_as_subtree_list(boundary, holes=None, tolerance=1e-5):
    """Get a straight skeleton as a list of Subtree source and sink points.

    Args:
        boundary: A ladybug-geometry Polygon2D for the boundary around the shape
            for which the straight skeleton will be computed.
        holes: An optional list of ladybug-geometry Polygon2D for the holes within
            the shape for which a straight skeleton will be computed. If None,
            it will be assumed that no holes exist in the shape. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-5).

    Returns:
        A list of Subtree objects that represent the straight skeleton.
    """
    # merge the holes into the boundary if they are specified
    if holes is not None and len(holes) != 0:
        bound_pts = list(boundary.vertices)
        hole_pts = [list(h.vertices) for h in holes]
        boundary = Polygon2D.from_shape_with_holes(bound_pts, hole_pts)

    # pre-process the boundary to be sure it is in the right order
    boundary = boundary.remove_duplicate_vertices(tolerance)
    if not boundary.is_clockwise:
        boundary = boundary.reverse()

    # skeletonize the objects
    skel_tol = tolerance / 100  # use a finer tolerance for actual skeleton
    subtree_list = _skeletonize(boundary.to_array(), skel_tol)

    # if holes were merged into the boundary, try to clean up the merge seams
    if holes is not None and len(holes) != 0:
        # gather the seams along which holes were merged
        all_segs, seam_segs = boundary.segments, []
        for i, seg in enumerate(all_segs):
            try:
                for o_seg in all_segs[i + 1:]:
                    if seg.p1.is_equivalent(o_seg.p2, tolerance) and \
                            seg.p2.is_equivalent(o_seg.p1, tolerance):
                        seam_segs.append(seg)
                        break
            except IndexError:
                pass  # we have reached the end of the list
        # replace the subtrees along the seam with the seam itself
        for seam in seam_segs:
            s_p1, s_p2 = seam.p1, seam.p2
            new_subtrees, connecting_pts, new_height = [], [], 0
            for sub_tree in subtree_list:
                sink_pts = sub_tree.sinks
                replace_sub = False
                if len(sink_pts) == 2:
                    if s_p1.is_equivalent(sink_pts[0], tolerance) or \
                            s_p1.is_equivalent(sink_pts[1], tolerance):
                        if s_p2.is_equivalent(sink_pts[0], tolerance) or \
                                s_p2.is_equivalent(sink_pts[1], tolerance):
                            new_height = sub_tree.height
                            replace_sub = True
                            connecting_pts.append(sub_tree.source)
                if not replace_sub:
                    new_subtrees.append(sub_tree)
            subtree_list = new_subtrees
            seam_sub = Subtree(seam.midpoint, new_height, connecting_pts + [s_p1, s_p2])
            subtree_list.append(seam_sub)

    return subtree_list


def skeleton_as_edge_list(boundary, holes=None, tolerance=1e-5, intersect=False):
    """Get a straight skeleton as list of LineSegment2D.

    Args:
        boundary: A ladybug-geometry Polygon2D for the boundary around the shape
            for which the straight skeleton will be computed.
        holes: An optional list of ladybug-geometry Polygon2D for the holes within
            the shape for which a straight skeleton will be computed. If None,
            it will be assumed that no holes exist in the shape. (Default: None).
        tolerance: Tolerance for point equivalence. (Default: 1e-5).
        intersect: A boolean to note whether the segments of the skeleton should
            be intersected with one another before being returned. This can
            help make the skeleton more usable in the event that its topology
            is not correct. (Default: False).

    Returns:
        A list of LineSegment2D that represent the straight skeleton.
    """
    # get the Subtree representation of the straight skeleton
    subtree_list = skeleton_as_subtree_list(boundary, holes, tolerance)

    # extract the LineSegment2D from the Subtree representation
    skeleton = []
    for subtree in subtree_list:
        source_pt = subtree.source
        for sink_pt in subtree.sinks:
            edge_seg = LineSegment2D.from_end_points(
                Point2D(source_pt.x, source_pt.y), Point2D(sink_pt.x, sink_pt.y))
            skeleton.append(edge_seg)

    if intersect:  # intersect skeleton segments to split them
        skeleton, _ = _intersect_skeleton_segments(skeleton, tolerance)
        skeleton = _remove_segments_outside_boundary(skeleton, boundary, tolerance)
    return skeleton
