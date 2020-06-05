# coding=utf-8
"""Classes for computing straight skeleton for 2D concave polygons."""
from __future__ import division

from ladybug_geometry_polyskel import polyskel
from ladybug_geometry.geometry2d.polygon import Polygon2D
from ladybug_geometry.geometry2d.pointvector import Point2D, Vector2D
from ladybug_geometry.geometry3d.pointvector import Vector3D, Point3D
from ladybug_geometry.geometry2d.line import LineSegment2D
from ladybug_geometry_polyskel.polygon_directed_graph import \
    PolygonDirectedGraph, _vector2hash

# from pprint import pprint as pp
# import math

# from astrobot.utils import *
# TODO: Move this to Polygon2D?
def interior_angles(polygon, radian=True):
    """Compute the interior angles of a polygon, accounting for reflex angles.

    Args:
        polygon: A Polygon2D.
        radian: Returns angles in radians if true else degree (default: True).

    Returns:
        generator of angles.
    """
    segments = list(polygon.segments)
    segments = segments + [segments[0]]

    for i in range(len(segments) - 1):
        seg1 = segments[i]
        seg2 = segments[i+1]

        # Get 3d direction vectors
        vec1 = Vector3D(*(seg1.p2 - seg1.p1).to_array(), 0)
        vec2 = Vector3D(*(seg2.p2 - seg2.p1).to_array(), 0)

        theta = math.pi if vec1.cross(vec2).z < 0 else 0.0

        # Flip vec2 to correct angle for calc
        theta += vec1.angle(vec2 * -1)

        if not radian:
            theta = theta * 180 / math.pi

        yield theta


def infinite_line2d_intersection(seg1, seg2):
    """Intersects two line2D assuming they are infinite.

    This method computes the linear equation from the segments.
    Then solves the linear system to find the intersection for both:
    A * x = C
    x = A^-1 * C
    """

    # Compute the A,B coefficients
    seg1_3d = Point3D(*(seg1.p2 - seg1.p1).to_array(), 0)
    seg2_3d = Point3D(*(seg2.p2 - seg2.p1).to_array(), 0)

    z = Point3D(0, 0, 1)
    norm1 = seg1_3d.cross(z)
    norm2 = seg2_3d.cross(z)
    norm1 = norm1 / norm1.magnitude
    norm2 = norm2 / norm2.magnitude

    # Simplify naming for the A matrix (linear equation coefficeints)
    _a, _b, _c, _d = norm1.x, norm1.y, norm2.x, norm2.y

    # Substitute point in each line to solve for C array
    c1 = (_a * seg1.p2.x) + (_b * seg1.p2.y)
    c2 = (_c * seg2.p2.x) + (_d * seg2.p2.y)

    # Check the determinate
    det = (_a * _d) - (_b * _c)
    if abs(det) < 1e-10:
        return

    # Solve for x my multiplying Ainv by C
    x = [( _d/det * c1) + (-_b/det * c2),
         (-_c/det * c1) + ( _a/det * c2)]
    intersection_pt = Point2D.from_array(x)

    return intersection_pt


def _offset_bisector(seg1, seg2, distance):
    """ Calculates the magnitude of the offset bisector."""

    p1 = infinite_line2d_intersection(seg1, seg2)

    # Normal offset outward
    seg1_3d = Point3D(*(seg1.p2 - seg1.p1).to_array(), 0)
    seg2_3d = Point3D(*(seg2.p2 - seg2.p1).to_array(), 0)
    z = Point3D(0, 0, 1)
    norm1 = seg1_3d.cross(z)
    norm2 = seg2_3d.cross(z)

    # Get offset
    off1 = seg1.move((norm1 / norm1.magnitude) * distance)
    off2 = seg2.move((norm2 / norm2.magnitude) * distance)
    p2 = infinite_line2d_intersection(off1, off2)
    if p1 is None or p2 is None:
        return

    return LineSegment2D.from_end_points(p1, p2)


def offset(polygon, distance, tol=1e-10):
    """Offset the polygon boundary by defined distance.

    Args:
        polygon: A Polygon2D object to offset.
        distance: Distance to offset. Only positive distance works in
            current implementation.
        tol: Tolerance for point equivalence.

    Returns:
        A list of offset contours as Polygon2D objects.
    """
    _buffer_angle_tol = math.pi/4  # 45 degrees

    # If positive do inward offset
    if distance > 0:
        _, offsets = polyskel.sub_polygons(polygon, distance, [], tol)
        return offsets

    # Init distance
    distance = abs(distance)
    sqrdist = distance * distance

    # Calculate buffer offset
    segments = list(polygon.segments)
    segments = segments + [segments[0]]

    buffer = 0.0
    for i in range(len(segments)-1):
        seg1, seg2 = segments[i], segments[i+1]
        bisector = _offset_bisector(seg1, seg2, distance)
        if bisector.length > buffer:
            buffer = bisector.length
    buffer += (distance/2.)

    # Make buffer
    buffer_frame = Polygon2D([
        Point2D(polygon.max.x + buffer, polygon.max.y + buffer),
        Point2D(polygon.min.x - buffer, polygon.max.y + buffer),
        Point2D(polygon.min.x - buffer, polygon.min.y - buffer),
        Point2D(polygon.max.x + buffer, polygon.min.y - buffer)])

    holes = [list(reversed(polygon.vertices))]
    perimeter_sub_polygon = polyskel.sub_polygons(buffer_frame, distance, holes)
    mispqll = ""

    # Make a graph from the perimeter (offset) polygons so we infer the core polygons
    g = PolygonDirectedGraph(tol)
    for poly in perimeter_sub_polygon:
        pts = poly.vertices
        for i in range(len(pts)-1):
            g.add_node(pts[i], [pts[i + 1]])
        g.add_node(pts[-1], [pts[0]])

    offset = Polygon2D([node.pt for node in g.exterior_cycles[-1]])

    return offset, perimeter_sub_polygon