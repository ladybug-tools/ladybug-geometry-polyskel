# coding=utf-8
"""Functions for offsetting a polygon using the straight skeleton."""
from __future__ import division
import math

from ladybug_geometry.geometry2d import Point2D, Ray2D, Polygon2D

from .polygraph import skeleton_as_directed_graph
from .polysplit import perimeter_core_subpolygons, _split_perimeter_subpolygons, \
    _exterior_cycles_as_polygons


def offset_polygon(polygon, distance, tolerance=1e-5):
    """Offset a Polygon2D by a given distance inwards or outwards.

    Assuming that the input is not self-intersecting, the resulting shapes
    should have all self-intersections resolved through the straight skeleton.
    However, if the distance is large enough, the result may be an empty list.

    Args:
        polygon: A Polygon2D to be offset inwards or outwards by a given distance.
        distance: The distance inwards that the polygon will be offset.
            Positive values will always be offset inwards while negative ones
            will be offset outwards.
        tolerance: Tolerance for point equivalence. (Default: 1e-5).
    """
    # determine the methods to use based on the offset direction
    if -tolerance <= distance <= tolerance:
        return [polygon]
    elif distance < 0:  # we are offsetting outwards
        offset_poly = polygon.offset(distance, check_intersection=False)
        if not offset_poly.is_self_intersecting:
            return [offset_poly]
        bound_rect = _bound_rect_for_ext_skeleton(polygon, distance)
        dg = skeleton_as_directed_graph(bound_rect, [polygon], tolerance)
        root_key = dg.hole_root_keys[0]
        _, _perimeter_sub_dg = _split_perimeter_subpolygons(
            dg, abs(distance), root_key, tolerance, None)
        offset_polygons = _exterior_cycles_as_polygons(_perimeter_sub_dg, tolerance)
    else:  # we are offsetting inwards
        _, offset_polygons = perimeter_core_subpolygons(
            polygon, distance, tolerance=tolerance, flat_core=True)
    return offset_polygons


def offset_perimeter_polygons(polygon, distance, tolerance=1e-5):
    """Get a list of Polygon2Ds along the perimeter of an input polygon.

    These can be used to split the input polygons into core/perimeter shapes.
    Assuming that the input is not self-intersecting, the resulting shapes
    should have all self-intersections resolved through the straight skeleton.

    Args:
        polygon: A Polygon2D to be offset inwards or outwards by a given distance.
        distance: The distance inwards that the polygon will be offset.
            Positive values will always be offset inwards while negative ones
            will be offset outwards.
        tolerance: Tolerance for point equivalence. (Default: 1e-5).
    """
    # determine the methods to use based on the offset direction
    if -tolerance <= distance <= tolerance:
        return []
    elif distance < 0:  # we are offsetting outwards
        offset_poly = polygon.offset(distance, check_intersection=False)
        if not offset_poly.is_self_intersecting:
            perimeter_polygons = []
            for i, (in_pt, out_pt) in enumerate(zip(polygon, offset_poly)):
                pts = (polygon[i - 1], offset_poly[i - 1], out_pt, in_pt)
                perimeter_polygons.append(Polygon2D(pts))
            return perimeter_polygons
        bound_rect = _bound_rect_for_ext_skeleton(polygon, distance)
        dg = skeleton_as_directed_graph(bound_rect, [polygon], tolerance)
        root_key = dg.hole_root_keys[0]
        perimeter_polygons, _ = _split_perimeter_subpolygons(
            dg, abs(distance), root_key, tolerance, None)
    else:  # we are offsetting inwards
        dg = skeleton_as_directed_graph(polygon, None, tolerance)
        root_key = dg.outer_root_key
        perimeter_polygons, _ = _split_perimeter_subpolygons(
            dg, distance, root_key, tolerance, None)
    return perimeter_polygons


def offset_skeleton(polygon, distance, tolerance=1e-5):
    """Get a list of LineSegment2D for the straight used to offset a polygon.

    Args:
        polygon: A Polygon2D to be offset inwards or outwards by a given distance.
        distance: The distance inwards that the polygon will be offset.
            Positive values will always be offset inwards while negative ones
            will be offset outwards.
        tolerance: Tolerance for point equivalence. (Default: 1e-5).
    """
    if -tolerance <= distance <= tolerance:
        return []
    elif distance < 0:  # we are offsetting outwards
        bound_rect = _bound_rect_for_ext_skeleton(polygon, distance)
        dg = skeleton_as_directed_graph(bound_rect, [polygon], tolerance)
    else:  # we are offsetting inwards
        dg = skeleton_as_directed_graph(polygon, None, tolerance)
    return dg.connection_segments


def _bound_rect_for_ext_skeleton(polygon, distance):
    """Get a bounding rectangle to be used for offsetting a polygon outwards.

    Args:
        polygon: A Polygon2D that is being offset outwards for which a bounding
            rectangle is needed.
        distance: The distance outwards that the polygon will be offset. This
            value should always be negative.

    Returns:
        A Polygon2D for the bounding rectangle around the shape.
    """
    # determine how far from the polygon the bounding rectangle will be drawn
    min_pt, max_pt = polygon.min, polygon.max
    poly_dim = min((max_pt.x - min_pt.x, max_pt.y - min_pt.y))
    offset_factor = 50 if abs(distance) < poly_dim / 10 else 10
    offset = abs(distance) * offset_factor

    # create a starting bounding rectangle from the min, max, and distance
    min_pt = Point2D(min_pt.x - offset, min_pt.y - offset)
    max_pt = Point2D(max_pt.x + offset, max_pt.y + offset)
    min_x_max_y, max_x_min_y = Point2D(min_pt.x, max_pt.y), Point2D(max_pt.x, min_pt.y)
    rect_pts = (min_pt, max_x_min_y, max_pt, min_x_max_y)
    rect_poly = Polygon2D(rect_pts)

    # find a corner point to be used as the base for the rectangle
    pts_dists, pts_info = [], []
    for corner_i, conrner_pt in enumerate(rect_poly):
        for pt_i, pt in enumerate(polygon.vertices):
            pts_dists.append(conrner_pt.distance_to_point(pt))
            pts_info.append((pt_i, corner_i))
    sort_info = [inf for _, inf in sorted(zip(pts_dists, pts_info),
                                          key=lambda pair: pair[0])]
    poly_pt_i, corner_pt_i = sort_info[0]

    # shoot a ray from the polygon point to the rectangle
    poly_pt = polygon.vertices[poly_pt_i]
    v1 = polygon.vertices[poly_pt_i - 1] - poly_pt
    end_i = poly_pt_i + 1 if poly_pt_i != len(polygon) - 1 else 0
    v2 = polygon.vertices[end_i] - poly_pt
    if not polygon.is_clockwise:
        ang = v1.angle_clockwise(v2) / 2
        if ang == 0:
            ang = math.pi / 2
        m_vec = v1.rotate(-ang).normalize()
        m_dist = -distance / math.sin(ang)
    else:
        ang = v1.angle_counterclockwise(v2) / 2
        if ang == 0:
            ang = math.pi / 2
        m_vec = v1.rotate(ang).normalize()
        m_dist = distance / math.sin(ang)
    m_vec = m_vec * m_dist
    corner_ray = Ray2D(poly_pt, m_vec)
    ray_int = rect_poly.intersect_line_ray(corner_ray)

    # if the ray intersection is too far from the corner, just use the start rectangle
    corner_pt = rect_pts[corner_pt_i]
    safe_distance = abs(distance) * (offset_factor - 3)
    if len(ray_int) == 0 or ray_int[0].distance_to_point(corner_pt) > safe_distance:
        return rect_poly

    # otherwise, adjust the bounding rectangle to give us a nicer skeleton
    bpt = ray_int[0]
    base, height = max_pt.x - min_pt.x, max_pt.y - min_pt.y
    if corner_pt_i == 0:  # lower left corner
        rect_pts = (
            bpt, Point2D(bpt.x + base, bpt.y), Point2D(bpt.x + base, bpt.y + height),
            Point2D(bpt.x, bpt.y + height))
    elif corner_pt_i == 1:  # lower right corner
        rect_pts = (
            Point2D(bpt.x - base, bpt.y), bpt, Point2D(bpt.x, bpt.y + height),
            Point2D(bpt.x - base, bpt.y + height))
    elif corner_pt_i == 2:  # upper right corner
        rect_pts = (
            Point2D(bpt.x - base, bpt.y - height), Point2D(bpt.x, bpt.y - height),
            bpt, Point2D(bpt.x - base, bpt.y))
    else:  # upper left corner
        rect_pts = (
            Point2D(bpt.x, bpt.y - height), Point2D(bpt.x + base, bpt.y - height),
            Point2D(bpt.x + base, bpt.y), bpt)
    return Polygon2D(rect_pts)
