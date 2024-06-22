# coding=utf-8
"""Functions for offsetting a polygon using the straight skeleton."""
from __future__ import division

from .polygraph import skeleton_as_directed_graph
from .polysplit import perimeter_core_subpolygons, _split_perimeter_subpolygons


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
    # if distance is negative, use the normal offset method and check for intersection
    if distance <= 0:
        offset_poly = polygon.offset(distance, check_intersection=False)
        if offset_poly.is_self_intersecting:
            off_polys = offset_poly.split_through_self_intersection(tolerance)
            off_polys = sorted(off_polys, key=lambda x: x.area, reverse=True)
            offset_poly = off_polys[0]
        return [offset_poly]
    # if the distance is positive, use the straight skeleton 
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
    # if distance is negative, use the normal offset method and check for intersection
    if distance == 0:
        return []
    elif distance < 0:  # we are offsetting outwards
        polygon = polygon.offset(distance, check_intersection=False)
        if polygon.is_self_intersecting:
            off_polys = polygon.split_through_self_intersection(tolerance)
            off_polys = sorted(off_polys, key=lambda x: x.area, reverse=True)
            polygon = off_polys[0]
        distance = abs(distance)
    # initialize sub polygon lists
    perimeter_polygons = []
    # compute the straight skeleton of the polygon as a directed graph
    dg = skeleton_as_directed_graph(polygon, None, tolerance)
    # traverse the directed graph to get all of the perimeter polygons
    _perimeter_sub_dg = None
    root_keys = [dg.outer_root_key] + dg.hole_root_keys
    for root_key in root_keys:
        # compute the polygons on the perimeter of the polygon
        _perimeter_sub_polys, _perimeter_sub_dg = _split_perimeter_subpolygons(
            dg, distance, root_key, tolerance, _perimeter_sub_dg)
        perimeter_polygons.extend(_perimeter_sub_polys)  # collect perimeter sub-polys
    return perimeter_polygons
