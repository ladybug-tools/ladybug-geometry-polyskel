"""ladybug-geometry-polyskel commands."""
import click
import sys
import logging
import json

from ladybug_geometry.geometry3d import Face3D
from ladybug_geometry_polyskel.polysplit import perimeter_core_subfaces as \
    lgp_perimeter_core_subfaces

_logger = logging.getLogger(__name__)


# command group for all geometry extension commands.
@click.group(help='ladybug geometry polyskel commands.')
@click.version_option()
def main():
    pass


@main.command('perimeter-core-subfaces')
@click.argument('face3d-file', type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.argument('perimeter-offset', type=float)
@click.option(
    '--tolerance', '-t', help='The maximum difference between x, y, and z '
    'values at which vertices are considered equivalent.',
    type=float, default=1e-5)
@click.option(
    '--output-file', help='Optional file to output the he string of the visualization '
    'file contents. By default, it will be printed out to stdout',
    type=click.File('w'), default='-', show_default=True)
def perimeter_core_subfaces_cli(face3d_file, perimeter_offset, tolerance, output_file):
    """Compute perimeter and core Face3Ds using the straight skeleton of input Face3Ds.

    \b
    Args:
        face3d_file: Full path to a JSON file containing an array of Face3D
            objects for which core/perimeter Face3Ds will be computed.
        perimeter_offset: Distance to offset perimeter sub-faces.
    """
    try:
        perimeter_core_subfaces(face3d_file, perimeter_offset, tolerance, output_file)
    except Exception as e:
        _logger.exception('Failed to compute core/perimeter Face3Ds.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)


def perimeter_core_subfaces(
        face3d_file, perimeter_offset, tolerance=1e-5, output_file=None
):
    """Compute perimeter and core Face3Ds using the straight skeleton of input Face3Ds.

    Args:
        face3d_file: Full path to a JSON file containing an array of Face3D
            objects for which core/perimeter Face3Ds will be computed.
        perimeter_offset: Distance to offset perimeter sub-faces.
        tolerance: The maximum difference between x, y, and z values at which
            vertices are considered equivalent.
        output_file: Optional file to output the string of the visualization
            file contents. If None, the string will simply be returned from
            this method.
    """
    # load the Face3D objects
    with open(face3d_file) as inf:
        data = json.load(inf)
    face_3ds = [Face3D.from_dict(geo) for geo in data]

    # compute core/perimeter Face3Ds
    core_per_geos = []
    for geo_obj in face_3ds:
        perimeter_sub_faces, core_sub_faces = \
            lgp_perimeter_core_subfaces(geo_obj, perimeter_offset, tolerance)
        core_per_geos.extend([f.to_dict() for f in perimeter_sub_faces])
        core_per_geos.extend([f.to_dict() for f in core_sub_faces])

    # write the resulting geometries into the output
    if output_file is None:
        return json.dumps(core_per_geos)
    else:
        output_file.write(json.dumps(core_per_geos))
