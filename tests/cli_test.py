"""Test cli."""
from click.testing import CliRunner
import json

from ladybug_geometry.geometry3d import Face3D
from ladybug_geometry_polyskel.cli import perimeter_core_subfaces_cli


def test_perimeter_core_subfaces():
    input_file = './tests/assets/sample_geometries.json'
    runner = CliRunner()
    result = runner.invoke(
        perimeter_core_subfaces_cli, [input_file, '1.0', '--tolerance', '0.01'])
    assert result.exit_code == 0

    result_dict = json.loads(result.output)
    new_faces = [Face3D.from_dict(f) for f in result_dict]
    assert len(new_faces) == 33
