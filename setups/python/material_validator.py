# monacoQC: An object-oriented Matlab-based device engineering tool for
# quantum cascade structures.
#
# Copyright (C) 2025, Computational Photonics Group, Technical University of
# Munich
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import json
import yaml
import sys

from jsonschema import Draft202012Validator
from pathlib import Path


# directory of this script file
script_dir = Path(__file__).parent

# directory of material files
material_dir = (script_dir / "../materials").resolve()

# load schema file
with open(script_dir / "material_schema.json", "r") as schemafile:
    schema = json.load(schemafile)

# initialize schema validator
validator = Draft202012Validator(schema)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # loop only over materials given on command line
        material_files = [material_dir / f"{n}.yaml" for n in sys.argv[1:]]
    else:
        # loop over all materials
        material_files = sorted(material_dir.glob("*.yaml"))

    for material_file in material_files:
        print(material_file)
        with open(material_file, "r") as matfile:
            material = yaml.safe_load(matfile)
        validator.validate(material)
