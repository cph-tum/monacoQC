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

# directory of setup files
setup_dir = (script_dir / "../devices").resolve()

# load schema file
with open(script_dir / "device_schema.json", "r") as schemafile:
    schema = json.load(schemafile)

# initialize schema validator
validator = Draft202012Validator(schema)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # loop only over setups given on command line
        setup_files = [setup_dir / f"{n}.yaml" for n in sys.argv[1:]]
    else:
        # loop over all setups
        setup_files = sorted(setup_dir.glob("*.yaml"))

    for setup_file in setup_files:
        print(setup_file)
        with open(setup_file, "r") as setupfile:
            setup = yaml.safe_load(setupfile)
        validator.validate(setup)
