#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
import sys
import os

#
if len(sys.argv) != 2:
    print("Usage: {0} <script_name.py>".format(sys.argv[0]))
    sys.exit(-1)
script_filename = sys.argv[1]
script_name, _ = os.path.splitext(script_filename)
script_dir = os.path.dirname(os.path.abspath(__file__))
pytoml_path = os.path.join(script_dir, "pyproject.toml.in")

#
try:
    import tomllib
    # Open and parse the TOML file
    with open(pytoml_path, "rb") as fid:
        data = tomllib.load(fid)
        scripts = data['project']['scripts']
except ImportError:
    with open(pytoml_path, "r") as fid:
        found_scripts = False
        scripts = {}
        for line in fid:
            if found_scripts:
                if line.strip() == "":
                    break
                scripts[line.split("=")[0].strip()] = line.split("=")[1].strip()
            else:
                if line.startswith("[project.scripts]"):
                    found_scripts = True

#
script_template = """#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
''' Stub wrapper for script {0}.{1} '''
from {0} import {1}


if __name__ == "__main__":
    {1}()
"""
try:
    import_path, function_name = scripts[script_name].split(":")
except KeyError:
    print("Error: script name '{0}' not found in pyproject.toml.in".format(script_name))
    sys.exit(-1)
with open(script_filename, "w") as fid:
    fid.write(script_template.format(import_path, function_name))
os.chmod(script_filename, 0o744)