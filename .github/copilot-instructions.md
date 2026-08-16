# Copilot Instructions for OpenFUSIONToolkit

This file contains Copilot-tailored instructions for the OpenFUSIONToolkit project. General information about the OpenFUSIONToolkit project can be found in the `AGENTS.md` file in the root of the repository.

## Common Pitfalls

1. **Source setup_env.sh**: Always source this before any build/test commands. It activates the Python venv. In CI, this file is generated during prerequisites setup.

2. **Build from `builds/` directory**: The `build_libs.py` script must run from the `builds/` directory. It creates `config_cmake.sh` there. CMake then creates `build_release/` and `install_release/` inside `builds/`.