# Copilot Instructions for OpenFUSIONToolkit

This file contains Copilot-tailored instructions for the OpenFUSIONToolkit project. General information about the OpenFUSIONToolkit project can be found in the `AGENTS.md` file in the root of the repository.

## Code review

When reviewing code, please consider and explicitly highlight any issues in the following:
 - **Documentation**: All new or modified functions, classes, and modules should be well-documented using Doxygen-compatible comments. Check for missing or outdated documentation in the codebase and associated documentation files within `src/docs`. Ensure that the code is well-documented and easy to understand. Significant new functionality should be accompanied by an example.

 - **Code Style:** Ensure that the code adheres to the project's coding standards and style guidelines. This includes consistent naming conventions, indentation, and formatting. Check for any deviations from the established code style and provide suggestions for improvement.

 - **Tests:** Significant new functionality or changes should be accompanied by appropriate tests. Check for the presence of relevant test cases.

## Common Pitfalls

1. **Source setup_env.sh**: Always source this before any build/test commands. It activates the Python venv. In CI, this file is generated during prerequisites setup.

2. **Build from `builds/` directory**: The `build_libs.py` script must run from the `builds/` directory. It creates `config_cmake.sh` there. CMake then creates `build_release/` and `install_release/` inside `builds/`.