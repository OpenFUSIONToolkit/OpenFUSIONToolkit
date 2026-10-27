# Copilot Instructions for OpenFUSIONToolkit

This file contains Copilot-tailored instructions for the OpenFUSIONToolkit project. General information about the OpenFUSIONToolkit project can be found in the `AGENTS.md` file in the root of the repository.

## Code review

When reviewing code, please consider and explicitly highlight any issues in the following:
 - **Documentation**: New or modified functions, classes, and modules in code that Doxygen processes should be well-documented using Doxygen-compatible comments. Check for missing or outdated documentation in the codebase and associated documentation files within `src/docs`. Ensure that the code is well-documented and easy to understand. Significant new functionality should be accompanied by an example.

   **Exception — test and example code:** follow "Test and example code" in `AGENTS.md`, which is the authoritative wording. In short: example sources under `src/examples/` carry documentation markup that is parsed into the published pages and must be left alone, while functions elsewhere under `src/tests/` and `src/examples/` are not parsed into the API documentation and do not need Doxygen-style documentation blocks — do not flag their absence.

 - **Code Style:** Ensure that the code adheres to the project's coding standards and style guidelines. This includes consistent naming conventions, indentation, and formatting. Check for any deviations from the established code style and provide suggestions for improvement.

 - **Tests:** Significant new functionality or changes should be accompanied by appropriate tests. Check for the presence of relevant test cases.

## Common Pitfalls

1. **Source setup_env.sh (ephemeral environments only)**: In the Copilot agent environment, which is provisioned by `copilot-setup-steps.yml`, `setup_env.sh` is generated during prerequisites setup and activates the Python venv at `oft_venv/` — always source it before any build/test command. Neither file is created by `build_libs.py`, and neither exists in a user's local checkout; see "Python environment" in `AGENTS.md` for that case.

2. **Build from `builds/` directory (ephemeral environments only)**: In the Copilot agent environment, `build_libs.py` runs from `builds/`, creating `config_cmake.sh` there, and CMake then creates `build_release/` and `install_release/` inside `builds/`. User checkouts commonly nest one build tree per configuration (eg. `builds/GCC/MPI/`); see "Locating the build tree" in `AGENTS.md` before assuming these paths.