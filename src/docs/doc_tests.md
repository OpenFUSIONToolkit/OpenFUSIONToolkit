Regression Tests     {#doc_tests}
================

[TOC]

The Open FUSION Toolkit includes a set of regression tests which are used to monitor the main API automatically. These
tests are located in sub-folders `src/*/tests` and documented in \ref testing "regression tests". Tests
are run automatically each time the master branch is updated and ensure that certain basic cases produce
the expected results.

\section doc_tests_run Running Tests

The tests are setup to be run as a batch using the Python testing framework
[pytest](https://docs.pytest.org/en/latest/). With pytest installed the tests may be built and run from the
`src` folder with the following commands. The flag `"-m 'not slow'"` filters the tests so that time consuming
tests, such as those for time-dependent MHD modules, are not run.

\verbatim
~$ make test
~$ pytest -m 'not slow' grid lin_alg fem physics
\endverbatim

\section doc_tests_write Writing Tests

Documentation needed
