# How to contribute to the Open FUSION Toolkit (OFT)

## Reporting issues and suggesting features

### Do you have questions about using OFT?
TODO: Create a suitable place for questions/discussion (eg. Discord).

### Did you find a bug?
Before opening an issue:
* **Ensure this is a bug!** While we make an effort to provide a robust and clear tool there are many nuances that still exist. Please view the [FAQs on the wiki](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/wiki) and other documentation to ensure the behavior is truly unexpected.
* **Ensure the bug was not already reported** by searching on GitHub under [Issues](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/issues).

If you're unable to find an open issue addressing the problem, [open a new one](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/issues/new). Be sure to include a **title and clear description** and as much relevant information as possible, but *at least* the following:
 * OS type and version where the bug occurs (eg. Ubuntu 20.04)
 * CPU for machine where bug occurs (eg. Apple M2)
 * Whether you are using a binary or built the code from source

Additionally, if at all possible please include **sample code** or input files that reproduce the bug and/or unexpected behavior. If you do not want to provide example code on GitHub but would be willing to share privately with a developer, please note this in the issue so a developer may contact you.

### Are you suggesting a feature?
Please note that in general feature suggestions/requests will not be able to be supported unless a developer is willing to volunteer or the feature is already within the scope of a developer's planned work. Therefore **features should only be suggested if the requesting party is planning to primarily develop the feature themselves**.

Before opening an issue:
* **Ensure the feature is not already planned** by searching on GitHub under [Issues](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/issues).

If you unable to find an open issue related to the suggested feature **and you plan to primarily develop the feature yourself**, [open a new issue](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/issues/new) with the `enhancement` tag.

## Contributing new code, examples, and documentation

All contributions to OFT should be made via pull requests from a fork of the main repository. Specific guidance on types of pull requests is given below, but in general all pull request should include the following:
 * A detailed description of the changes being made and why they are being made
 * Any ancillary changes to the core code that were required for these changes
 * State explicitly whether or not the changes modify any APIs or input files

### Did you write a patch that fixes a bug?

Before requesting a merge of your fix, please ensure the following:
 * You have run the regression tests on your local platform and the fix did not cause any new failures
 * If the new code adds new module variables, adds new subroutines/functions, or changes a subroutine/function signature, ensure that all changes are documented for inclusion in the Doxygen-generated documentation.
 * If the bug was a result of an untested edge (but reasonable) case and/or your code adds or changes functionality, ensure that it has a regression test that covers its functionality in at least some way

Once the above is satisfied, 
* Open a new GitHub pull request with the patch.
* Ensure the PR description clearly describes the problem and solution. Include the relevant issue number if applicable.


### Do you intend to add a new feature or change an existing one?

Suggest your change as a [new issue](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/issues) with the `enhancement` tag. Although you are free to start writing code and testing ahead of time, if you'd like feedback first please provide at least one week for response from one of the main developers.

Before requesting a merge of your fix, please ensure the following:
 * Ensure that only changes **necessary for the feature or change to function properly** are included in the patch
 * You have run the regression tests on your local platform and the fix did not cause any new failures
 * Please ensure that all changes are documented for inclusion in the Doxygen-generated documentation.
 * Please ensure that it has a regression test that covers its functionality in at least some way

When ready to merge, open a new GitHub pull request with the patch. Be sure to monitor the result of the CI builds for your new pull request and, if possible, fix any issues that arise in the initial builds.

### Do you want to contribute to the documentation and examples?

* Please read [writing documentation for OFT](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/wiki/Writing-documentation-for-OFT) or [adding an example to OFT](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/wiki/Adding-an-example-to-OFT).

### Did you fix whitespace, format code, or make a purely cosmetic patch?

Changes that are cosmetic in nature and do not add anything substantial to the stability, functionality, or testability of the Open FUSION Toolkit will generally not be accepted.

## Thank You!

The Open FUSION Toolkit is a community effort. We encourage you to pitch in in any way you can!

Thanks! :sunny: :heart:

The Open FUSION Toolkit Team

## Contributors
Daniel Burgess (@d-burg)\
Sophia Guizzo (@sguizzo)\
Chris Hansen (@hansec)\
Julia Kirby (@juliagoestoikea)\
Francois Logak\
Matthew Pharr (@matt-pharr)

### Contributors to predecessor codes
Thomas Benedett\
Alan Kaptanoglu\
George Marklin\
Derek Sutherland