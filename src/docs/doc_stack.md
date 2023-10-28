Internal Stack and Profiling     {#doc_stack}
============================

[TOC]

The Open FUSION Toolkit (OFT) has the ability to maintain its own call stack during execution. This feature can be used to
aid in debugging by narrowing down the location of crashes as well as provide some basic profiling
information. For more information on using this functionality and interpreting output see the \ref
doc_stack_usage section below.

In order to enable stack tracing OFT must be built with supplemental preprocessor
definitions, these can be added and removed using an external Python script `generate_stack.py`. For
more information on how to enable tracing and profiling see \ref doc_stack_build below.

\section doc_stack_usage Using Stack Tracing

When stack tracing is active OFT sets up custom error handlers which will dump the current call stack
on an error, internal abort, and completetion of execution. The stack indicates, from left to right, the
MPI task id, call depth, and subroutine name for each call. An example stack dump is included for reference
below.

\verbatim
Stacktrace
----------------------------------------------
[    0]   4                                         sortarray::oft_sorta4
[    0]   3                                    fem_base::fem_self_linkage
[    0]   2                                           fem_base::fem_setup
[    0]   1                                  oft_lag_basis::oft_lag_setup
\endverbatim

\note A stack trace is dumped for every MPI task which catches and abort signal, this may produce a large
amount of output for parallel runs.

\subsection doc_stack_usage_prof Interpretting Profiling Output

Docs needed

\section doc_stack_build Enabling Tracing

In order to use the stack and profiling functionality OFT must be built/rebuilt with the necessary
definitions. Module and subroutine context information is provided by preprocessor definitions which
identify regions of code. These definitions are created automatically by the Python script `generate_stack.py`.
The script must be run from the base `src` directory as shown in the example below.

\verbatim
~$ python utilities/generate_stack.py
\endverbatim

\note This script must be run whenever subroutines are added or removed, in order to maintain consistent
stack definitions.

Once the definitions have been created/updated OFT must be rebuilt to include the updated definitions.
The `OFT_STACK` directive in the `make_inc.mk` file must also be modified to ensure that the stack mechanics
are compiled into the source. The two example below illustrate the required compiler flags to enable stack
and stack+profiling mechanics.

\verbatim
#---------------------------------------
# Use stack tracing only
#---------------------------------------
STACK_TRACE = -DOFT_STACK

#---------------------------------------
# Use stack tracing and create profiling information
#---------------------------------------
STACK_TRACE = -DOFT_STACK -DOFT_PROFILE
\endverbatim

\subsection doc_stack_build_clean Removing Stack Definitions

Before changes are commited to a local or remote repository stack definitions should be removed from the
source files. This can be done by running the `generate_stack.py` script with the `-c`flag to indicating
definition cleaning should be performed. The script should be run from the main `src` directory as above.

\verbatim
~$ python utilities/generate_stack.py -c
\endverbatim
