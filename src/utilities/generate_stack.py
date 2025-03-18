#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Script for handling creation of Open FUSION Toolkit (OFT) internal stack representation for
# usage with debug and profiling.
#
# Note: Script must be run from the root "src" directory.
#
#------------------------------------------------------------------------------
from __future__ import print_function
import sys
import os
import pathlib
import re
# Regular expressions
TEST_LINE_REGEX = re.compile(r'[ ]*(PURE|IMPURE|ELEMENTAL|RECURSIVE|PROGRAM|MODULE|SUBROUTINE|FUNCTION)+[ ]+', re.I)
SUB_MOD_REGEX = re.compile(r'[ ]*(PURE|IMPURE|ELEMENTAL|RECURSIVE)+', re.I)
SUB_REGEX = re.compile(r'[ ]*(SUBROUTINE|FUNCTION)[ ]+([a-z0-9_]+)', re.I)
END_SUB_REGEX = re.compile(r'[ ]*END[ ]*(SUBROUTINE|FUNCTION)\s', re.I)
MOD_REGEX = re.compile(r'[ ]*(PROGRAM|MODULE)[ ]+([a-z0-9_]+)', re.I)
DEF_REGEX = re.compile(r'(#define|#undef) DEBUG_STACK_', re.I)
#----------------------------------------------------------------
# Module class for FORTRAN subroutine containers
#----------------------------------------------------------------
class module:
    def __init__(self,name,ind):
        self.name=name
        self.ind=ind

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        return str(self.name)

    def __eq__(self,other):
        return self.name==other.name
#----------------------------------------------------------------
# Function class for FORTRAN subroutines
#----------------------------------------------------------------
class function:
    def __init__(self,name,module,ind):
        self.name=name
        self.module=module
        self.ind=ind

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        if self.module == "":
            return str(self.name)
        else:
            return str(self.name)

    def __eq__(self,other):
        return (self.name==other.name and self.module==other.module)
#----------------------------------------------------------------
# Parser to analyze FORTRAN source and determine modules and subroutines
#----------------------------------------------------------------
def parse_fortran_file(fid,modules,functions,debug=False):
    def check_for_obj(line):
        test_match = TEST_LINE_REGEX.match(line)
        if test_match is None:
            return None, None
        # Detect Programs/Modules
        mod_match = MOD_REGEX.match(line)
        if mod_match is not None:
            name = mod_match.group(2)
            if name.lower() != 'procedure':
                return name.lower(), 'mod'
        # Remove modifiers
        modifiers_match = SUB_MOD_REGEX.match(line)
        while modifiers_match is not None:
            line = line[modifiers_match.end(0):]
            modifiers_match = SUB_MOD_REGEX.match(line)
        # Detect Subroutines/Functions
        sub_match = SUB_REGEX.match(line)
        if sub_match is not None:
            name = sub_match.group(2)
            return name.lower(), 'sub'
        return None, None
    # Main body
    file_buffer = ""
    current_module = ""
    current_function = "default"
    fun_depth = 0
    line_number = 0
    entry_count = 0
    mod_id = 1
    for line in fid:
        line_number += 1
        if (line.count("#define DEBUG_STACK_") > 0) or (line.count("#undef DEBUG_STACK_") > 0):
            continue
        if line.startswith("DEBUG_STACK_PUSH"):
            entry_count += 1
            if entry_count > 1:
                print('  WARNING: Multiple entry points to subroutine "{0}" at line {1}'.format(current_function, line_number))
        file_buffer += line
        if fun_depth > 0:
            end_match = END_SUB_REGEX.match(line)
            if end_match is not None:
                if fun_depth==1 and entry_count==0:
                    functions.pop()
                    prev_def_start = file_buffer.rindex("#undef DEBUG_STACK_SUB_IND")
                    prev_def_end = prev_def_start + 2 + file_buffer[prev_def_start:].index("\n")
                    prev_def_end += file_buffer[prev_def_end:].index("\n")
                    file_buffer = file_buffer[:prev_def_start] + file_buffer[prev_def_end+1:]
                fun_depth -= 1
                continue
        name, type = check_for_obj(line)
        if name is None:
            continue
        if type == "mod":
            if debug:
                print('  Module "{0}" at line {1}'.format(name,line_number))
            current_module = name
            ind = len(modules) + 1
            new_module = module(current_module,ind)
            if new_module in modules:
                print('  WARNING: Module "{0}" declared more than once at line {1}'.format(current_module, line_number))
            modules.append(new_module)
            mod_id = ind
            file_buffer += "#undef DEBUG_STACK_MOD_IND\n#define DEBUG_STACK_MOD_IND " + str(ind) + "\n"
        elif type=="sub":
            if debug:
                if fun_depth > 0:
                    print('  Skipping contained subroutine "{0}" at line {1}'.format(name,line_number))
                else:
                    print('  Subroutine "{0}" at line {1}'.format(name,line_number))
            if fun_depth > 0:
                fun_depth += 1
                continue
            current_function = name
            ind = len(functions) + 1
            new_function = function(current_function,mod_id,ind)
            if new_function in functions:
                print('  WARNING: Subroutine "{0}" declared more than once at line {1}'.format(current_function,line_number))
            functions.append(new_function)
            file_buffer += "#undef DEBUG_STACK_SUB_IND\n#define DEBUG_STACK_SUB_IND " + str(ind) + "\n"
            entry_count = 0
            fun_depth = 1
    return modules, functions, file_buffer
#----------------------------------------------------------------
# Create stack variables in FORTRAN syntax
#----------------------------------------------------------------
def create_debug_list(modules,functions):
    file_buffer = "! Stack definitions for OFT\n"
    #
    file_buffer = file_buffer + "INTEGER(i4), PARAMETER :: stack_len = " + str(stack_len) + "\n"
    #
    nmods = len(modules)
    file_buffer = file_buffer + "INTEGER(i4), PARAMETER :: stack_nmods = " + str(nmods) + "\n"
    file_buffer = file_buffer + "CHARACTER(LEN=stack_len), PARAMETER :: stack_mods(stack_nmods) = [ & \n"
    for i in range(nmods):
        mod_name = str(modules[i])
        if len(mod_name)>stack_len:
            print("ERROR: Module name exceeds allocated space!")
            print("MODULE: {0}".format(mod_name))
            print("LENGTH: {0} > {1}".format(len(mod_name),stack_len))
            sys.exit(1)

        if i<nmods-1:
            file_buffer = file_buffer + "'" + mod_name.ljust(stack_len) + "', & \n"
        else:
            file_buffer = file_buffer + "'" + mod_name.ljust(stack_len) + "'] \n"
    #
    nfuns = len(functions)
    file_buffer = file_buffer + "INTEGER(i4), PARAMETER :: stack_nfuns = " + str(nfuns) + "\n"
    file_buffer = file_buffer + "CHARACTER(LEN=stack_len), PARAMETER :: stack_funs(stack_nfuns) = [ & \n"
    for i in range(nfuns):
        fun_name = str(functions[i])
        if len(fun_name)>stack_len:
            print("ERROR: Function name exceeds allocated space!")
            print("FUNCTION: {0}".format(fun_name))
            print("LENGTH: {0} > {1}".format(len(fun_name),stack_len))
            sys.exit(1)

        if i<nfuns-1:
            file_buffer = file_buffer + "'" + fun_name.ljust(stack_len) + "', & \n"
        else:
            file_buffer = file_buffer + "'" + fun_name.ljust(stack_len) + "'] \n"
    #
    file_buffer = file_buffer + "INTEGER(i4), PARAMETER :: stack_fun_mods(stack_nfuns) = [ & \n"
    for i in range(nfuns):
        if i<nfuns-1:
            file_buffer = file_buffer + str(functions[i].module) + ", & \n"
        else:
            file_buffer = file_buffer + str(functions[i].module) + "] \n"
    #
    with open("include/stack_defs.h","w+") as fid:
        fid.write(file_buffer)
#----------------------------------------------------------------
# Remove STACK preprocessor definitions from a FORTRAN source file
#----------------------------------------------------------------
def clean_fortran_file(fid):
    file_buffer = ""
    for line in fid:
        if line.count("#define DEBUG_STACK_") > 0 or line.count("#undef DEBUG_STACK_") > 0:
            continue
        file_buffer = file_buffer + line
    return file_buffer
#----------------------------------------------------------------
# Script driver
#----------------------------------------------------------------
if __name__ == '__main__':
    test_run = False
    debug = False
    clean = False
    #
    import argparse
    parser = argparse.ArgumentParser(description='Generate in source stack and profile declarations for OFT.')
    parser.add_argument('-t', '--test_run', help='Run without writing any changes.', action="store_true")
    parser.add_argument('-d', '--debug', help='Display debug information.', action="store_true")
    parser.add_argument('-c', '--clean', help='Clean stack declarations from sources.', action="store_true")
    args=parser.parse_args()
    #
    if args.test_run:
        test_run = True
    if args.debug:
        debug = True
    if args.clean:
        clean = True
    # Check for correct run path
    correct_path = os.path.isfile("base/oft_local.F90")
    if not(correct_path):
        print("Invalid Run Directory!!!")
        print("Must be run from root source directory.")
        sys.exit(1)
    # Set source directories for search
    obj_dirs = ["base", "grid", "lin_alg", "fem", "physics", "."]
    stack_len = 40
    # Initialize containers
    skip_files = []
    modules = [module("dummy",1)]
    functions = [function("dummy",1,1)]
    # If clean is set remove all preprocessor defs
    if clean:
        print("\n%========================================%")
        print("Cleaning Source Files")
    # Otherwise parse source files and add STACK definitions
    else:
        print("\n%========================================%")
        print("Parsing Source Files")
    # Look in each source directory
    for dir in obj_dirs:
        files = os.listdir(dir)
        # Check each file in this directory
        for filename in files:
            # Check if file is a FORTRAN source file
            extension = pathlib.Path(filename).suffix
            # If file is a FORTRAN file parse and remove definitions
            if ".F90" == extension:
                path = os.path.join(dir,filename)
                print(path)
                with open(path,'r') as fid:
                    if clean:
                        new_file = clean_fortran_file(fid)
                    else:
                        modules, functions, new_file = parse_fortran_file(fid,modules,functions,debug)
                # If in test mode do not replace file
                if test_run:
                    continue
                with open(path,"w+") as fid:
                    fid.write(new_file)
    if not clean:
        # Print module statistics from search
        print("\n%========================================%")
        print("Found {0} Modules".format(len(modules)))
        if debug:
            for mod in modules:
                print(mod)
        # Print subroutine statistics from search
        print("\n%========================================%")
        print("Found {0} Functions".format(len(functions)))
        if debug:
            for func in functions:
                print(repr(func))
        # Create header file with STACK variables
        if not test_run:
            create_debug_list(modules,functions)
