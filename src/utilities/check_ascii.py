#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Lint check for non-ASCII characters in tracked source and Markdown files.
#
# Source files (Fortran, C/C++, Python, CMake, shell, and "*.in" templates)
# containing non-ASCII characters are reported as errors and cause a non-zero
# exit status. Markdown files are reported as warnings only.
#
# Files are found with "git ls-files", so the script can be run from anywhere
# inside the repository. When run under GitHub Actions, results are also
# emitted as workflow annotations.
#
#------------------------------------------------------------------------------
import os
import subprocess
import sys

# Extensions (lowercase) treated as source files
SOURCE_EXTS = (
    '.f90', '.f', '.f03', '.f08', '.inc',
    '.c', '.h', '.cxx', '.cpp', '.cc', '.hpp', '.hxx',
    '.py', '.cmake', '.sh', '.in',
)
# File names treated as source files regardless of extension
SOURCE_NAMES = ('CMakeLists.txt',)
# Extensions (lowercase) treated as Markdown files
MARKDOWN_EXTS = ('.md',)


def classify(path):
    '''Return "error" for source files, "warning" for Markdown, or None to skip'''
    name = os.path.basename(path)
    ext = os.path.splitext(name)[1].lower()
    if (name in SOURCE_NAMES) or (ext in SOURCE_EXTS):
        return 'error'
    if ext in MARKDOWN_EXTS:
        return 'warning'
    return None


def find_non_ascii(path):
    '''Return a list of (line, column, character) for each non-ASCII character in a file'''
    with open(path, 'rb') as fid:
        raw_lines = fid.read().splitlines()
    hits = []
    for i, raw_line in enumerate(raw_lines):
        if max(raw_line, default=0) < 128:
            continue
        line = raw_line.decode('utf-8', errors='replace')
        for j, char in enumerate(line):
            if ord(char) > 127:
                hits.append((i+1, j+1, char))
    return hits


def describe(char):
    '''Printable description of a character, eg. "'²' (U+00B2)"'''
    if char == '�':
        return 'invalid UTF-8 byte'
    return "'{0}' (U+{1:04X})".format(char, ord(char))


def report(level, path, line, col, chars, annotate):
    message = 'Non-ASCII character(s) found: {0}'.format(', '.join(describe(char) for char in chars))
    if annotate:
        print('::{0} file={1},line={2},col={3}::{4}'.format(level, path, line, col, message))
    else:
        print('{0}:{1}:{2}: {3}: {4}'.format(path, line, col, level, message))


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Check tracked source and Markdown files for non-ASCII characters.')
    parser.add_argument('files', nargs='*', help='Files to check (default: all files tracked by git)')
    args = parser.parse_args()

    # Work from the repository root so paths match those in annotations
    root = subprocess.check_output(['git', 'rev-parse', '--show-toplevel']).decode().strip()
    if len(args.files) > 0:
        files = [os.path.relpath(os.path.abspath(f), root) for f in args.files]
    else:
        files = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode('utf-8').split('\0')
    os.chdir(root)
    annotate = (os.environ.get('GITHUB_ACTIONS', '') == 'true')

    counts = {'error': 0, 'warning': 0}
    for path in sorted(files):
        level = classify(path)
        if (level is None) or (not os.path.isfile(path)):
            continue
        # Group characters by line so each line is reported once
        by_line = {}
        for line, col, char in find_non_ascii(path):
            if line not in by_line:
                by_line[line] = [col, []]
            if char not in by_line[line][1]:
                by_line[line][1].append(char)
        for line, (col, chars) in by_line.items():
            report(level, path, line, col, chars, annotate)
            counts[level] += 1

    print('')
    print('Non-ASCII check: {0} source line(s) with errors, {1} Markdown line(s) with warnings'.format(
        counts['error'], counts['warning']))
    if counts['error'] > 0:
        print('Source files must contain only ASCII characters')
        sys.exit(1)


if __name__ == '__main__':
    main()
