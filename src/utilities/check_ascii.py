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
# Only files under SOURCE_DIR are checked, excluding those under EXCLUDE_DIRS.
# Source files (Fortran, C/C++, and Python) containing non-ASCII characters
# are reported as errors and cause a non-zero exit status. Markdown files
# under MARKDOWN_DIR are reported as warnings only; other Markdown files are
# not checked.
#
# Files are found with "git ls-files", so the script can be run from anywhere
# inside the repository. When run under GitHub Actions, errors are also
# emitted as workflow annotations, along with a single summary warning for
# Markdown files.
#
#------------------------------------------------------------------------------
import os
import subprocess
import sys

# Extensions (lowercase) treated as source files
SOURCE_EXTS = (
    '.f90', '.f', '.f03', '.f08',
    '.c', '.h', '.cxx', '.cpp', '.cc', '.hpp', '.hxx',
    '.py',
)
# Extensions (lowercase) treated as Markdown files
MARKDOWN_EXTS = ('.md',)
# Directories below are relative to the repository root
# Directory containing all checked files
SOURCE_DIR = 'src'
# Directory containing the checked Markdown files
MARKDOWN_DIR = 'src/docs'
# Directories that are not checked
EXCLUDE_DIRS = ('src/utilities', 'src/ext_libs')


def in_dir(path, directory):
    '''Return True if a repository-relative path is inside a directory'''
    return path.startswith(directory + '/')


def classify(path):
    '''Return "error" for source files, "warning" for Markdown, or None to skip'''
    path = path.replace(os.sep, '/')
    if not in_dir(path, SOURCE_DIR):
        return None
    for exclude_dir in EXCLUDE_DIRS:
        if in_dir(path, exclude_dir):
            return None
    ext = os.path.splitext(path)[1].lower()
    if ext in SOURCE_EXTS:
        return 'error'
    if (ext in MARKDOWN_EXTS) and in_dir(path, MARKDOWN_DIR):
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
        # Invalid bytes are kept as U+DC80-U+DCFF so describe() can report them
        line = raw_line.decode('utf-8', errors='surrogateescape')
        for j, char in enumerate(line):
            if ord(char) > 127:
                hits.append((i+1, j+1, char))
    return hits


def describe(char):
    '''Printable description of a character, eg. "'<char>' (U+00B2)" or "invalid UTF-8 byte 0xE9"'''
    if 0xDC80 <= ord(char) <= 0xDCFF:
        return 'invalid UTF-8 byte 0x{0:02X}'.format(ord(char) - 0xDC00)
    return "'{0}' (U+{1:04X})".format(char, ord(char))


def escape_data(value):
    '''Escape a GitHub Actions workflow command message'''
    return value.replace('%', '%25').replace('\r', '%0D').replace('\n', '%0A')


def escape_property(value):
    '''Escape a GitHub Actions workflow command property value'''
    return escape_data(value).replace(':', '%3A').replace(',', '%2C')


def report(level, path, line, col, chars, annotate):
    message = 'Non-ASCII character(s) found: {0}'.format(', '.join(describe(char) for char in chars))
    print('{0}:{1}:{2}: {3}: {4}'.format(path, line, col, level, message))
    # Only errors are annotated; GitHub shows few annotations per step and
    # Markdown warnings would crowd out the errors
    if annotate and (level == 'error'):
        print('::error file={0},line={1},col={2}::{3}'.format(
            escape_property(path), line, col, escape_data(message)))


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

    if annotate and (counts['warning'] > 0):
        print('::warning::{0} line(s) in Markdown files contain non-ASCII characters (see log for details)'.format(
            counts['warning']))
    print('')
    print('Non-ASCII check: {0} source line(s) with errors, {1} Markdown line(s) with warnings'.format(
        counts['error'], counts['warning']))
    if counts['error'] > 0:
        print('Source files must contain only ASCII characters')
        sys.exit(1)


if __name__ == '__main__':
    main()
