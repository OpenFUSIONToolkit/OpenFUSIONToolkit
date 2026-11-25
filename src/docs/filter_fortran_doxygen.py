#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Doxygen input filter for Fortran sources that instantiate module body templates
# (`*.inc` files included with `#include` and parameterized by `#define` macros).
#
# - Fortran files: each `#include "name.inc"` of a template in the same directory is
#   replaced by the template, with the macros defined at that point substituted in
#   its code lines (preprocessor lines are kept, so Doxygen still evaluates `#ifdef`).
#   Every instantiated module is then documented with its procedures. Files without
#   such includes are passed through unchanged.
# - Template files (`*.inc`): only a file description is emitted, built from the
#   header comment, so the templated procedures are documented once per
#   instantiation and not again as part of the template.
#------------------------------------------------------------------------------
import os
import re
import sys

DIRECTIVE = re.compile(r'^\s*#')
DEFINE = re.compile(r'^\s*#\s*define\s+(\w+)(\(([^)]*)\))?\s*(.*)$')
UNDEF = re.compile(r'^\s*#\s*undef\s+(\w+)')
INCLUDE = re.compile(r'^\s*#\s*include\s+"([^"]+\.inc)"')


def expand_call(line, name, args, body):
    '''Expand calls of a function-like macro (arguments may contain parentheses)'''
    pattern = re.compile(r'\b' + name + r'\s*\(')
    out = ''
    while True:
        match = pattern.search(line)
        if match is None:
            return out + line
        depth = 1
        pos = match.end()
        values = []
        start = pos
        while pos < len(line) and depth > 0:
            if line[pos] == '(':
                depth += 1
            elif line[pos] == ')':
                depth -= 1
                if depth == 0:
                    values.append(line[start:pos].strip())
            elif line[pos] == ',' and depth == 1:
                values.append(line[start:pos].strip())
                start = pos + 1
            pos += 1
        if depth > 0:  # Unbalanced (eg. continued line), leave unchanged
            return out + line
        result = body
        for arg, value in zip(args, values):
            result = re.sub(r'\b' + arg + r'\b', lambda m, value=value: value, result)
        out += line[:match.start()] + result
        line = line[pos:]


def substitute(line, macros):
    '''Substitute object-like and function-like macros in a code line'''
    for name, (args, body) in macros.items():
        if args is None:
            line = re.sub(r'\b' + name + r'\b', lambda m, body=body: body, line)
        else:
            line = expand_call(line, name, args, body)
    return line


class Expander:
    '''Expand template includes, evaluating conditionals on the template macros

    Template macro definitions are removed from the output (Doxygen does not parse
    the following Fortran correctly when they are present) and `#ifdef`/`#ifndef`
    on them are evaluated here. Other preprocessor lines are passed to Doxygen.
    '''
    def __init__(self, names):
        self.names = names  # Macros defined in the file (template parameters)
        self.macros = {}
        self.stack = []  # (handled, active) for each open conditional

    def active(self):
        return all(active for _, active in self.stack)

    def directive(self, line):
        '''Process a preprocessor line, returning the text to emit'''
        match = DEFINE.match(line)
        if match:
            if self.active():
                args = None
                if match.group(2) is not None:
                    args = [a.strip() for a in match.group(3).split(',') if a.strip()]
                self.macros[match.group(1)] = (args, match.group(4).strip())
            return '\n'
        match = UNDEF.match(line)
        if match:
            if self.active():
                self.macros.pop(match.group(1), None)
            return '\n'
        match = re.match(r'^\s*#\s*(ifdef|ifndef)\s+(\w+)', line)
        if match and match.group(2) in self.names:
            defined = match.group(2) in self.macros
            self.stack.append((True, defined if match.group(1) == 'ifdef' else not defined))
            return '\n'
        if re.match(r'^\s*#\s*if', line):
            self.stack.append((False, True))
        elif re.match(r'^\s*#\s*else', line) and self.stack and self.stack[-1][0]:
            self.stack[-1] = (True, not self.stack[-1][1])
            return '\n'
        elif re.match(r'^\s*#\s*endif', line) and self.stack:
            handled, _ = self.stack.pop()
            if handled:
                return '\n'
            return line if self.active() else '\n'
        return line if self.active() else '\n'

    def code(self, line):
        '''Process a Fortran line, returning the text to emit'''
        if not self.active():
            return '\n'
        if line.lstrip().startswith('!'):
            # An apostrophe in a comment is taken as an unterminated character literal by
            # Doxygen's preprocessor, which then drops the following directives and code
            return substitute(line, self.macros).replace("'", '\u2019')
        return substitute(line, self.macros)


def filter_fortran(path):
    '''Expand template includes of a Fortran file with the macros in effect'''
    folder = os.path.dirname(os.path.abspath(path))
    with open(path, 'r') as fid:
        lines = fid.readlines()
    templates = [m.group(1) for m in (INCLUDE.match(line) for line in lines)
                 if m and os.path.isfile(os.path.join(folder, m.group(1)))]
    if len(templates) == 0:
        sys.stdout.write(''.join(lines))
        return
    expander = Expander({m.group(1) for m in (DEFINE.match(line) for line in lines) if m})
    for i, line in enumerate(lines):
        module = re.match(r'^\s*MODULE\s+(\w+)\s*$', line, re.IGNORECASE)
        if module and i+1 < len(lines):
            match = INCLUDE.match(lines[i+1])
            if match and match.group(1) in templates:
                header = template_header(os.path.join(folder, match.group(1)))
                sys.stdout.write(module_doc(module.group(1), match.group(1), header, expander.macros))
        match = INCLUDE.match(line)
        if match and match.group(1) in templates:
            with open(os.path.join(folder, match.group(1)), 'r') as fid:
                for tline in fid:
                    if DIRECTIVE.match(tline):
                        sys.stdout.write(expander.directive(tline))
                    elif tline.lstrip().startswith('!') and not tline.lstrip().startswith(('!>', '!!', '!<', '!$')):
                        sys.stdout.write('\n')  # Plain comments (template header) are not needed here
                    else:
                        sys.stdout.write(expander.code(tline))
            continue
        if DIRECTIVE.match(line):
            sys.stdout.write(expander.directive(line))
        else:
            sys.stdout.write(expander.code(line))


def template_header(path):
    '''Doxygen comment lines (without the leading markers) from the header of a template'''
    header = []
    with open(path, 'r') as fid:
        for line in fid:
            stripped = line.strip()
            if not stripped.startswith('!'):
                break
            text = stripped.lstrip('!').rstrip()
            if text.startswith('---') or 'SPDX-License-Identifier' in text or 'Open FUSION Toolkit' in text:
                continue
            header.append(text)
    while header and header[0].strip() == '':
        header.pop(0)
    out = []
    in_list = False
    for text in header:
        # Macro lists ("NAME - description") are rendered as a list
        item = re.match(r'^\s+(\w+)\s+-\s+(.*)$', text)
        if item:
            out.append('- `{0}`: {1}'.format(item.group(1), item.group(2)))
            in_list = True
            continue
        if in_list and text.startswith('    '):
            out.append('  ' + text.strip())
            continue
        in_list = False
        out.append(text.strip().replace("'", '\u2019'))
    return out


def filter_template(path):
    '''Emit a file description for a template from its header comment'''
    sys.stdout.write('!> @file {0}\n!!\n'.format(os.path.basename(path)))
    for text in template_header(path):
        sys.stdout.write('!! {0}\n'.format(text))
    sys.stdout.write('!!\n!! The procedures are documented in each module instantiating this template.\n')


def module_doc(name, template, header, macros):
    '''Documentation block for a module instantiated from a template'''
    out = ['!> Instance of the `{0}` template'.format(template), '!!']
    values = ['`{0}` = `{1}`'.format(k, v[1] if v[0] is None else '{0}({1})'.format(k, ','.join(v[0])) + ' -> ' + v[1])
              for k, v in macros.items()]
    if values:
        out.append('!! Template parameters: ' + ', '.join(values))
        out.append('!!')
    out += ['!! ' + text for text in header]
    return '\n'.join(out) + '\n'


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python filter_fortran_doxygen.py <input_file>')
        sys.exit(1)
    if sys.argv[1].endswith('.inc'):
        filter_template(sys.argv[1])
    else:
        filter_fortran(sys.argv[1])
