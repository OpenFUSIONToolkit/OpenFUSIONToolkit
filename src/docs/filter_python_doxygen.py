#!/usr/bin/env python3
import re
import sys

# Get input file
if len(sys.argv) != 2:
    print("Usage: python filter_python_doxygen.py <input_file>")
    sys.exit(1)

debug = False
indent_level = -1
skip_next_line = False
start_doc = False
read_doc = False
doc_sentinal = None
doc_lines = []
var_def = None
with open(sys.argv[1], 'r') as fid:
    for line in fid:
        # Skip line if desired
        if skip_next_line:
            skip_next_line = False
            if debug:
                print(' '*(indent_level) + '#' + line[indent_level:], end='')
            continue
        # Start reading documentation for @property function
        if start_doc:
            name = re.search(r'^\s*def\s+(\w+)\s*\(', line).group(1)
            doc_lines.append('')
            var_def = ' '*(indent_level) + 'self.' + name + ' = None'
            if debug:
                print(' '*(indent_level) + '#' + line[indent_level:], end='')
            start_doc = False
            read_doc = True
            continue
        # Continue reading documentation for @property function
        if read_doc:
            if doc_sentinal is None:
                loc = line.find('"""!')
                if loc >= 0:
                    doc_sentinal = '"""'
                else:
                    loc = line.find("'''!")
                    if loc >= 0:
                        doc_sentinal = "'''"
                if doc_sentinal is None:
                    read_doc = False # Abandon doc parsing if no docstring is found
                    if debug:
                        print(' '*(indent_level) + '#' + line[indent_level:], end='')
                    doc_lines.append(var_def)
                else:
                    loc2 = line[loc+4:].find(doc_sentinal)
                    if loc2 >= 0:
                        doc_lines.append(' '*(indent_level) + '## ' + line[loc+4:loc+4+loc2].strip())
                        doc_lines.append(var_def)
                        read_doc = False
                        doc_sentinal = None
                    else:
                        doc_lines.append(' '*(indent_level) + '## ' + line[loc+4:].strip())
                continue
            else:
                sentintal_loc = line.find(doc_sentinal)
                if sentintal_loc >= 0:
                    read_doc = False
                    doc_sentinal = None
                    doc_lines.append(' '*(indent_level) + '# ' + line[indent_level:sentintal_loc].strip())
                    doc_lines.append(var_def)
                else:
                    doc_lines.append(' '*(indent_level) + '# ' + line[indent_level:].strip())
                continue
        # Detect @property, start doc parsing, and remove from output
        if re.search(r'^\s*@property', line):
            indent_level = line.find('@')
            start_doc = True
            if debug:
                print(' '*(indent_level) + '#' + line[indent_level:], end='')
            continue
        # Detect @*.setter and remove from output
        if re.search(r'^\s*@.+\.setter', line):
            indent_level = line.find('@')
            skip_next_line = True # Skip setter definitions
            if debug:
                print(' '*(indent_level) + '#' + line[indent_level:], end='')
            continue
        # Detect end of property definition and reset
        if indent_level > 0:
            if re.search(r'^\s*def\s+', line) or re.search(r'^\s*class\s+', line):
                if len(doc_lines) > 0:
                    print(' '*(indent_level) + 'def _zz_doxygen_dummy(self):')
                    print('\n'.join(['    ' + doc_line for doc_line in doc_lines]))
                    print()
                    doc_lines = []
                indent_level = -1
            else:
                if debug:
                    if len(line) <= 1:
                        print(' '*(indent_level) + '#')
                    else:
                        print(' '*(indent_level) + '#' + line[indent_level:], end='')
                continue
        # Print regular line
        print(line, end='')