#!/usr/bin/env python
#---------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#---------------------------------------------------------------------------
#
# Script for handling creation of Doxygen documentation from structured examples.
#
# Note: Script must be run from the root "src" directory.
#
#---------------------------------------------------------------------------
from __future__ import print_function
import sys
import os
import glob
import shutil
import subprocess
import re
#----------------------------------------------------------------
# Comment seperator template
#----------------------------------------------------------------
sep_template = """!---------------------------------------------------------------------------
! {0}
!---------------------------------------------------------------------------
"""
output_reg = re.compile(r'    [ ]*[\S]+')

def run_command(command, timeout=10):
    # Run shell command
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Wait for process to complete or timeout
    try:
        outs, _ = pid.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        pid.kill()
        outs, _ = pid.communicate()
        print("WARNING: Command timeout")
    errcode = pid.poll()
    result = outs.decode()
    return result, errcode

#----------------------------------------------------------------
# Parse structured example file and generate documentation
#----------------------------------------------------------------
def parse_fortran_file(fid):
    # Initialize local variables
    file_buffer = ""
    code_buffer = ""
    doc_buffer = ""
    full_code = ""
    read_full = False
    incode = False
    indoc = True
    include_full = True
    line_number = 0
    entry_count = 0
    mod_id = 1
    for line in fid:
        line_number = line_number + 1
        # Determine documentation prefix
        if line_number == 1:
            doc_buffer = line[2:]
            tmp = line.split('{')[1]
            tmp = tmp.split('}')[0]
            ex_prefix = tmp[1:]
            continue
        # Check if full source version should be included
        if line.find("! START SOURCE") == 0:
          read_full = True
          continue
        if line.find("! STOP SOURCE") == 0:
          read_full = False
          continue
        # Found documentation block
        if line.find("!!") == 0:
            # Start documentation section
            if incode:
                file_buffer = file_buffer + code_buffer
                file_buffer = file_buffer + "~~~~~~~~~\n\n"
                incode = False
                indoc = True
                doc_buffer = ""
            # Add to existing documentation section
            doc_buffer = doc_buffer + line[2:]
            if read_full:
              # Search for section break
              if line.find(r'\subsection') >= 0:
                splits = line.split(" ")
                term_found = False
                sec_title = None
                for (i,term) in enumerate(splits):
                  if(term_found):
                    sec_title = ' '.join(splits[i+1:])
                    break
                  if term.find(r'\subsection') >= 0:
                    term_found = True
                if sec_title != None:
                  full_code = full_code + sep_template.format(sec_title[:-1])
        # Code block
        else:
            # Start code section
            if indoc:
                file_buffer = file_buffer + doc_buffer + "\n"
                incode = True
                indoc = False
                code_buffer = "~~~~~~~~~{.F90}\n"
            # Add to existing code section
            code_buffer = code_buffer + line
            if read_full:
              full_code = full_code + line
    # Terminate current section
    if incode:
        file_buffer = file_buffer + code_buffer
        file_buffer = file_buffer + "~~~~~~~~~\n"
    elif indoc:
        file_buffer = file_buffer + doc_buffer + "\n"
    # Add full source if needed
    if full_code != "":
        file_buffer = file_buffer + "\n" + r'\section ' + ex_prefix + "_full Complete Source\n"
        file_buffer = file_buffer + "~~~~~~~~~{.F90}\n"
        file_buffer = file_buffer + full_code
        file_buffer = file_buffer + "~~~~~~~~~\n"
    return file_buffer
#----------------------------------------------------------------
# Script driver
#----------------------------------------------------------------
if __name__ == '__main__':
    # Check for correct run path
    correct_path = os.path.isfile("base/oft_local.F90")
    if not(correct_path):
        print("Invalid Run Directory!!!")
        print("Must be run from root source directory.")
        sys.exit(1)
    # Create docs folder if necessary
    if not os.path.isdir("docs/generated"):
        os.mkdir("docs/generated")
    if not os.path.isdir("docs/generated/images"):
        os.mkdir("docs/generated/images")
    # Loop over all example files
    files = glob.glob("examples/*.F90") + glob.glob("examples/*/*/*.F90")
    print("\n==========================================")
    print("Parsing Example Files")
    for filename in files:
        basename = os.path.basename(filename)
        print(filename,basename)
        with open(filename,'r') as fid:
            new_file = parse_fortran_file(fid)
        # Write documentation file to doc folder
        path = "docs/generated/doc_" + basename + ".md"
        with open(path,"w+") as fid:
            fid.write(new_file)
    # Add Jupyter notebooks to documentation
    print()
    print("\n==========================================")
    print("Converting Jupyter notebooks")
    eq_reg = re.compile(r'\$(.*?)\$')
    files = glob.glob("examples/*/*/*.ipynb")
    for filename in files:
        print(filename)
        base_path = filename.split('.')[0]
        full_name = os.path.basename(filename)
        file_name, _ = os.path.splitext(full_name)
        _, errcode = run_command("jupyter nbconvert --to markdown {0}".format(filename))
        if errcode != 0:
           print("Jupyter notebook->markdown conversion failed for {0}".format(filename))
           continue
        # Copy files to doc directory
        with open(base_path+".md", 'r') as fid:
           contents = fid.read()
        # Update image paths
        contents = contents.replace("{0}_files".format(file_name), "images")
        contents = contents.replace("[png]", "[]")
        # Convert notes to note blocks
        contents = contents.replace("**Note:**", r'@note')
        contents = contents.replace("**Warning:**", r'@warning')
        # Convert code block style
        contents_split = contents.split('```')
        for i, content_segment in enumerate(contents_split):
            if (i % 2) == 0:
                contents_split[i] = re.sub(eq_reg,r'\\f$\1\\f$',content_segment)
        contents = '```'.join(contents_split)
        contents = contents.replace('```python','~~~~~~~~~~~~~{.py}') 
        contents = contents.replace('```','~~~~~~~~~~~~~')
        # Write updated markdown file to doc directory
        with open("docs/generated/doc_{0}.md".format(file_name), 'w+') as fid:
           fid.write(contents)
        # Copy images to img directory
        imgs = glob.glob("{0}_files/*".format(base_path))
        for img in imgs:
           shutil.copy(img, "docs/generated/images/{0}".format(os.path.basename(img)))
        # Remove temporary files
        if os.path.isdir("{0}_files".format(base_path)):
            shutil.rmtree("{0}_files".format(base_path))
        os.remove(base_path+".md")
