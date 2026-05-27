#!/bin/bash

# Setup python environment
source $OFT_PYENV/bin/activate
export PYTHONPATH=$OFT_ROOTPATH/python:$PYTHONPATH

print_help() {
  echo "Usage: container_entrypoint.sh [options]
  
  -h Print this help message
  -s Path to the script to run
  -r Path to requirements.txt file to install dependencies before running the script or notebook"
  return 0  # Success
}

# Parse arguments
LAUNCH_JUPYTER=false
RUN_NOTEBOOK=false
while getopts ":hs:r:" opt; do
  case $opt in
    h) print_help; exit 0;;
    s) SCRIPT=$OPTARG;;
    r) REQUIREMENTS=$OPTARG;;
    \?) echo "Invalid option"; echo ""; print_help; exit 1;;
  esac
done

# Install dependencies if requirements.txt provided
if [[ -n "$REQUIREMENTS" ]]; then
  pip install --no-cache-dir -r "$REQUIREMENTS"
fi

# Execute desired function
if [[ -z "$SCRIPT" ]]; then
  echo "No options provided"
  echo ""
  print_help
  exit 1
fi
# Run script passed as argument
exec python "$SCRIPT"
