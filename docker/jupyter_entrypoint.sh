#!/bin/bash

# Setup python environment
source $OFT_PYENV/bin/activate
export PYTHONPATH=$OFT_ROOTPATH/python:$PYTHONPATH

print_help() {
  echo "Usage: container_entrypoint.sh [options]
  
  -h Print this help message
  -j Run Jupyter lab server on port 8888
  -s Path to the script to run
  -n Path to the Jupyter notebook to execute
  -r Path to requirements.txt file to install dependencies before running the script or notebook"
  return 0  # Success
}

# Parse arguments
LAUNCH_JUPYTER=false
RUN_NOTEBOOK=false
while getopts ":hjs:n:r:" opt; do
  case $opt in
    h) print_help; exit 0;;
    j) LAUNCH_JUPYTER=true;;
    s) SCRIPT=$OPTARG;;
    n) SCRIPT=$OPTARG; RUN_NOTEBOOK=true;;
    r) REQUIREMENTS=$OPTARG;;
    \?) echo "Invalid option"; echo ""; print_help; exit 1;;
  esac
done

# Install dependencies if requirements.txt provided
if [[ -n "$REQUIREMENTS" ]]; then
  pip install --no-cache-dir -r "$REQUIREMENTS"
fi

# Execute desired function
if [ "$LAUNCH_JUPYTER" = true ]; then
  exec jupyter lab --no-browser --ip=''
else
  if [ "$RUN_NOTEBOOK" = true ]; then
    exec jupyter nbconvert --execute --to notebook --inplace --ExecutePreprocessor.kernel_name=python3 "$SCRIPT"
  else
    if [[ -z "$SCRIPT" ]]; then
      echo "No options provided"
      echo ""
      print_help
      exit 1
    fi
    # Run script passed as argument
    exec python "$SCRIPT"
  fi
fi
