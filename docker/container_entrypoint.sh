#!/bin/bash
# Setup python environment
source /usr/local/OFT_venv/bin/activate
export PYTHONPATH=$OFT_ROOTPATH/python:$PYTHONPATH

# Parse arguments
LAUNCH_JUPYTER=true
RUN_NOTEBOOK=false
while getopts "hsj:" opt; do
  case $opt in
    h) echo "Usage: script -h (help) or -s (script)"; exit 0;;
    s) SCRIPT=$OPTARG; LAUNCH_JUPYTER=false;;
    j) SCRIPT=$OPTARG; LAUNCH_JUPYTER=false; RUN_NOTEBOOK=true;;
    \?) echo "Invalid option\n\nUsage: script -h (help) or -s (script)"; exit 1;;
  esac
done

if [ "$LAUNCH_JUPYTER" = true ]; then
  jupyter lab --no-browser --ip=''
else
  if [ "$RUN_NOTEBOOK" = true ]; then
    jupyter nbconvert --execute --to notebook --inplace --ExecutePreprocessor.kernel_name=Python3 "$SCRIPT"
  else
    # Run script passed as argument
    python "$SCRIPT"
  fi
fi
