#!/bin/bash

print_help() {
  echo "Usage: build_container.sh [options]
  
  -h Print this help message
  -d Use docker instead of podman
  -n Force no-cache build
  -f Path to the Dockerfile
  -t Image tag"
  return 0  # Success
}

# Parse arguments
USE_DOCKER=false
USE_CACHE=true
while getopts ":hdnf:t:" opt; do
  case $opt in
      h) print_help; exit 0;;
      d) USE_DOCKER=true;;
      n) USE_CACHE=false;;
      f) FILE_PATH=$OPTARG;;
      t) IMAGE_TAG=$OPTARG;;
      \?) echo "Invalid option"; echo ""; print_help; exit 1;;
  esac
done

if [[ -z "$FILE_PATH" || -z "$IMAGE_TAG" ]]; then
    echo "Error: -f and -t are required."
    echo "Usage: $0 -f <Dockerfile_path> -t <image_tag>"
    exit 1
fi

HAVE_GIT=$(git status 1>/dev/null 2>&1; echo $?)
if [[ "${HAVE_GIT}" -eq 0 ]]; then
  GITBRANCH=$(git branch | grep "*" | sed "s/* //" | sed "s/ /_/g" | sed -e "s/[()]//g" | sed 's/\r//g')
  GITVERSION=$(git rev-parse --short HEAD | sed 's/\r//g')
else
  if [[ -z "$GITBRANCH" ]]; then
    GITBRANCH=unknown
  fi
  if [[ -z "$GITVERSION" ]]; then
    GITVERSION=unknown
  fi
fi

HAVE_PODMAN=$(podman --version 1>/dev/null 2>&1; echo $?)
if [[ "${HAVE_PODMAN}" -ne 0 ]]; then
  if [ "$USE_DOCKER" = false ]; then
    echo "podman not available, falling back to docker"
  fi
  USE_DOCKER=true
fi
if [ "$USE_DOCKER" = true ]; then
  HAVE_DOCKER=$(docker --version 1>/dev/null 2>&1; echo $?)
  if [[ "${HAVE_DOCKER}" -ne 0 ]]; then
    echo "docker not available, exiting"
    exit 1
  fi
  if [ "$USE_CACHE" = true ]; then
    docker build --build-arg GITBRANCH=$GITBRANCH --build-arg GITVERSION=$GITVERSION -f $FILE_PATH -t $IMAGE_TAG .
  else
    docker build --no-cache --build-arg GITBRANCH=$GITBRANCH --build-arg GITVERSION=$GITVERSION -f $FILE_PATH -t $IMAGE_TAG .
  fi
else
  if [ "$USE_CACHE" = true ]; then
    podman build --build-arg GITBRANCH=$GITBRANCH --build-arg GITVERSION=$GITVERSION -f $FILE_PATH -t $IMAGE_TAG .
  else
    podman build --no-cache --build-arg GITBRANCH=$GITBRANCH --build-arg GITVERSION=$GITVERSION -f $FILE_PATH -t $IMAGE_TAG .
  fi
fi
