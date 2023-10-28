#!/bin/bash
HAVE_GIT=$(git status 1>/dev/null 2>&1; echo $?)
if [[ "${HAVE_GIT}" -eq 0 ]]
then
  GITBRANCH=$(git branch | grep "*" | sed "s/* //" | sed "s/ /_/g" | sed -e "s/[()]//g")
  GITVERSION=$(git rev-parse --short HEAD)
else
  GITBRANCH=unknown
  GITVERSION=unknown
fi
echo "#define GITBRANCH '${GITBRANCH}'"
echo "#define GITVERSION '${GITVERSION}'"