#!/bin/bash
HAVE_GIT=$(git status 1>/dev/null 2>&1; echo $?)
if [[ "${HAVE_GIT}" -eq 0 ]]
then
  GITBRANCH=$(git branch | grep "*" | sed "s/* //" | sed "s/ /_/g" | sed -e "s/[()]//g" | sed 's/\r//g')
  GITVERSION=$(git rev-parse --short HEAD | sed 's/\r//g')
else
  GITBRANCH=unknown
  GITVERSION=unknown
fi
echo "#define GITBRANCH '${GITBRANCH}'"
echo "#define GITVERSION '${GITVERSION}'"