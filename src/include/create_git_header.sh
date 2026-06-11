#!/bin/bash
HAVE_GIT=$(git status 1>/dev/null 2>&1; echo $?)
if [[ "${HAVE_GIT}" -eq 0 ]]; then
  GITBRANCH=$(git branch | grep "*" | sed "s/* //" | sed "s/ /_/g" | sed -e "s/[()]//g" | sed 's/\r//g')
  GITVERSION=$(git rev-parse --short HEAD | sed 's/\r//g')
  GITTAG=$(git describe --tags --abbrev=0)
  IS_TAGGED="$(git tag --points-at HEAD | head -n 1)"
  if [[ -n "$IS_TAGGED" ]]; then
    IS_RELEASE=true
  else
    IS_RELEASE=false
  fi
else
  if [[ -z "$GITBRANCH" ]]; then
    GITBRANCH=unknown
  fi
  if [[ -z "$GITVERSION" ]]; then
    GITVERSION=unknown
  fi
  if [[ -z "$GITTAG" ]]; then
    GITTAG=unknown
  fi
  if [[ -z "$IS_RELEASE" ]]; then
    IS_RELEASE=false
  fi
fi
echo '#define GITTAG "'$GITTAG'"'
if [[ "$IS_RELEASE" != true ]]; then
  echo '#define GITBRANCH "'$GITBRANCH'"'
  echo '#define GITVERSION "'$GITVERSION'"'
fi
