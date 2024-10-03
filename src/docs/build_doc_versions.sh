#!/bin/bash

# Setup version tracking file and menu patch string
> doc_versions.txt
MENU_PATCH='menudata.children.push({text:"Versions",url:"index.html",children:['

# Get information about git environment (commit, release)
COMMIT_HASH=$(git rev-parse --short HEAD)
CURR_RELEASE=$(git describe --tags --abbrev=0)

# Build current commit
mkdir latest && cd latest
cmake -DOFT_BUILD_DOCS:BOOL=ON -DOFT_DOCS_ONLY:BOOL=ON -DOFT_PACKAGE_NIGHTLY:BOOL=ON ../../src
make docs
cd ..
echo "latest" >> doc_versions.txt
MENU_PATCH+='{text:"'$COMMIT_HASH'",url:"../latest/index.html"},'

# Build most recent release
git checkout $CURR_RELEASE
mkdir $CURR_RELEASE && cd $CURR_RELEASE
cmake -DOFT_BUILD_DOCS:BOOL=ON -DOFT_DOCS_ONLY:BOOL=ON -DOFT_PACKAGE_NIGHTLY:BOOL=OFF ../../src
make docs
cd ..
echo $CURR_RELEASE >> doc_versions.txt
MENU_PATCH+='{text:"'$CURR_RELEASE'",url:"'../$CURR_RELEASE'/index.html"},'

# Build previous releases (explicitly set for now)
for verTag in v1.0.0-beta3 v1.0.0-beta2 ;
do
    git checkout $verTag
    mkdir $verTag && cd $verTag
    cmake -DOFT_BUILD_DOCS:BOOL=ON -DOFT_DOCS_ONLY:BOOL=ON -DOFT_PACKAGE_NIGHTLY:BOOL=OFF ../../src
    make docs
    cd ..
    echo $verTag >> doc_versions.txt
    MENU_PATCH+='{text:"'$verTag'",url:"'../$verTag'/index.html"},'
done

# Reset to head for following steps
git checkout $COMMIT_HASH

# Create common output directory and pointer page
mkdir -p doc/html
echo '<html><head><meta http-equiv="Cache-Control" content="no-cache, no-store, must-revalidate" />' > doc/html/index.html
echo '<meta http-equiv="Pragma" content="no-cache" />' >> doc/html/index.html
echo '<meta http-equiv="Expires" content="0" />' >> doc/html/index.html
echo '<meta http-equiv="REFRESH" content="0;URL='$CURR_RELEASE'/index.html"></head></html>' >> doc/html/index.html

# Copy other versions to output directory and patch menu
while read verTag; do
  mv $verTag/doc/html doc/html/$verTag
  echo ${MENU_PATCH%?}"]})" >> doc/html/$verTag/menudata.js
done <doc_versions.txt