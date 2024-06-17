#!/bin/bash

# Setup version tracking file and menu patch string
> doc_versions.txt
MENU_PATCH='menudata.children.push({text:"Versions",url:"index.html",children:['

# Get information about git environment (commit, releases, etc.)
GIT_DESC=$(git describe --tags --match 'v*')
CURR_RELEASE=$(git describe --tags --abbrev=0)
RESET_HASH=$(git rev-parse HEAD)

# Build current commit if not a release commit
if [[ "$GIT_DESC" =~ -([0-9]+)-([a-z0-9]{8}) ]]
then
    COMMIT_HASH=$(git rev-parse --short HEAD)
    mkdir $COMMIT_HASH && cd $COMMIT_HASH

    # Build current version
    cmake -DOFT_BUILD_DOCS:BOOL=ON -DOFT_DOCS_ONLY:BOOL=ON -DOFT_PACKAGE_NIGHTLY:BOOL=ON ../../src
    make docs
    cd ..
    echo $COMMIT_HASH >> doc_versions.txt
    MENU_PATCH+='{text:"'$COMMIT_HASH'",url:"'../$COMMIT_HASH'/index.html"},'

    # Checkout most recent release
    git checkout $CURR_RELEASE
fi

# Build most recent release
mkdir $CURR_RELEASE && cd $CURR_RELEASE
cmake -DOFT_BUILD_DOCS:BOOL=ON -DOFT_DOCS_ONLY:BOOL=ON -DOFT_PACKAGE_NIGHTLY:BOOL=OFF ../../src
make docs
cd ..
echo $CURR_RELEASE >> doc_versions.txt
MENU_PATCH+='{text:"'$CURR_RELEASE'",url:"'../$CURR_RELEASE'/index.html"},'

# Build previous releases (explicitly set for now)
for verTag in v1.0.0-beta2 v1.0.0-beta3 ;
do
    git checkout $verTag
    mkdir $verTag && cd $verTag
    cmake -DOFT_BUILD_DOCS:BOOL=ON -DOFT_DOCS_ONLY:BOOL=ON -DOFT_PACKAGE_NIGHTLY:BOOL=OFF ../../src
    make docs
    cd ..
    echo $verTag >> doc_versions.txt
    MENU_PATCH+='{text:"'$verTag'",url:"'../$verTag'/index.html"},'
done

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