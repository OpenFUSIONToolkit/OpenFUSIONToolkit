#!/bin/bash

ROOT_PATH=$(pwd)
rm -f $ROOT_PATH/example_test.err
touch $ROOT_PATH/example_test.err

cd TokaMaker
./test_examples.sh $ROOT_PATH/example_test.err
cd ..

cd ThinCurr
./test_examples.sh $ROOT_PATH/example_test.err
cd ..

cd Marklin
./test_examples.sh $ROOT_PATH/example_test.err
cd ..

if [ $(wc -c < "$ROOT_PATH/example_test.err") -gt 0 ]; then
    cat $ROOT_PATH/example_test.err
    exit 1
else
    exit 0
fi
