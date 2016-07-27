#!/bin/bash
for dir in ./*/
do
    dir=${dir%*/}
    cd ${dir##*/}
    if [ -f "test_main.cpp" ];  
    then
        make examples
    fi
    if [ -f "test2d" ];
    then
        echo "************Running "test2d"************"
        ./test2d input2d
    fi
    if [ -f "test3d" ];
    then
        echo "************Running "test3d"************"
        ./test2d input3d
    fi
    cd ..
done
