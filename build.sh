#!/bin/bash

rm -r build
mkdir build; cd build 
cmake ../cmake
PKGS="-D PKG_GRANULAR=on -D PKG_VTK=on -D PKG_SRD=on -D PKG_MPI=on"
cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake $PKGS ../cmake
#cmake -D CMAKE_BUILD_TYPE=Debug -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake $PKGS ../cmake
cmake --build . -- -j 10
