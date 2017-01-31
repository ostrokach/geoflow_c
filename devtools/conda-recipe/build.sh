#!/bin/bash

mkdir build && cd build
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -DCMAKE_PREFIX_PATH:PATH=$PREFIX \
    -DCMAKE_INCLUDE_PATH:PATH="$PREFIX/include" \
    -DCMAKE_LIBRARY_PATH:PATH="$PREFIX/lib" \
    -DBUILD_SHARED_LIBS=ON \
    -DENABLE_GEOFLOW_APBS=ON \
    ..
make
make install
