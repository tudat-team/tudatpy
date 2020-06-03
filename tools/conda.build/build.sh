#!/usr/bin/env bash

mkdir build
cd build

if [[ "$(uname)" == "Darwin" ]]; then
    export ENABLE_TESTS=no
else
    LDFLAGS="-lrt ${LDFLAGS}"
    export ENABLE_TESTS=yes
fi

cmake \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_PREFIX_PATH=$PREFIX \
    -DPREFIX=$PREFIX \
    -DTUDAT_CONDA_BUILD=on \
    -DUSE_TUDAT_EXAMPLE_APPLICATIONS=off \
    ..

make -j2

make install
