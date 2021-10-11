#!/usr/bin/env bash

cd ../..

#mkdir install
#cd install
#export TUDATPY_INSTALL_DIR=`pwd`
#cd ..

mkdir build
cd build
export TUDATPY_BUILD_DIR=`pwd`

cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$TUDATPY_INSTALL_DIR \
    -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
    -DCMAKE_CXX_STANDARD=14 \
    -DTUDATPY_CONDA_BUILD=on \
    -D_GLIBCXX_USE_CXX11_ABI=1 \
    ..

make -j4
