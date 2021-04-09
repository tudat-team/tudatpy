#!/usr/bin/env bash

# build directory
mkdir build

cd build || {
  echo 'entry into build folder failed'
  exit 1
}

# configuration step
cmake -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_CXX_STANDARD=14 \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  ..

# build step
cmake --build .
