#!/usr/bin/env bash

# env control arguments
BUILD_DIR="${build_dir:-build}"
RUN_TESTS="${run_tests:-0}"
BUILD_TESTS="${build_tests:-0}"

# build directory
mkdir "${BUILD_DIR}"

cd "${BUILD_DIR}" || {
  echo 'entry into build folder failed'
  exit 1
}

# configuration step
cmake -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_CXX_STANDARD=14 \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DTUDAT_BUILD_TESTS="${BUILD_TESTS}" \
  ..

# build step
cmake --build .
