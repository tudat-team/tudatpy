#!/usr/bin/env bash

# allow input to set compilation config
while getopts d:j:c:t flag
do
    case "${flag}" in
        d) build_dir=${OPTARG};;
        j) number_of_processors=${OPTARG};;
        c) clean_build=${OPTARG};;
        t) build_type=${OPTARG};;
    esac
done

# env control arguments
BUILD_DIR="${build_dir:-build}"
RUN_TESTS="${run_tests:-1}"
BUILD_TESTS="${build_tests:-1}"
NUMBER_OF_PROCESSORS=${number_of_processors:-1}
CLEAN_BUILD=${clean_build:-0}
BUILD_TYPE=${build_type:-Release}

# build directory
mkdir "${BUILD_DIR}"

cd "${BUILD_DIR}" || {
  echo 'entry into build folder failed'
  exit 1
}

# configuration step
cmake -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
  -DCMAKE_CXX_STANDARD=14 \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DTUDAT_BUILD_TESTS="${BUILD_TESTS}" \
  ..

# if required by user, clean
if [ "${CLEAN_BUILD}" = "true" ]; then
    cmake --build . --target clean -j"${NUMBER_OF_PROCESSORS}"
    printf "\nClean finished\n\n"
fi
# build
cmake --build . -j"${NUMBER_OF_PROCESSORS}"
