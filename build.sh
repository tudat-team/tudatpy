#!/usr/bin/env bash

##################################################################
##################################################################
# This script is used to build the TU Delft Astrodynamics Toolbox
# Build options
#    -d: build directory
#    -j: number of processors
#    -t: build type (Release, Debug)
#    -c: clean build
##################################################################
##################################################################

# allow input to set compilation config
while getopts d:j:t:c flag
do
    case "${flag}" in
        d) build_dir=${OPTARG};;
        j) number_of_processors=${OPTARG};;
        t) build_type=${OPTARG};;
        c) 
            echo -en "\n\n\e[33m[USER INPUT REQUIRED] You have requested a clean build. Are you sure? [y/N]\e[0m "
            read -r response
            if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
            then
                clean_build=1
                echo -e "\n\n\e[33m   ---> Proceeding with a clean build!\e[0m\n\n"
            else
                clean_build=0
                echo -e "\n\n\e[33m   ---> Clean build cancelled. Proceeding with regular build.\e[0m\n\n"
            fi
            ;;
    esac
done

##################################################################
##################################################################
# Declare build options
##################################################################
##################################################################

# env control arguments
BUILD_DIR="${build_dir:-build}"
RUN_TESTS="${run_tests:-0}"
BUILD_TESTS="${build_tests:-1}"
NUMBER_OF_PROCESSORS=${number_of_processors:-1}
CLEAN_BUILD=${clean_build:-0}
BUILD_TYPE=${build_type:-Release}

# build directory
mkdir -p "${BUILD_DIR}"

cd "${BUILD_DIR}" || {
  echo 'entry into build folder failed'
  exit 1
}

##################################################################
##################################################################
# Configure CMake
##################################################################
##################################################################

# configuration step
cmake -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
  -DCMAKE_CXX_STANDARD=14 \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DTUDAT_BUILD_TESTS="${BUILD_TESTS}" \
  ..

##################################################################
##################################################################
# Build
##################################################################
##################################################################

# if required by user, clean
if [ "${CLEAN_BUILD}" = "1" ]; then
    cmake --build . --target clean -j"${NUMBER_OF_PROCESSORS}"
    printf "\nClean finished\n\n"
fi

# build
cmake --build . -j"${NUMBER_OF_PROCESSORS}"
