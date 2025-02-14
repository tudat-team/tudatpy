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

# construct conda environment path
export PYTHONPATH=".:$CONDA_PREFIX:$PYTHONPATH"

# move to tudatpy module inside build directory
cd $BUILD_DIR/tudatpy

# generate stubs
PYTHONPATH=".:$PYTHONPATH" stubgen -p tudatpy.kernel -o .

# move stubs out from tudatpy/kernel and into tudatpy/
rsync -avh --remove-source-files --ignore-times --exclude '/__init__.pyi' tudatpy/kernel/* .

# remove kernel stub directory
rm -r tudatpy/kernel/*