
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

# retrieve root directory (the parent directory of this file)
root=$(dirname "$(realpath "$0")")

export CMAKE_COLOR_DIAGNOSTICS=ON
cd "${root}/${BUILD_DIR}"
cmake \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
    -DCMAKE_CXX_STANDARD=14 \
    -DBoost_NO_BOOST_CMAKE=ON \
    -G 'CodeBlocks - Unix Makefiles' \
    -S "${root}" \
    -B "${root}/${BUILD_DIR}"