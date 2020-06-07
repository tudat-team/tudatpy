#!/usr/bin/env bash

# Core deps.
sudo apt-get install build-essential cmake

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
    -DTUDAT_BUILD_TUDAT_TUTORIALS=off \
    -DTUDAT_BUILD_WITH_SOFA_INTERFACE=on \
    -DTUDAT_BUILD_WITH_SPICE_INTERFACE=on \
    ..

make  # -j2: Kills the remote build server.

# Log CircleCI results in JUnit in ~/tmp/test-reports
if [ -z "$CIRCLECI"  ]; then
  mkdir -p "$HOME/tmp/test-reports"
  touch "$HOME/tmp/test-reports/test.log"
  for filepath in ./tests/test_*; do
    filename="$(basename -- $filepath)"
    $filepath --log_format=JUNIT --log_level=all --log_sink="$HOME/tmp/test-reports/${filename}.xml" >> "$HOME/tmp/test-reports/test.log"
  done
fi

#ctest -T Test
#
## Extract test.xml results to $HOME/tmp/test-reports
#mkdir -p $HOME/tmp/test-reports
#find ./Testing -name \*.xml -exec cp {} $HOME/tmp/test-reports/test.xml \;

make install
