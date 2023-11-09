#!/usr/bin/env bash

# env control arguments
build_dir="cmake-build-release-wsl"
export build_dir

# build kernel
bash build.sh -j10          # regular build
#bash build.sh -j10 -c true  # clean build

# expose kernel modules
cd tudatpy

declare -a KERNEL_MODULES_TO_BE_EXPOSED=(
    "astro"
    "trajectory_design"
    "constants"
    "interface"
    "math"
    "numerical_simulation"
)

for KERNEL_MODULE in "${KERNEL_MODULES_TO_BE_EXPOSED[@]}"; do
    conda run -n tudat-bundle python expose_kernel_module.py "${KERNEL_MODULE}" --build-dir "${build_dir}" --init
done
