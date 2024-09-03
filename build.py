import sys
from pathlib import Path
from contextlib import chdir
import subprocess
import os


def usage() -> None:
    """Print usage information for the script."""

    print("Usage: python build.py [OPTIONS]", end="\n\n")
    print("Options:")
    print("  -h | --help                 Show this message")
    print("  -j                          Number of processors to use [Default: 1]")
    print("  -c | --clean                Clean after build")
    print("  --build-dir                 Build directory [Default: build]")
    print("  --no-tests                  Don't build tests")
    print("  --cxx-std                   C++ standard for compilation [Default: 14]")
    print(
        "  --build-type                Release, Debug, RelWithDebInfo [Default: Release]"
    )


if __name__ == "__main__":

    ARGUMENTS = {
        "NUMBER_OF_PROCESSORS": 1,
        "CLEAN_BUILD": False,
        "BUILD_DIR": "build",
        "BUILD_TESTS": True,
        "RUN_TESTS": True,
        "CXX_STANDARD": "14",
        "BUILD_TYPE": "Release",
    }
    CONDA_PREFIX = os.environ["CONDA_PREFIX"]

    # Parse input
    args = iter(sys.argv[1:])
    for arg in args:
        if arg in ("-h", "--help"):
            usage()
            exit(0)
        elif arg == "-j":
            ARGUMENTS["NUMBER_OF_PROCESSORS"] = next(args)
        elif arg in ("-c", "--clean"):
            ARGUMENTS["CLEAN_BUILD"] = True
        elif arg == "--build-dir":
            ARGUMENTS["BUILD_DIR"] = next(args)
        elif arg == "--no-tests":
            ARGUMENTS["BUILD_TESTS"] = False
        elif arg == "--cxx-std":
            ARGUMENTS["CXX_STANDARD"] = next(args)
        elif arg == "--build-type":
            ARGUMENTS["BUILD_TYPE"] = next(args)
        else:
            print("Invalid command")
            usage()
            exit(1)

    # Ensure build directory exists
    build_dir = Path(ARGUMENTS["BUILD_DIR"]).resolve()
    build_dir.mkdir(parents=True, exist_ok=True)

    # Build
    with chdir(build_dir):
        outcome = subprocess.run(
            [
                "cmake",
                f"-DCMAKE_PREFIX_PATH={CONDA_PREFIX}",
                f"-DCMAKE_INSTALL_PREFIX={CONDA_PREFIX}",
                f'-DCMAKE_CXX_STANDARD={ARGUMENTS["CXX_STANDARD"]}',
                "-DBoost_NO_BOOST_CMAKE=ON",
                f'-DCMAKE_BUILD_TYPE={ARGUMENTS["BUILD_TYPE"]}',
                f'-DTUDAT_BUILD_TESTS={ARGUMENTS["BUILD_TESTS"]}',
                "..",
            ]
        )
        if outcome.returncode:
            exit(outcome.returncode)

        build_command = ["cmake", "--build", "."]
        if ARGUMENTS["CLEAN_BUILD"]:
            build_command.append("--target")
            build_command.append("clean")
        build_command.append(f"-j{ARGUMENTS['NUMBER_OF_PROCESSORS']}")
        outcome = subprocess.run(build_command)
        if outcome.returncode:
            exit(outcome.returncode)
