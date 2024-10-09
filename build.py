import os
import sys
import subprocess
from pathlib import Path
from contextlib import chdir

def usage():
    print("Usage: build.py [options]")
    print("Options:")
    print("  -h, --help            Show this help message and exit")
    print("  -j N                  Number of processors to use")
    print("  -c, --clean           Clean build")
    print("  --build-dir DIR       Build directory")
    print("  --no-tests            Do not build tests")
    print("  --cxx-std STD         C++ standard to use")
    print("  --build-type TYPE     Build type (e.g., Release, Debug)")
    print("  --output-to-file      Output logs to file instead of terminal")

if __name__ == "__main__":

    ARGUMENTS = {
        "NUMBER_OF_PROCESSORS": 1,
        "CLEAN_BUILD": False,
        "BUILD_DIR": "build",
        "BUILD_TESTS": True,
        "RUN_TESTS": True,
        "CXX_STANDARD": "14",
        "BUILD_TYPE": "Release",
        "OUTPUT_TO_FILE": False,
        "OUTPUT_FILE": "build/build_output.txt",
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
        elif arg == "--output-to-file":
            ARGUMENTS["OUTPUT_TO_FILE"] = True
        else:
            print("Invalid command")
            usage()
            exit(1)

    # Determine output destination
    if ARGUMENTS["OUTPUT_TO_FILE"]:
        output_dest = open(ARGUMENTS["OUTPUT_FILE"], "w")
    else:
        output_dest = None

    try:
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
                ],
                stdout=output_dest,
                stderr=output_dest
            )
            if outcome.returncode:
                exit(outcome.returncode)

            build_command = ["cmake", "--build", "."]
            if ARGUMENTS["CLEAN_BUILD"]:
                build_command.append("--target")
                build_command.append("clean")
            build_command.append(f"-j{ARGUMENTS['NUMBER_OF_PROCESSORS']}")
            outcome = subprocess.run(build_command, stdout=output_dest, stderr=output_dest)
            if outcome.returncode:
                exit(outcome.returncode)
    finally:
        if output_dest:
            output_dest.close()
