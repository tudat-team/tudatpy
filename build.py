import argparse
from pathlib import Path
import shutil
from contextlib import chdir
import subprocess
import os
import sys

TUDATPY_ROOT = Path("tudatpy/src/tudatpy").resolve()
CONDA_PREFIX = os.environ["CONDA_PREFIX"]
PYLIB_PREFIX = (
    Path(sys.exec_prefix)
    / sys.platlibdir
    / f"python{sys.version_info.major}.{sys.version_info.minor}"
    / "site-packages"
).resolve()


class BuildParser(argparse.ArgumentParser):
    """Argument parser for build script"""

    def __init__(self) -> None:

        # Default initialization
        super().__init__(prog="build.py", description="Build tools for Tudat")

        # Control basic behavior of the build script
        basic_group = self.add_argument_group("Basic configuration")
        basic_group.add_argument(
            "--tests",
            dest="build_tests",
            action="store_true",
            help="Build with tests [Default: False]",
        )
        basic_group.add_argument(
            "--clean",
            dest="clean_build",
            action="store_true",
            help="Remove pre-existing build directory [Default: False]",
        )
        basic_group.add_argument(
            "--skip-setup",
            dest="force_setup",
            action="store_false",
            help="Skip the execution of the CMake setup command. [Default: True]",
        )
        basic_group.add_argument(
            "-d", "--docs",
            dest="build_api_docs",
            action="store_true",
            help="Build API documentation. Output will be stored in <build-dir>/api-docs. [Default: False]",
        )

        # Control CMake behavior
        cmake_group = self.add_argument_group("CMake configuration")
        cmake_group.add_argument(
            "-j",
            metavar="<cores>",
            type=int,
            default=1,
            help="Number of processors to use",
        )
        cmake_group.add_argument(
            "--cxx-standard",
            metavar="<std>",
            default="14",
            help="C++ standard to use",
        )
        cmake_group.add_argument(
            "--build-dir",
            metavar="<dir>",
            default="build",
            help="Build directory",
        )
        cmake_group.add_argument(
            "--build-type",
            metavar="<type>",
            default="Release",
            help="Build type (e.g., Release, Debug)",
        )

        # Misc
        misc_group = self.add_argument_group("Miscellaneous")
        misc_group.add_argument(
            "--output-to-file",
            dest="output_to_file",
            action="store_true",
            help="Output logs to file instead of terminal [Default: False]",
        )
        misc_group.add_argument(
            "--verbose",
            action="store_true",
            help="Verbose output during build",
        )

        return None


class Builder:

    def __init__(self) -> None:

        # Command line arguments to control the build
        self.args = BuildParser().parse_args()

        # Build configuration attributes
        self.build_dir = Path(self.args.build_dir).resolve()
        # self.skip_tudat = "OFF" if self.args.build_tudat else "ON"
        # self.skip_tudatpy = "OFF" if self.args.build_tudatpy else "ON"
        self.build_tests = "ON" if self.args.build_tests else "OFF"

        return None

    def build_libraries(self) -> None:

        # If clean build is requested, delete build directory
        if self.build_dir.exists() and self.args.clean_build:

            # TODO: Logger
            print(f"Removing pre-existing build directory: {self.build_dir}")
            shutil.rmtree(self.build_dir)

        # If output to file is requested, set up output redirection
        if self.args.output_to_file:
            _output_dest = self.build_dir / "build_output.txt"
        else:
            _output_dest = None

        # If build directory still exists, skip setup
        if not self.build_dir.exists() or self.args.force_setup:

            # Set up build directory
            self.build_dir.mkdir(parents=True, exist_ok=True)
            with chdir(self.build_dir):

                if _output_dest is None:
                    outcome = subprocess.run(
                        [
                            "cmake",
                            # f"-DSKIP_TUDAT={self.skip_tudat}",
                            # f"-DSKIP_TUDATPY={self.skip_tudatpy}",
                            f"-DCMAKE_PREFIX_PATH={CONDA_PREFIX}",
                            f"-DCMAKE_INSTALL_PREFIX={CONDA_PREFIX}",
                            f"-DCMAKE_CXX_STANDARD={self.args.cxx_standard}",
                            "-DBoost_NO_BOOST_CMAKE=ON",
                            f"-DCMAKE_BUILD_TYPE={self.args.build_type}",
                            f"-DTUDAT_BUILD_TESTS={self.build_tests}",
                            "-B",
                            f"{self.build_dir}",
                            "-S",
                            "..",
                        ]
                    )
                else:
                    with _output_dest.open("w") as output_dest:
                        outcome = subprocess.run(
                            [
                                "cmake",
                                # f"-DSKIP_TUDAT={self.skip_tudat}",
                                # f"-DSKIP_TUDATPY={self.skip_tudatpy}",
                                f"-DCMAKE_PREFIX_PATH={CONDA_PREFIX}",
                                f"-DCMAKE_INSTALL_PREFIX={CONDA_PREFIX}",
                                f"-DCMAKE_CXX_STANDARD={self.args.cxx_standard}",
                                "-DBoost_NO_BOOST_CMAKE=ON",
                                f"-DCMAKE_BUILD_TYPE={self.args.build_type}",
                                f"-DTUDAT_BUILD_TESTS={self.build_tests}",
                                "-B",
                                f"{self.build_dir}",
                                "-S",
                                "..",
                            ],
                            stdout=output_dest,
                            stderr=output_dest,
                        )
            if outcome.returncode:
                shutil.rmtree(self.build_dir)
                raise RuntimeError("CMake setup failed")

        # Run build command
        with chdir(self.build_dir):

            build_command = ["cmake", "--build", ".", f"-j{self.args.j}"]
            if self.args.verbose:
                build_command.append("--verbose")
            outcome = subprocess.run(build_command)
            if outcome.returncode:
                exit(outcome.returncode)

        if self.args.build_api_docs:

            print("Building API documentation...")

            api_docs_build_dir = Path(self.build_dir) / "api-docs"
            api_docs_build_dir.mkdir(parents=True, exist_ok=True)

            # -E is added for a clean build, as changes in the compiled kernel are not detected otherwise
            api_docs_build_command = ["sphinx-build", "-b", "html", "tudatpy/docs/source", str(api_docs_build_dir), "-E"]

            if self.args.verbose:
                api_docs_build_command.append("--verbose")

            outcome = subprocess.run(api_docs_build_command)

            if outcome.returncode:
                exit(outcome.returncode)

        return None


if __name__ == "__main__":

    Builder().build_libraries()
