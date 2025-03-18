import sys
from pathlib import Path
import os

# New
######################################

import argparse
from build import BuildParser, setup_build_dir

# Globals
CONDA_PREFIX = Path(os.environ["CONDA_PREFIX"]).resolve()


# Argument parser
class InstallParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__(
            prog="install.py",
            description="Install Tudat and TudatPy in active conda environment",
        )

        self.add_argument(
            "--build-dir",
            metavar="<path>",
            type=str,
            default="build",
            help="Build directory",
        )

        # Choose what to install
        self.add_argument(
            "--no-tudat",
            dest="install_tudat",
            action="store_false",
            help="Do not install Tudat",
        )
        self.add_argument(
            "--no-tudatpy",
            dest="install_tudatpy",
            action="store_false",
            help="Do not install TudatPy",
        )
        self.add_argument(
            "--no-stubs",
            dest="install_stubs",
            action="store_false",
            help="Do not install TudatPy stubs",
        )

        return None


# Editable installation
def install(src: Path, dst: Path, manifest, iterate: bool = False) -> None:
    """Create symbolic link from source to destination and log in manifest.

    :param src: Source directory or file
    :param dst: Destination directory or file
    :param manifest: Writable installation manifest
    :param iterate: Whether to iterate over source directory
    """

    if iterate:
        for src_file in src.iterdir():
            dst_file = dst / src_file.name

            if not dst_file.is_symlink():
                dst_file.symlink_to(src_file)
            else:
                print(f"Already installed: {src_file}")
            manifest.write(str(dst_file) + "\n")
    else:
        if not dst.is_symlink():
            dst.symlink_to(src, target_is_directory=src.is_dir())
        else:
            print(f"Already installed: {src}")
        manifest.write(str(dst) + "\n")

    return None


if __name__ == "__main__":

    # Parse arguments
    args = InstallParser().parse_args()

    # Check if TudatPy is already installed
    tudatpy_installed = True
    try:
        from tudatpy import __version__
    except ImportError:
        tudatpy_installed = False

    if tudatpy_installed:
        print("TudatPy is already installed!")
        exit(0)

    # Source and destination directories
    build_dir = Path(args.build_dir).resolve()
    if not build_dir.exists():
        setup_build_dir(BuildParser().parse_args(None), build_dir)

    tudat_dir = (build_dir.parent / "tudat").resolve()
    tudatpy_dir = (build_dir.parent / "tudatpy").resolve()
    pylib_prefix = (
        Path(sys.exec_prefix)
        / sys.platlibdir
        / f"python{sys.version_info.major}.{sys.version_info.minor}"
        / "site-packages"
    ).resolve()

    # Perform installation
    with open(build_dir / "custom-manifest.txt", "w") as manifest:

        # Tudat
        if args.install_tudat:
            # Tudat static libraries
            install(
                build_dir / "lib", CONDA_PREFIX / "lib", manifest, iterate=True
            )

            # Tudat headers
            # (conda_prefix / "include/tudat").mkdir(parents=True, exist_ok=True)
            # No need to create include/tudat. It is created by tudat-resources
            install(
                build_dir / "tudat/include/tudat/config.hpp",
                CONDA_PREFIX / "include/tudat/config.hpp",
                manifest,
            )
            for item in (tudat_dir / "include/tudat").iterdir():
                install(
                    item, CONDA_PREFIX / "include/tudat" / item.name, manifest
                )

            # Tudat cmake files
            (CONDA_PREFIX / "lib/cmake/tudat").mkdir(
                parents=True, exist_ok=True
            )
            for item in (build_dir / "tudat").iterdir():
                if item.suffix == ".cmake" and "tudat" in item.name.lower():
                    install(
                        item,
                        CONDA_PREFIX / "lib/cmake/tudat" / item.name,
                        manifest,
                    )
            manifest.write(str(CONDA_PREFIX / "lib/cmake/tudat") + "\n")

        # Tudatpy
        if args.install_tudatpy:
            install(
                tudatpy_dir / "src/tudatpy",
                pylib_prefix / "tudatpy",
                manifest,
            )

        # Tudatpy stubs
        if args.install_stubs:
            if Path(tudatpy_dir / "src/tudatpy-stubs").exists():
                install(
                    tudatpy_dir / "src/tudatpy-stubs",
                    pylib_prefix / "tudatpy-stubs",
                    manifest,
                )

#####################################


# def usage() -> None:
#     """Print usage information."""
#     print("Usage: python install.py [OPTIONS]", end="\n\n")
#     print("Options:")
#     print("  --build-dir <path>  Select build directory")
#     print("  --help, -h          Display this information")
#     return None


# if __name__ == "__main__":

#     ENVIRONMENT = os.environ
#     ARGUMENTS = {"BUILD_DIR": "build"}

#     # Define arguments
#     args = iter(sys.argv[1:])
#     for arg in args:
#         if arg == "--build-dir":
#             ARGUMENTS["BUILD_DIR"] = next(args)
#         elif arg in ("--help", "-h"):
#             usage()
#             exit(0)
#         else:
#             usage()
#             raise ValueError("Invalid argument")

#     # Source and destination directories
#     build_dir = Path(ARGUMENTS["BUILD_DIR"]).resolve()
#     tudat_dir = (build_dir.parent / "tudat").resolve()
#     tudatpy_dir = (build_dir.parent / "tudatpy").resolve()
#     conda_prefix = Path(ENVIRONMENT["CONDA_PREFIX"]).resolve()
#     pylib_prefix = (
#         Path(sys.exec_prefix)
#         / sys.platlibdir
#         / f"python{sys.version_info.major}.{sys.version_info.minor}"
#         / "site-packages"
#     ).resolve()

#     # Perform installation
#     with open(build_dir / "custom-manifest.txt", "w") as manifest:

#         # Tudat static libraries
#         install(build_dir / "lib", conda_prefix / "lib", manifest, iterate=True)

#         # Tudat headers
#         # (conda_prefix / "include/tudat").mkdir(parents=True, exist_ok=True)
#         # No need to create include/tudat. It is created by tudat-resources
#         install(
#             build_dir / "tudat/include/tudat/config.hpp",
#             conda_prefix / "include/tudat/config.hpp",
#             manifest,
#         )
#         for item in (tudat_dir / "include/tudat").iterdir():
#             install(item, conda_prefix / "include/tudat" / item.name, manifest)

#         # Tudat cmake files
#         (conda_prefix / "lib/cmake/tudat").mkdir(parents=True, exist_ok=True)
#         for item in (build_dir / "tudat").iterdir():
#             if item.suffix == ".cmake" and "tudat" in item.name.lower():
#                 install(
#                     item, conda_prefix / "lib/cmake/tudat" / item.name, manifest
#                 )
#         manifest.write(str(conda_prefix / "lib/cmake/tudat") + "\n")

#         # Tudatpy
#         (pylib_prefix / "tudatpy").mkdir(parents=True, exist_ok=True)

#         install(
#             build_dir / "tudatpy/tudatpy/_version.py",
#             pylib_prefix / "tudatpy/_version.py",
#             manifest,
#         )

#         install(
#             build_dir / "tudatpy/tudatpy/kernel.so",
#             pylib_prefix / "tudatpy/kernel.so",
#             manifest,
#         )

#         for item in (tudatpy_dir / "tudatpy").iterdir():
#             if item.name not in (
#                 "_version.py.in",
#                 "CMakeLists.txt",
#                 "__pycache__",
#             ):
#                 install(item, pylib_prefix / "tudatpy" / item.name, manifest)
