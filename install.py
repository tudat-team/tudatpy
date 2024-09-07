import sys
from pathlib import Path
import os
import argparse
from build import setup_build_dir, BuildParser

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
        from tudatpy import __version__  # type: ignore
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

        # Tudatpy
        install(
            tudatpy_dir / "src/tudatpy",
            pylib_prefix / "tudatpy",
            manifest,
        )

        # Tudatpy stubs
        if Path(tudatpy_dir / "src/tudatpy-stubs").exists():
            install(
                tudatpy_dir / "src/tudatpy-stubs",
                pylib_prefix / "tudatpy-stubs",
                manifest,
            )
