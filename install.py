import sys
from pathlib import Path
import os


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


def usage() -> None:
    """Print usage information."""
    print("Usage: python install.py [OPTIONS]", end="\n\n")
    print("Options:")
    print("  --build-dir <path>  Select build directory")
    print("  --help, -h          Display this information")
    return None


if __name__ == "__main__":

    ENVIRONMENT = os.environ
    ARGUMENTS = {"BUILD_DIR": "build"}

    # Define arguments
    args = iter(sys.argv[1:])
    for arg in args:
        if arg == "--build-dir":
            ARGUMENTS["BUILD_DIR"] = next(args)
        elif arg in ("--help", "-h"):
            usage()
            exit(0)
        else:
            usage()
            raise ValueError("Invalid argument")

    # Source and destination directories
    build_dir = Path(ARGUMENTS["BUILD_DIR"]).resolve()
    tudat_dir = (build_dir.parent / "tudat").resolve()
    tudatpy_dir = (build_dir.parent / "tudatpy").resolve()
    conda_prefix = Path(ENVIRONMENT["CONDA_PREFIX"]).resolve()
    pylib_prefix = (
        Path(sys.exec_prefix)
        / sys.platlibdir
        / f"python{sys.version_info.major}.{sys.version_info.minor}"
        / "site-packages"
    ).resolve()

    # Perform installation
    with open(build_dir / "custom-manifest.txt", "w") as manifest:

        # Tudat static libraries
        install(build_dir / "lib", conda_prefix / "lib", manifest, iterate=True)

        # Tudat headers
        # (conda_prefix / "include/tudat").mkdir(parents=True, exist_ok=True)
        # No need to create include/tudat. It is created by tudat-resources
        install(
            build_dir / "tudat/include/tudat/config.hpp",
            conda_prefix / "include/tudat/config.hpp",
            manifest,
        )
        for item in (tudat_dir / "include/tudat").iterdir():
            install(item, conda_prefix / "include/tudat" / item.name, manifest)

        # Tudat cmake files
        (conda_prefix / "lib/cmake/tudat").mkdir(parents=True, exist_ok=True)
        for item in (build_dir / "tudat").iterdir():
            if item.suffix == ".cmake" and "tudat" in item.name.lower():
                install(item, conda_prefix / "lib/cmake/tudat" / item.name, manifest)
        manifest.write(str(conda_prefix / "lib/cmake/tudat") + "\n")

        # Tudatpy
        (pylib_prefix / "tudatpy").mkdir(parents=True, exist_ok=True)

        install(
            build_dir / "tudatpy/tudatpy/_version.py",
            pylib_prefix / "tudatpy/_version.py",
            manifest,
        )

        for item in (tudatpy_dir / "src/tudatpy").iterdir():
            if item.name not in ("_version.py.in", "CMakeLists.txt", "__pycache__"):
                install(item, pylib_prefix / "tudatpy" / item.name, manifest)
