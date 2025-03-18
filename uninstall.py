import sys
from pathlib import Path
import os
import argparse


def usage() -> None:
    """Print usage information."""
    print("Usage: python uninstall.py [OPTIONS]", end="\n\n")
    print("Options:")
    print("  --build-dir <path>  Select build directory")
    print("  --help, -h          Display this information")
    return None


# Globals
CONDA_PREFIX = Path(os.environ["CONDA_PREFIX"]).resolve()

# Argument parser
parser = argparse.ArgumentParser(
    prog="uninstall.py",
    description="Uninstall Tudat and TudatPy from active conda environment",
)
parser.add_argument(
    "--build-dir",
    metavar="<path>",
    type=str,
    default="build",
    help="Build directory",
)

if __name__ == "__main__":

    # Parse command line arguments
    args = parser.parse_args()

    # Ensure that build dir exists
    build_dir = Path(args.build_dir).resolve()
    if not build_dir.exists():
        raise FileNotFoundError("Failed to cd into build directory")

    # Look for installation manifest
    if not (build_dir / "custom-manifest.txt").exists():
        raise FileNotFoundError("Installation manifest not found")

    # Uninstall
    with (build_dir / "custom-manifest.txt").open() as manifest:
        for line in manifest:
            path = Path(line.strip())
            if path.exists():
                try:
                    path.unlink()
                except PermissionError:
                    path.rmdir()
            else:
                print(f"Not installed: {path}")
