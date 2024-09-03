import sys
from pathlib import Path
import os


def usage() -> None:
    """Print usage information."""
    print("Usage: python uninstall.py [OPTIONS]", end="\n\n")
    print("Options:")
    print("  --build-dir <path>  Select build directory")
    print("  --help, -h          Display this information")
    return None


if __name__ == "__main__":

    ARGUMENTS = {"BUILD_DIR": "build"}
    ENVIRONMENT = os.environ

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

    # Ensure that build dir exists
    build_dir = Path(ARGUMENTS["BUILD_DIR"]).resolve()
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
