from .build import ChangeDir, StubGenerator, TUDATPY_ROOT, BuildParser
import subprocess
from pathlib import Path

if __name__ == "__main__":

    # Parse arguments
    args = BuildParser().parse_args()

    # Build directory
    build_dir = Path(args.build_dir).resolve()

    # Generate stubs
    with ChangeDir("src"):

        try:
            from tudatpy import __version__  # type: ignore
        except ImportError:
            print("Skipping stub generation: tudatpy not installed")

        stub_generator = StubGenerator(clean=args.stubs_clean)
        stub_generator.generate_stubs(TUDATPY_ROOT)

    # Install stubs
    print("Installing stubs...")
    install_command = ["cmake", "--install", f"{build_dir}"]
    outcome = subprocess.run(install_command)
    if outcome.returncode:
        exit(outcome.returncode)
