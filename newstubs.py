import argparse
import ast
from pathlib import Path
import os
import sys
import subprocess
import tempfile
import shutil


class StubParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__(prog="stubs.py")

        self.add_argument(
            "--root-suffix",
            dest="root_suffix",
            default="",
            help="To be appended to the name of the root directory of the stubs",
        )

        self.add_argument(
            "--build-dir",
            metavar="<path>",
            type=str,
            default="build",
            help="Build directory",
        )

        return None


class StubGenerator:

    def __init__(self) -> None:

        # Initialize attributes from user input
        self.args = StubParser().parse_args()
        self.root_suffix = self.args.root_suffix

        # Resolve build directory
        self.build_dir = Path(self.args.build_dir)
        if not self.build_dir.exists():
            raise FileNotFoundError(
                f"Build directory {self.build_dir} does not exist."
            )

        # Resolve conda prefix
        self.conda_prefix = Path(os.environ["CONDA_PREFIX"])
        if not self.conda_prefix.exists():
            raise FileNotFoundError(
                f"Conda prefix {self.conda_prefix} does not exist."
            )

        # Resolve pylib destination directory
        self.pylib_dir = (
            Path(sys.exec_prefix)
            / sys.platlibdir
            / f"python{sys.version_info.major}.{sys.version_info.minor}"
            / "site-packages"
        )
        if not self.pylib_dir.exists():
            raise FileNotFoundError(
                f"Python library directory {self.pylib_dir} does not exist."
            )

        return None

    def generate_kernel_stubs(self) -> None:

        # Create directory for stubs in build
        stubs_dir = self.build_dir / "tudatpy-stubs"
        if stubs_dir.exists():
            shutil.rmtree(stubs_dir)
        stubs_dir.mkdir(exist_ok=False, parents=False)

        with tempfile.TemporaryDirectory() as tmp_str:

            # Generate default stubs
            outcome = subprocess.run(
                [
                    "pybind11-stubgen",
                    "tudatpy.kernel",
                    "-o",
                    tmp_str,
                    f"--root-suffix={self.root_suffix}",
                    "--numpy-array-remove-parameters",
                ]
            )
            if outcome.returncode:
                raise RuntimeError("Failed to generate stub for tudatpy.kernel")

            # Show content of directory
            tmp = Path(tmp_str)
            for file in (tmp / "tudatpy/kernel").iterdir():
                dst = stubs_dir / file.name
                if dst.is_dir():
                    shutil.copytree(file, dst)
                else:
                    shutil.copyfile(file, dst)


if __name__ == "__main__":

    stubs = StubGenerator()
    stubs.generate_kernel_stubs()
