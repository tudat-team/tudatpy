from pathlib import Path
import sys
import os
import argparse
import subprocess


class InstallParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__(
            prog="install.py",
            usage="Install tudat and tudatpy in the active conda environment.\n Note: The use without the -e flag, which performs an editable installation, is currently discouraged",
        )

        self.add_argument(
            "-e",
            dest="editable",
            action="store_true",
            help="Install in editable mode",
        )

        self.add_argument(
            "--build-dir",
            metavar="<path>",
            type=str,
            default="build",
            help="Build directory",
        )

        self.add_argument(
            "--skip-tudat",
            dest="install_tudat",
            action="store_false",
            help="Skip installation of tudat",
        )

        self.add_argument(
            "--skip-tudatpy",
            dest="install_tudatpy",
            action="store_false",
            help="Skip installation of tudatpy",
        )

        return None


class Installer:

    def __init__(self) -> None:

        self.args = InstallParser().parse_args()

        # Resolve build directory
        self.build_dir = Path(self.args.build_dir)
        if not self.build_dir.exists():
            raise FileNotFoundError(
                f"Build directory {self.build_dir} does not exist."
            )

        # Resolve base directories for tudat and tudatpy
        self.base_tudat = Path(__file__).parent / "tudat"
        if not self.base_tudat.exists():
            raise FileNotFoundError(
                f"Directory {self.base_tudat} does not exist."
            )
        self.base_tudatpy = self.base_tudat.parent / "tudatpy"
        if not self.base_tudatpy.exists():
            raise FileNotFoundError(
                f"Directory {self.base_tudatpy} does not exist."
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

        # Define path to installation manifest
        self.manifest_dir = self.build_dir / "manifests"
        self.manifest_dir.mkdir(parents=True, exist_ok=True)
        self.manifest = self.manifest_dir / f"{self.conda_prefix.name}.txt"

        # Backwards compatibility: Check if the old manifest exists
        old_manifest = self.build_dir / "custom-manifest.txt"
        if old_manifest.exists():
            print(
                "IMPORTANT WARNING:\n"
                "The installation manifest from an older version of this script"
                " was found.\nThis probably means that you installed tudatpy "
                "using an older version of\nthis script, and you did not "
                "uninstall it before updating tudat-bundle.\n\n"
                "Please, have a look at this document before proceeding:\n"
                "https://github.com/tudat-team/tudat-bundle/wiki/backwards-compatibility"
            )
            exit(0)

        if self.manifest.exists():
            raise RuntimeError("Delete current manifest before installation")
        # if self.manifest.exists():
        #     print(
        #         "WARNING: Installation manifest already exists."
        #         "Uninstalling libraries before performing new installation."
        #     )
        #     for line in self.manifest.read_text().splitlines():
        #         path = Path(line.strip())
        #         if path.exists():
        #             try:
        #                 path.unlink()
        #             except PermissionError:
        #                 path.rmdir()
        #         else:
        #             print(f"Not installed: {path}")
        #     self.manifest.unlink()

        self.directories = []
        self.manifest_list = []

        return None

    def link_single(self, src: Path, dst: Path) -> None:

        # Create symbolic link
        if dst.exists() or dst.is_symlink():
            print(f"Already installed: {dst}")
            return None

        # Ensure that source file exists
        src = src.absolute()
        if not src.exists():
            print(f"Skipping {src} because it does not exist")
            return None

        dst.symlink_to(src, target_is_directory=src.is_dir())

        # Update manifest list
        self.manifest_list.append(str(dst))

        return None

    def link_content(
        self,
        src_dir: Path,
        dst_dir: Path,
        extensions: str | list[str],
        _skip: list[str] | None = None,
        recursive: bool = False,
    ) -> None:

        # Get list of files to skip
        skip = _skip if _skip is not None else []

        # Turn extension into list if needed
        if isinstance(extensions, str):
            extensions = [extensions]
        assert isinstance(extensions, list)

        for ext in extensions:

            if recursive:
                pool = src_dir.rglob(f"*{ext}")
            else:
                pool = src_dir.glob(f"*{ext}")

            for src_file in pool:

                skip_flag = False
                for skip_item in skip:
                    if skip_item in str(src_file):
                        skip_flag = True
                        break
                if skip_flag:
                    continue

                # File to be created
                dst_file = dst_dir / src_file.relative_to(src_dir)

                # Create directory if it does not exist
                if not dst_file.parent.exists():

                    # Find the first directory that does not exist
                    base_dir = dst_file.parent
                    while not base_dir.parent.exists():
                        base_dir = base_dir.parent

                    # Figure out if it should be added to list of directories
                    is_subdir = False
                    for directory in self.directories:
                        if directory in base_dir.parents:
                            is_subdir = True
                            break

                    # If not a subdirectory of a created directory, add to list
                    if not is_subdir:
                        self.directories.append(base_dir)

                    # Create the original parent directory
                    dst_file.parent.mkdir(parents=True, exist_ok=False)

                # Create symbolic link
                self.link_single(src_file, dst_file)

        return None

    def editable_install(self) -> None:

        if self.args.install_tudat:

            # Install tudat static libraries
            self.link_content(
                self.build_dir / "lib",
                self.conda_prefix / "lib",
                [".a", ".lib"],
            )

            # Install tudat headers
            self.link_content(
                self.build_dir / "tudat/include/tudat",
                self.conda_prefix / "include/tudat",
                ".hpp",
            )
            self.link_content(
                self.base_tudat / "include/tudat",
                self.conda_prefix / "include/tudat",
                "",
            )

            # Install tudat CMake files
            self.link_content(
                self.build_dir / "tudat",
                self.conda_prefix / "lib/cmake/tudat",
                ".cmake",
                _skip=["cmake_install.cmake", "CTestTestfile.cmake"],
            )

        if self.args.install_tudatpy:

            # Install python files of tudatpy
            self.link_content(
                self.base_tudatpy / "src/tudatpy",
                self.pylib_dir / "tudatpy",
                [".py", ".typed"],
                recursive=True,
                _skip=["__pycache__"],
            )

            # Install kernel
            self.link_content(
                self.build_dir / "tudatpy/src/tudatpy",
                self.pylib_dir / "tudatpy",
                [".so", ".dll", ".dylib"],
            )

        # Filter elements inside created directories out of manifest
        manifest_list = []
        for item in self.manifest_list:
            item_path = Path(item)
            in_directory = False
            for directory in self.directories:
                if directory in item_path.parents:
                    in_directory = True
                    break
            if not in_directory:
                manifest_list.append("000 " + item)  # File label 000

        # Add created directories to manifest [With directory label 999 ]
        for directory in self.directories:
            manifest_list.append("999 " + str(directory))

        # Write manifest
        with self.manifest.open("w") as manifest:
            for item in manifest_list:
                manifest.write(f"{item}\n")

        return None

    def regular_install(self) -> None:

        print("WARNING")
        print("-----------------------------------------------------")
        print("Regular installation is performed via `cmake --install`")
        print("The current version of CMake does not include an `uninstall`")
        print("command, meaning that the user is responsible for uninstalling")
        print("the package manually.")
        print("Call this script with the `-e` flag to perform an editable")
        print("installation. In this case, you will be able to use the ")
        print("`uninstall.py` script to uninstall the package.")
        print("-----------------------------------------------------")
        input_request = (
            "Do you want to continue with the regular installation? [y/N] "
        )
        proceed = True if input(input_request).lower() == "y" else False
        if not proceed:
            print("Installation aborted.")
            return None

        outcome = subprocess.run(
            [
                "cmake",
                "--install",
                str(self.build_dir),
            ]
        )
        if outcome.returncode != 0:
            raise RuntimeError(
                f"Installation failed with error code {outcome.returncode}"
            )

        return None

    def install(self) -> None:

        if self.args.editable:
            self.editable_install()
        else:
            self.regular_install()

        return None


if __name__ == "__main__":

    Installer().install()
