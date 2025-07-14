from pathlib import Path
import sys
import os
import argparse
import shutil


class RemoveParser(argparse.ArgumentParser):

    def __init__(self) -> None:
        super().__init__(
            prog="uninstall.py",
        )

        self.add_argument(
            "--build-dir",
            metavar="<path>",
            type=str,
            default="build",
            help="Build directory",
        )

        return None


class Remover:

    def __init__(self) -> None:

        self.args = RemoveParser().parse_args()

        # Resolve build directory
        self.build_dir = Path(self.args.build_dir).resolve()
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

        # Resolve installation manifest
        self.manifest_dir = self.build_dir / "manifests"
        self.manifest = self.manifest_dir / f"{self.conda_prefix.name}.txt"
        if not self.manifest.exists():

            # Backwards compatibility: Try to look for the old manifest
            self.manifest = self.build_dir / "custom-manifest.txt"
            if not self.manifest.exists():
                raise FileNotFoundError(
                    "Installation manifest not found for active "
                    f"conda environment: {self.conda_prefix.name}"
                )

            # Check if the format of the old manifest is compatible
            content = self.manifest.read_text().splitlines()
            for line in content:
                if len(line.split(" ")) != 2 and line != "":
                    print(
                        "WARNING\n"
                        "It looks like your current installation of tudatpy was"
                        "\ninstalled using a version of `tudat-bundle` that is"
                        "\nincompatible with this script.\n"
                        "You can find instructions on how to remove it here:\n"
                        "https://github.com/tudat-team/tudat-bundle/wiki/backwards-compatibility"
                    )
                    exit(1)

        return None

    def preserve_remaining(self) -> list[str]:

        preserve = []
        for line in self.manifest.read_text().splitlines():

            code, _element = line.split(" ")
            element = Path(_element.strip())

            if element.exists():
                preserve.append(f"{code} {element}")

        # Write remaining entries to manifest
        with self.manifest.open("w") as manifest:
            manifest.write("\n".join(preserve))

        return preserve

    def remove(self) -> None:

        links = []
        directories = []
        remaining = {}
        for line in self.manifest.read_text().splitlines():

            code, _element = line.split(" ")
            element = Path(_element.strip())

            # Parse manifest for invalid entries before starting removal
            match code:
                case "000":
                    if not element.is_symlink():
                        raise ValueError(
                            "Aborting uninstall: "
                            f"Not a symlink with code 000: {element}"
                        )
                    links.append(element)
                case "999":
                    if not element.is_dir():
                        raise ValueError(
                            "Aborting uninstall: "
                            f"Not a directory with code 999: {element}"
                        )
                    if element.name not in ["tudat", "tudatpy"]:
                        raise ValueError(
                            "Aborting uninstall: Attempted to delete unexpected"
                            f" directory: {element}\n"
                            "This error is a safeguard and it should not occur."
                            "Please, contact the tudatpy team for support. "
                        )
                    directories.append(element)
                case _:
                    raise ValueError(
                        "Aborting uninstall: "
                        f"Unknown code {code} for {element}"
                    )

        # Remove links
        for link in links:
            try:
                link.unlink()
            finally:
                self.preserve_remaining()

        # Remove directories
        for directory in directories:
            try:
                shutil.rmtree(directory)
            finally:
                self.preserve_remaining()

        # Remove manifest
        remaining = self.preserve_remaining()
        if len(remaining) == 0:
            self.manifest.unlink()
        else:
            for item in remaining:
                print(f"Failed to remove: {item.split(' ')[1]}")


if __name__ == "__main__":

    Remover().remove()
