from contextlib import chdir
from pathlib import Path
import subprocess
import sys

if __name__ == "__main__":

    root = Path(__file__).parent
    tudatpy_root = root / "src"

    # Extensions
    targets = [
        "astro",
        "constants",
        "interface",
        "math",
        "trajectory_design",
        # "numerical_simulation",
    ]
    ignore = ["__pycache__", "tests"]

    with chdir(tudatpy_root):

        # Loop over extensions
        for item in (tudatpy_root).rglob("**/*.so"):
            relpath = item.relative_to(tudatpy_root)
            if relpath.parts[1] in targets:
                extension_name = relpath.stem.split(".")[0]
                submodule_path = ".".join(relpath.parts[:-1])
                import_path = f"{submodule_path}.{extension_name}"
                subprocess.run(["stubgen", "-m", import_path, "-o", "."])

        for item in (tudatpy_root).rglob("**/*.py"):
            relpath = item.relative_to(tudatpy_root)
            if relpath.parts[1] in targets:
                subprocess.run(["stubgen", str(relpath), "-o", "."])
