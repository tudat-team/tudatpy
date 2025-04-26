import argparse
import ast
from pathlib import Path
import os
import sys
import subprocess
import tempfile
import shutil
from typing import TypeGuard


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

    # Default indentation length in pybind11-stubgen
    indentation: str = " " * 4

    # Ignored modules and methods
    ignored_modules: list[str] = ["temp", "io"]
    ignored_methods: list[str] = ["_pybind11_conduit_v1_"]

    def __init__(self) -> None:

        # Initialize attributes from user input
        self.args = StubParser().parse_args()
        self.root_suffix = self.args.root_suffix

        # Resolve build directory
        self.build_dir = Path(self.args.build_dir).absolute()
        if not self.build_dir.exists():
            raise FileNotFoundError(
                f"Build directory {self.build_dir} does not exist."
            )

        # Resolve source directory of tudatpy
        self.tudatpy_source_dir = Path(__file__).parent / "tudatpy/src/tudatpy"
        if not self.tudatpy_source_dir.exists():
            raise FileNotFoundError(
                f"Source directory {self.tudatpy_source_dir} does not exist."
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

        # Define the path to the stub directory
        self.stubs_dir = self.build_dir / "tudatpy-stubs"

        return None

    def __create_mock_environment(
        self, tmp: Path
    ) -> tuple[dict[str, str], Path]:
        """Mock installation of tudatpy in a temporary directory

        pybind11-stubgen can only generate stubs for packages that are installed in the Python environment. To allow for stub generation without having to install the library, this function creates a mock python environment by installing tudatpy in a temporary directory, copying the current environment (os.environ) and adding the temporary directory to the PYTHONPATH. This fools pybind11-stubgen into thinking that the library is installed and makes it possible to integrate stub generation into the build process.
        """

        # Create mock installation directory in tmp
        mock_install_directory = tmp / "mock_install"
        mock_install_directory.mkdir(exist_ok=False, parents=False)

        # Install tudatpy in the mock installation directory
        shutil.copytree(
            self.build_dir.parent / "tudatpy/src/tudatpy",
            mock_install_directory / "tudatpy",
        )
        shutil.copy(
            self.build_dir / "tudatpy/src/tudatpy/kernel.so",
            mock_install_directory / "tudatpy/kernel.so",
        )

        # Create copy of the current environment with mock installation directory in PYTHONPATH
        env = os.environ.copy()
        if "PYTHONPATH" not in env:
            env["PYTHONPATH"] = str(mock_install_directory)
        else:
            env["PYTHONPATH"] += f":{mock_install_directory}"

        return env, mock_install_directory

    def __parse_script(self, stub_path: Path) -> ast.Module:

        # Check that stub exists and parse content
        if not stub_path.exists():
            raise FileNotFoundError(f"Stub {stub_path} does not exist.")

        # Read text and replace backslashes with a placeholder
        # This is necessary to handle special characters in docstrings
        text = stub_path.read_text().splitlines()
        normalized_text = "\n".join([line.replace("\\", "¿") for line in text])

        return ast.parse(normalized_text)

    def __unparse_script(self, content: ast.Module) -> str:

        # Unparse content and replace placeholder with backslashes
        unparsed_lines = ast.unparse(content).splitlines()
        unparsed_content = "\n".join(
            [line.replace("¿", "\\") for line in unparsed_lines]
        )

        return unparsed_content

    def fix_tudatpy_imports(
        self, module: ast.Module, stub: Path, tmp: Path
    ) -> ast.Module:
        """Fix imports from tudatpy

        - Replace kernel imports with submodule imports
        - Make tudatpy imports relative
        """

        # Find path of stub relative to base directory
        relative_path = stub.relative_to(self.stubs_dir)
        level = len(relative_path.parts)

        # Container for replacements of absolute imports in code
        import_replacements: dict[str, str] = {}

        # Process all statements in the module
        updated_body: list[ast.stmt] = []
        for statement in module.body:

            # Pybind11-stubgen only does absolute imports, not import from
            if isinstance(statement, ast.Import):

                # Multiple items can be imported in one statement
                # Loop over them and mark them as external or internal
                # Internal imports will be replaced with relative imports
                imported_items = {}
                for alias in statement.names:

                    if "tudatpy.kernel." not in alias.name:
                        imported_items[alias.name] = "external"
                    else:

                        # Remove tudatpy.kernel. prefix from import
                        relative_alias = alias.name.replace(
                            "tudatpy.kernel.", ""
                        )
                        imported_items[relative_alias] = "internal"

                # If multiple items are imported in the statement, turn each
                # into an individual import statement
                for name, import_type in imported_items.items():

                    if import_type == "external":
                        # Keep original import statement
                        updated_stmt = ast.Import([ast.alias(name)])
                    else:
                        # Separate module from imported submodule
                        _components = name.split(".")
                        if len(_components) == 1:
                            _module = ""
                            _name = name
                        else:
                            _module = ".".join(name.split(".")[:-1])
                            _name = name.split(".")[-1]

                        # Create import..from statement
                        updated_stmt = ast.ImportFrom(
                            module=_module,
                            names=[ast.alias(_name)],
                            level=level,
                        )

                        # Add submodule name to import replacements
                        import_replacements[f"tudatpy.kernel.{name}"] = _name

                    # Add import statement to updated body
                    updated_body.append(updated_stmt)
            else:
                updated_body.append(statement)

        # Update body of module
        module.body = updated_body

        # Replace calls to kernel imports with submodule imports
        for original, replacement in import_replacements.items():

            text = ast.unparse(module).replace(original, replacement)
            module = ast.parse(text)

        return module

    def fix_external_imports(self, module: ast.Module) -> ast.Module:
        """Fix import section of stub

        - Remove `from __future__ import annotations`, added automatically by pybind11-stubgen
        - Ensure that `typing` is imported [Don't remember why this is needed]
        """

        # Initialize containers and flags
        updated_body: list[ast.stmt] = []
        includes_typing: bool = False

        # Process all statements in the module
        for statement in module.body:

            # Import from statements
            if isinstance(statement, ast.ImportFrom):

                # Remove `from __future__ import annotations`
                if (
                    ast.unparse(statement)
                    == "from __future__ import annotations"
                ):
                    continue

                # # Replace kernel imports with submodule imports
                # if (
                #     statement.module is not None
                #     and "kernel." in statement.module
                # ):
                #     statement.module = statement.module.replace("kernel.", "")

            # Regular import statements
            if isinstance(statement, ast.Import):

                # Check if `import typing` is present
                if ast.unparse(statement) == "import typing":
                    includes_typing = True

                # # Replace kernel imports with submodule imports
                # for alias in statement.names:
                #     if "kernel." in alias.name:
                #         alias.name = alias.name.replace("kernel.", "")

            # Add statement to updated body
            updated_body.append(statement)

        # Update body of module
        module.body = updated_body

        # Import `typing` if not already imported
        if not includes_typing:
            module.body.insert(0, ast.Import([ast.alias("typing")]))

        return module

    def remove_autogenerated_methods(self, module: ast.Module) -> ast.Module:
        """Remove autogenerated methods from the stub

        Pybind11 automatically adds methods to wrapped classes. This function removes them from the stub.
        """

        # Initialize updated body of the module
        updated_body: list[ast.stmt] = []

        # Process all statements in the module
        for statement in module.body:

            # Autogenerated methods are added to classes
            if not isinstance(statement, ast.ClassDef):
                updated_body.append(statement)
                continue

            # Initialize updated class body
            updated_class_body: list[ast.stmt] = []
            for item in statement.body:

                # Skip if item is an autogenerated method
                if (
                    isinstance(item, ast.FunctionDef)
                    and item.name in self.ignored_methods
                ):
                    continue

                # Add item to updated class body
                updated_class_body.append(item)

            # If class body is empty, add an ellipsis
            if len(updated_class_body) == 0:
                updated_class_body.append(ast.Expr(ast.Constant(Ellipsis)))

            # Update class body
            statement.body = updated_class_body
            updated_body.append(statement)

        # Update module and return
        module.body = updated_body
        return module

    def adjust_docstring_indentation(self, module: ast.Module) -> ast.Module:
        """Adjust indentation of docstring in a module

        From https://docs.python.org/3/library/ast.html#ast.get_docstring, the types that can have a docstring are: ast.Module, ast.ClassDef, ast.FunctionDef, and ast.AsyncFunctionDef.
        """

        # Define container for updated body
        updated_body: list[ast.stmt] = []

        # Process all statements in the module
        for statement in module.body:

            # Skip if statement cannot have a docstring
            # if not self.__can_have_docstring(statement):
            if not isinstance(
                statement,
                (
                    ast.ClassDef,
                    ast.FunctionDef,
                    ast.AsyncFunctionDef,
                    ast.Module,
                ),
            ):
                updated_body.append(statement)
                continue

            # Get docstring and skip if not present
            docstring = ast.get_docstring(statement)
            if docstring is None:
                updated_body.append(statement)
                continue

            # Get size of docstring indentation
            if isinstance(statement, ast.Module):
                indentation_level = 0
            else:
                indentation_level = statement.col_offset
            docstring_indentation = (indentation_level + 1) * self.indentation

            # Adjust indentation
            indented_lines: list[str] = [
                docstring_indentation + line for line in docstring.split("\n")
            ]
            indented_lines[0] = indented_lines[0].lstrip()
            docstring = "\n".join(indented_lines)

            # Update docstring
            statement.body[0] = ast.Expr(ast.Constant(docstring))
            updated_body.append(statement)

        # Update body of module
        module.body = updated_body
        return module

    def __create_stubs_directory_structure(self) -> None:

        for item in self.tudatpy_source_dir.rglob("*"):

            # Skip if not a directory or if it is a cache directory
            if not item.is_dir() or item.name == "__pycache__":
                continue

            # Make path relative to source directory
            item = item.relative_to(self.tudatpy_source_dir)

            # Create directory if the module should not be ignored
            if item.parts[0] not in self.ignored_modules:
                (self.stubs_dir / item).mkdir(exist_ok=True, parents=True)

        return None

    def __incorrect_init(self, stub: Path) -> bool:
        """Check if an __init__.pyi stub contains definitions"""

        # Check if the stub is an __init__.pyi file
        if stub.name != "__init__.pyi":
            raise ValueError(
                "Attempted to check for definitions in a non-__init__ stub: "
                f"{stub}"
            )

        # Parse content
        module = self.__parse_script(stub)
        for statement in module.body:

            # Should only contain imports and __all__ or __version__ assignments
            if not isinstance(
                statement, (ast.Import, ast.ImportFrom, ast.Assign)
            ):
                return True

            # If statement is not an assignment, it is an import -> No action
            if not isinstance(statement, ast.Assign):
                continue

            # Valid assignments have only one target
            if len(statement.targets) > 1:
                return True
            _target = statement.targets[0]

            # In valid assignments, the target must be of ast.Name type
            if not isinstance(_target, ast.Name):
                return True

            # Check if name is __all__ or __version__
            if _target.id not in ["__all__", "__version__"]:
                return True

        return False

    def __generate_default_kernel_stubs(
        self, tmp: Path, env: dict[str, str]
    ) -> None:

        # Generate stubs for tudatpy.kernel
        outcome = subprocess.run(
            [
                "pybind11-stubgen",
                "tudatpy.kernel",
                "-o",
                str(tmp),
                f"--root-suffix={self.root_suffix}",
                "--numpy-array-remove-parameters",
            ],
            env=env,
        )
        if outcome.returncode:
            raise RuntimeError("Failed to generate stub for tudatpy.kernel")

        # Relocate stubs in the build directory
        tmp_stubs_dir: Path = tmp / "tudatpy/kernel"
        for stub in tmp_stubs_dir.rglob("*.pyi"):

            # Get path relative to output directory of pybind11-stubgen
            relative_path = stub.relative_to(tmp_stubs_dir)

            # Handle special case for __init__.pyi
            if stub.name == "__init__.pyi":

                # If the stub includes definitions, it should not be considered
                # as an __init__.pyi file. Rename as extension.pyi
                if self.__incorrect_init(stub):
                    shutil.copy(
                        stub,
                        self.stubs_dir / relative_path.parent / "extension.pyi",
                    )
                # Regular __init__.pyi files are generated manually at a later
                # state -> Skip
                else:
                    continue
            # If the stub is not __init__.pyi, copy to the final directory
            # with name extension.pyi
            else:
                shutil.copy(
                    stub,
                    self.stubs_dir
                    / relative_path.with_suffix("")
                    / "extension.pyi",
                )

        return None

    def __generate_default_python_stubs(
        self, tmp: Path, mock_install: Path
    ) -> None:

        # Create directory for python stubs in tmp
        python_stubs = tmp / "python_stubs"
        python_stubs.mkdir(exist_ok=False, parents=False)

        # Generate stubs for Python source code
        for item in (mock_install / "tudatpy").rglob("*.py"):

            # Skip if __init__.py or if it is in __pycache__
            if "__pycache__" in str(item):
                continue

            # Skip if it belongs to an ignored module
            if (
                item.relative_to(mock_install / "tudatpy").parts[0]
                in self.ignored_modules
            ):
                continue

            # Generate stub with stubgen
            outcome = subprocess.run(
                [
                    "stubgen",
                    item,
                    "-o",
                    str(python_stubs),
                    "--include-docstrings",
                    "--parse-only",
                ]
            )
            if outcome.returncode:
                raise RuntimeError(f"Failed to generate stub for {item}")

        # Move stubs to build directory
        for stub in (python_stubs / "tudatpy").rglob("*.pyi"):
            shutil.copy(
                stub,
                self.stubs_dir / stub.relative_to(python_stubs / "tudatpy"),
            )

        return None

    def __fix_autogenerated_stubs(self, tmp: Path) -> None:

        for stub in self.stubs_dir.rglob("*.pyi"):

            module = self.__parse_script(stub)
            module = self.fix_external_imports(module)
            module = self.fix_tudatpy_imports(module, stub, tmp)
            module = self.remove_autogenerated_methods(module)
            module = self.adjust_docstring_indentation(module)
            content = self.__unparse_script(module)
            stub.write_text(content)

        return None

    def generate_stubs(self) -> None:

        # Create directory for stubs in build
        if self.stubs_dir.exists():
            shutil.rmtree(self.stubs_dir)
        self.stubs_dir.mkdir(exist_ok=False, parents=False)

        # Create directory directory structure for stubs
        self.__create_stubs_directory_structure()

        # Create temporary directory for stub generation
        with tempfile.TemporaryDirectory() as tmp_str:

            # Create path object for the tmp directory
            tmp = Path(tmp_str)

            # Generate mock environment with tudatpy installed
            env, mock_install = self.__create_mock_environment(tmp)

            # Generate default stubs for tudatpy.kernel
            self.__generate_default_kernel_stubs(tmp, env)

            # Generate Python stubs
            # self.__generate_default_python_stubs(tmp, mock_install)

            # Fix autogenerated stubs
            self.__fix_autogenerated_stubs(tmp)

        # Missing: Fix __init__.pyi stubs to expand star imports

        return None


if __name__ == "__main__":

    StubGenerator().generate_stubs()

    # StubGenerator().generate_stubs()
    # StubGenerator().relocate_stubs()

    # stubs.generate_kernel_stubs()
    # stubs.generate_kernel_stubs()
    # stubs.generate_from_build()
    # stubs.generate_kernel_stubs()

# class PostInstallCommand:

#     def run(self, install_lib):
#         self.py_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'mpool'))
#         self.build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build'))
#         self.target_dir = os.path.join(install_lib, 'mpool')
#         self._install_lib = install_lib
#         self.copy_python_file()
#         self.copy_shared_library()
#         self.generate_stub()

#     def generate_stub(self):
#         env = os.environ.copy()
#         if 'PYTHONPATH' not in env:
#             env['PYTHONPATH'] = os.path.abspath(self._install_lib)
#         else:
#             env['PYTHONPATH'] += ':' + os.path.abspath(self._install_lib)
#         subprocess.check_call(['pybind11-stubgen', 'mpool', '-o', self._install_lib, '--ignore-all-errors'], env=env)
