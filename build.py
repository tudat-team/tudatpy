import argparse
from pathlib import Path
import shutil
from contextlib import chdir, contextmanager
import subprocess
import os
import ast
import tempfile
from dataclasses import dataclass
from typing import Generator


@dataclass
class Environment:
    """Mock environment

    :param variables: Dictionary of environment variables
    :param tmp: Path to the temporary directory
    :param prefix: Path to the mock installation prefix
    """

    variables: dict[str, str]
    tmp: Path
    prefix: Path


class BuildParser(argparse.ArgumentParser):
    """Argument parser for build script"""

    def __init__(self) -> None:

        # Default initialization
        super().__init__(prog="build.py", description="Build tools for Tudat")

        # Control basic behavior of the build script
        basic_group = self.add_argument_group("Basic configuration")
        basic_group.add_argument(
            "--tests",
            dest="build_tests",
            action="store_true",
            help="Build with tests [Default: False]",
        )
        basic_group.add_argument(
            "--clean",
            dest="clean_build",
            action="store_true",
            help="Remove pre-existing build directory [Default: False]",
        )
        basic_group.add_argument(
            "--skip-setup",
            dest="force_setup",
            action="store_false",
            help="Skip the execution of the CMake setup command. [Default: True]",
        )
        basic_group.add_argument(
            "--skip-build",
            dest="skip_build",
            action="store_true",
            help="Don't build the libraries [Default: False]",
        )
        basic_group.add_argument(
            "--skip-stubs",
            dest="skip_stubs",
            action="store_true",
            help="Skip the generation of stubs [Default: False]",
        )
        basic_group.add_argument(
            "--docs",
            dest="build_api_docs",
            action="store_true",
            help="Build API documentation. Output will be stored in <build-dir>/api-docs. [Default: False]",
        )
        basic_group.add_argument(
            "--stubs-only",
            dest="skip_build",
            action="store_true",
            help="Only update the stubs [Default: False]",
        )

        # Control CMake behavior
        cmake_group = self.add_argument_group("CMake configuration")
        cmake_group.add_argument(
            "-j",
            metavar="<cores>",
            type=int,
            default=1,
            help="Number of processors to use",
        )
        cmake_group.add_argument(
            "--cxx-standard",
            metavar="<std>",
            default="17",
            help="C++ standard to use",
        )
        cmake_group.add_argument(
            "--build-dir",
            metavar="<dir>",
            default="build",
            help="Build directory",
        )
        cmake_group.add_argument(
            "--build-type",
            metavar="<type>",
            default="Release",
            help="Build type (e.g., Release, Debug)",
        )
        cmake_group.add_argument(
            "--with-mcd",
            dest="build_with_mcd",
            action="store_true",
            help="Build with Mars Climate Database support [Default: False]",
        )

        # Misc
        misc_group = self.add_argument_group("Miscellaneous")
        misc_group.add_argument(
            "--output-to-file",
            dest="output_to_file",
            action="store_true",
            help="Output logs to file instead of terminal [Default: False]",
        )
        misc_group.add_argument(
            "--verbose",
            action="store_true",
            help="Verbose output during build",
        )

        # CI options
        ci_group = self.add_argument_group("CI options")
        ci_group.add_argument(
            "--build-github-actions",
            dest="build_github_actions",
            action="store_true",
            help="Build tudatpy with GitHub Actions",
        )

        return None


def is_star_import(statement: ast.stmt) -> bool:
    """Check if a statement is a star import

    :param statement: The statement to be checked
    :return: Whether the statement is a star import or not
    """

    # All star imports are of type: from X import *
    if not isinstance(statement, ast.ImportFrom):
        return False

    # The number of imported elements must be equal to 1 (the star)
    if len(statement.names) != 1:
        return False

    # The imported item must be *
    if statement.names[0].name != "*":
        return False

    return True


def module_has_init_stub(module_stubs_path: Path) -> bool:
    """Checks if an __init__.pyi stub exists for a module

    :param module_stubs_path: Path to the module directory in the stubs tree
    :return: Whether __init__.pyi exists in the directory of the module in the stubs tree
    """

    # Path to __init__.pyi
    init_stub_path = module_stubs_path / "__init__.pyi"
    return init_stub_path.exists()


class StubGenerator:

    # Default indentation length in pybind11-stubgen
    indentation: str = " " * 4

    # Ignored modules and methods
    ignored_modules: list[str] = ["temp", "io", "numerical_simulation"]
    ignored_methods: list[str] = ["_pybind11_conduit_v1_"]

    def __init__(self, build_dir: Path, mock_env: "Environment") -> None:

        # Ensure that build directory exists
        build_dir = Path(build_dir).absolute()
        if not build_dir.exists():
            raise FileNotFoundError(
                f"Failed to generate stubs: "
                f"Build directory {build_dir} does not exist."
            )

        # Source directory of tudatpy
        self.python_source_dir = Path(__file__).parent / "src/tudatpy"
        if not self.python_source_dir.exists():
            raise FileNotFoundError(
                f"Failed to generate stubs: "
                f"Source directory {self.python_source_dir} does not exist."
            )

        # Source directory for compiled extensions
        self.extension_source_dir = build_dir / "src/tudatpy"
        if not self.extension_source_dir.exists():
            raise FileNotFoundError(
                f"Failed to generate stubs: "
                f"Source directory {self.extension_source_dir} does not exist."
            )

        # Mock environment
        self.mock_env = mock_env

        # Define the path to the stub directory
        self.stubs_dir = build_dir / "tudatpy-stubs"

        return None

    def __parse_script(self, stub_path: Path) -> ast.Module:
        """Parse a python script

        Some of the docstrings in our code base contain special characters,
        which are not well handled by the ast module. This function reads the
        script as text, replaces backslashes with a placeholder (¿), and then
        parses the text. The complementary method __unparse_script reverses
        this process.

        :param stub_path: Path to the script
        :return: Parsed script as an ast.Module
        """

        # Check that stub exists and parse content
        if not stub_path.exists():
            raise FileNotFoundError(f"Stub {stub_path} does not exist.")

        # Read text and replace backslashes with a placeholder
        # This is necessary to handle special characters in docstrings
        text = stub_path.read_text().splitlines()
        normalized_text = "\n".join([line.replace("\\", "¿") for line in text])

        return ast.parse(normalized_text)

    def __unparse_script(self, content: ast.Module) -> str:
        """Unparse a python script

        Reverses the actions of __parse_script by replacing the placeholder
        (¿) with backslashes, and returning the content of the module as text.

        :param content: Parsed script as an ast.Module
        :return: Unparsed script as a string
        """

        # Unparse content and replace placeholder with backslashes
        unparsed_lines = ast.unparse(content).splitlines()
        unparsed_content = "\n".join(
            [line.replace("¿", "\\") for line in unparsed_lines]
        )

        return unparsed_content

    def __fix_tudatpy_imports(self, module: ast.Module, stub: Path) -> ast.Module:
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
                internal_imports: list[str] = []
                external_imports: list[ast.alias] = []
                for alias in statement.names:

                    if "tudatpy.kernel." not in alias.name:
                        external_imports.append(alias)
                    else:

                        # Remove tudatpy.kernel. prefix from import
                        relative_alias = alias.name.replace("tudatpy.kernel.", "")
                        internal_imports.append(relative_alias)

                # Keep original external import statements
                for alias in external_imports:
                    updated_body.append(statement)

                # Make internal imports relative
                for name in internal_imports:

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
                    updated_body.append(updated_stmt)

                    # Add submodule name to import replacements
                    import_replacements[f"tudatpy.kernel.{name}"] = _name

            else:
                updated_body.append(statement)

        # Update body of module
        module.body = updated_body

        # Replace calls to kernel imports with submodule imports
        for original, replacement in import_replacements.items():

            text = ast.unparse(module).replace(original, replacement)
            module = ast.parse(text)

        return module

    def __fix_external_imports(self, module: ast.Module) -> ast.Module:
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
                if ast.unparse(statement) == "from __future__ import annotations":
                    continue

            # Regular import statements
            if isinstance(statement, ast.Import):

                # Check if `import typing` is present
                if ast.unparse(statement) == "import typing":
                    includes_typing = True

            # Add statement to updated body
            updated_body.append(statement)

        # Update body of module
        module.body = updated_body

        # Import `typing` if not already imported
        if not includes_typing:
            module.body.insert(0, ast.Import([ast.alias("typing")]))

        return module

    def __remove_autogenerated_methods(self, module: ast.Module) -> ast.Module:
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

    def __adjust_docstring_indentation(self, module: ast.Module) -> ast.Module:
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

        for item in self.python_source_dir.rglob("*"):

            # Skip if not a directory or if it is a cache directory
            if not item.is_dir() or item.name == "__pycache__":
                continue

            # Make path relative to source directory
            item = item.relative_to(self.python_source_dir)

            # Create directory if the module should not be ignored
            if item.parts[0] not in self.ignored_modules:
                (self.stubs_dir / item).mkdir(exist_ok=True, parents=True)

        return None

    def __generate_default_kernel_stubs(self) -> None:

        print("Generating stubs for tudatpy.kernel...")

        # Generate stubs for tudatpy.kernel
        outcome = subprocess.run(
            [
                "pybind11-stubgen",
                "tudatpy.kernel",
                "-o",
                str(self.mock_env.tmp),
                "--numpy-array-wrap-with-annotated",
            ],
            env=self.mock_env.variables,
        )
        if outcome.returncode:
            raise RuntimeError("Failed to generate stub for tudatpy.kernel")

        # Relocate stubs in the build directory
        tmp_stubs_dir: Path = self.mock_env.tmp / "tudatpy/kernel"
        for stub in tmp_stubs_dir.rglob("*.pyi"):

            # Get path relative to output directory of pybind11-stubgen
            relative_path = stub.relative_to(tmp_stubs_dir)

            # Handle special case for __init__.pyi
            if stub.name == "__init__.pyi":

                shutil.copy(
                    stub,
                    self.stubs_dir / relative_path.parent / "extension.pyi",
                )

            # If the stub is not __init__.pyi, copy to the final directory
            # with name extension.pyi
            else:
                shutil.copy(
                    stub,
                    self.stubs_dir / relative_path.with_suffix("") / "extension.pyi",
                )

        return None

    def __generate_default_python_stubs(self) -> None:

        print("Generating stubs for Python source code...")

        # Create directory for python stubs in tmp
        python_stubs = self.mock_env.tmp / "python_stubs"
        python_stubs.mkdir(exist_ok=False, parents=False)

        # Generate stubs for Python source code
        for item in (self.mock_env.prefix / "tudatpy").rglob("*.py"):

            # Skip if __init__.py or if it is in __pycache__
            if item.name == "__init__.py" or "__pycache__" in str(item):
                continue

            # Skip if it belongs to an ignored module
            if (
                item.relative_to(self.mock_env.prefix / "tudatpy").parts[0]
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
                    "--quiet",
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

    def __fix_autogenerated_stubs(self) -> None:

        for stub in self.stubs_dir.rglob("*.pyi"):

            module = self.__parse_script(stub)
            module = self.__fix_external_imports(module)
            module = self.__fix_tudatpy_imports(module, stub)
            module = self.__remove_autogenerated_methods(module)
            module = self.__adjust_docstring_indentation(module)
            content = self.__unparse_script(module)
            stub.write_text(content)

        return None

    def __retrieve_items_in_all(
        self, statement: ast.Assign | ast.AnnAssign
    ) -> list[str]:
        """Retrieve items in __all__ = [item1, item2, ...] statement

        :param statement: __all__ statement
        :return: List of items in __all__
        """

        # If __all__ = list(), replace list() with []
        if isinstance(statement.value, ast.Call):
            if ast.unparse(statement.value.func) == "list":
                statement.value = ast.List(elts=[], ctx=ast.Load())

        if not isinstance(statement.value, ast.List):
            raise NotImplementedError(
                f"Failed to expand {ast.unparse(statement)}: "
                f"Expected __all__ = [...]"
            )

        # Update __all__ in the stub and generate expanded import statement
        all_contents: list[str] = []
        for _name in statement.value.elts:

            # All items in __all__ should be strings
            if not (isinstance(_name, ast.Constant) and isinstance(_name.value, str)):
                raise NotImplementedError(
                    f"Failed to expand {ast.unparse(statement)}: "
                    f"Unexpected item in __all__."
                )

            # Update __all__ list
            all_contents.append(_name.value)

        return all_contents

    def __find_all_statement(self, script: Path) -> ast.Assign | ast.AnnAssign:
        """Find __all__ statement in a script

        :param script: Path to the script
        :return: __all__ statement or None if not found
        """

        module = self.__parse_script(script)
        module_all: ast.Assign | ast.AnnAssign | None = None
        for statement in module.body:
            if isinstance(statement, ast.Assign) or isinstance(
                statement, ast.AnnAssign
            ):
                if "__all__" in ast.unparse(statement):
                    module_all = statement
                    break
        if module_all is None:
            raise ValueError(f"Script {script} does not contain an __all__ statement.")

        return module_all

    def __expand_kernel_star_import(
        self, statement: ast.ImportFrom, all_list: list[str]
    ) -> tuple[ast.ImportFrom | None, list[str]]:
        """Expand star import from kernel

        :param statement: Original star import statement
        :param all_list: Current list of items in __all__
        :return updated_statement: Expanded import statement or None if the kernel submodule is empty
        :return all_list: Updated list of items in __all__
        """

        # Ensure that the extension.pyi file exists for the module
        assert statement.module is not None  # Sanity
        extension_path = (
            self.stubs_dir
            / statement.module.replace("tudatpy.kernel.", "").replace(".", "/")
            / "extension.pyi"
        )
        if not extension_path.exists():
            raise FileNotFoundError(
                f"Failed to expand {ast.unparse(statement)}: "
                f"Stub {extension_path} does not exist."
            )

        # Find __all__ in extension.pyi
        extension_all = self.__find_all_statement(extension_path)

        # Create list of imported items and update __all__ list
        imported_items: list[ast.alias] = []
        for item in self.__retrieve_items_in_all(extension_all):
            all_list.append(item)
            imported_items.append(ast.alias(item))

        # If nothing is imported, skip the import statement
        if len(imported_items) == 0:
            return None, all_list

        # Create import statement
        import_all = ast.ImportFrom(
            module="extension",
            names=imported_items,
            level=1,
        )
        return import_all, all_list

    def __generate_single_init_stub(self, module_path: Path) -> None:
        """Generate stub for an __init__.py file

        The __init__.py files in tudatpy are subject to a series of rules:
        - If a module contains an extension, it must include a star import
        from `tudatpy.kernel.module_name`
        - Non-star import from `tudatpy.kernel` are not allowed
        - Relative kernel imports are not allowed
        - Imports from python scripts should be relative: `from .script import
        item, item...`
        - Submodules should only be imported if they are pure-python, otherwise
        the import statement is already included in the star import from
        `tudatpy.kernel`
        - The only allowed assign statements are __all__ and __version__. Other
        assign statements, and any other kind of statement, will be ignored
        during stub generation.

        This function generates the stub for the __init__.py file, under the
        assumption that these rules are followed. Not doing it might result in
        unexpected behavior during stub generation.

        :param module_path: Path to the module
        """

        # Initialize fresh stub file
        stub_path = module_path / "__init__.pyi"
        if stub_path.exists():
            stub_path.unlink()
        stub_path.touch()
        stub = self.__parse_script(stub_path)

        # Parse __init__.py of module
        init_path = (
            self.python_source_dir
            / module_path.relative_to(self.stubs_dir)
            / "__init__.py"
        )
        init = self.__parse_script(init_path)

        # Initialize containers for the contents of the stub
        stub_body = []
        stub_all: list[str] = []

        # Fill containers based on content of __init__.py
        for statement in init.body:

            # Regular import statements
            # e.g. "import numpy as np" or "import subprocess"
            if isinstance(statement, ast.Import):

                # Multiple items can be imported in a single statement
                for alias in statement.names:

                    # Regular imports from kernel are not supported
                    if "kernel" in alias.name:
                        raise ValueError(
                            f"Failed to generate {stub_path}: "
                            f"Regular imports from kernel are not supported. "
                            f"Requested: {ast.unparse(statement)}"
                        )

                    # External import
                    # TODO: Regular submodule imports should not be accepted
                    # TODO: How to check if an import is external?
                    stub_body.append(statement)
                    continue

            # ImportFrom statements
            # e.g. "from matplotlib import pyplot as plt", "from . import X"
            if isinstance(statement, ast.ImportFrom):

                # Submodule import from .
                # e.g. "from . import <submodule>"
                if statement.module is None:

                    # Check that statement follows rules
                    for alias in statement.names:

                        # Item should be a submodule
                        if not (module_path / alias.name).is_dir():
                            raise ValueError(
                                f"Failed to generate {stub_path}: "
                                f"Submodule {alias.name} not found in "
                                f"{module_path}"
                            )

                        # Submodule should not have an alias
                        if alias.asname is not None:
                            raise ValueError(
                                f"Failed to generate {stub_path}: "
                                f"Attempted to import submodule {alias.name} "
                                f"with alias {alias.asname}. "
                            )

                        # Add submodule to __all__
                        stub_all.append(alias.name)

                    # Add import statement to the stub
                    stub_body.append(statement)
                    continue

                # Relative import
                # e.g. "from .<submodule> import X"
                if statement.level > 0:

                    # Relative kernel imports are not supported
                    if "kernel" in statement.module:
                        raise ValueError(
                            f"Failed to generate {stub_path}: "
                            f"Relative kernel imports are not allowed. "
                            f"Requested: {ast.unparse(statement)}"
                        )

                    # Check for star import from submodule
                    if is_star_import(statement):

                        # Get submodule stubs path from statement
                        submodule_stubs_path = module_path / statement.module.replace(
                            ".", "/"
                        )

                        # Generate __init__.pyi if it does not exist
                        if not module_has_init_stub(submodule_stubs_path):
                            self.__generate_single_init_stub(submodule_stubs_path)

                        # Get contents in __all__ statement from __init__.pyi
                        submodule_all_statement = self.__find_all_statement(
                            submodule_stubs_path / "__init__.pyi"
                        )
                        submodule_contents = self.__retrieve_items_in_all(
                            submodule_all_statement
                        )

                        # If there is nothing to import, skip statement
                        if len(submodule_contents) == 0:
                            continue

                        # Replace star with imported items
                        __aliases = [ast.alias(item) for item in submodule_contents]
                        statement.names = __aliases

                    # Add all imported items, or their aliases, to __all__
                    for alias in statement.names:
                        if alias.asname is not None:
                            stub_all.append(alias.asname)
                        else:
                            stub_all.append(alias.name)

                    # Add import statement to the stub
                    stub_body.append(statement)
                    continue

                # Absolute import from kernel
                if "tudatpy.kernel" in statement.module:

                    if is_star_import(statement):

                        # Process star import
                        statement, stub_all = self.__expand_kernel_star_import(
                            statement, stub_all
                        )
                    else:

                        # Import from extension instead of kernel
                        if statement.module.split(".")[-1] == module_path.name:
                            statement.module = "extension"
                            statement.level = 1

                        # Add imported items to stub all
                        for item in statement.names:
                            stub_all.append(item.name)

                    # If the import statement is not empty, update stub body
                    if statement is not None:
                        stub_body.append(statement)
                    continue

                # External import statement
                stub_body.append(statement)

            # Assign statements
            if isinstance(statement, ast.Assign) or isinstance(
                statement, ast.AnnAssign
            ):

                # __all__ statement
                if "__all__ =" in ast.unparse(statement):

                    # Get items in __all__ statement
                    all_items = self.__retrieve_items_in_all(statement)

                    # Update __all__ list
                    stub_all.extend(all_items)
                    continue

                # __version__ statement
                if "__version__ =" in ast.unparse(statement):
                    stub_body.append(statement)
                    continue

                # Any other assign statement is unexpected
                raise NotImplementedError(
                    f"Failed to generate {stub_path}: "
                    f"Unexpected assign statement: "
                    f"{ast.unparse(statement)}"
                )

            # Other statements (We add them without modification)
            stub_body.append(statement)

        # Add __all__ to the stub
        if len(stub_all) > 0:
            all_stmt = ast.parse(
                "__all__ = [" + ", ".join([f"'{_name}'" for _name in stub_all]) + "]"
            ).body[0]
            stub_body.append(all_stmt)

        # Write stub to file
        stub.body = stub_body
        content = self.__unparse_script(stub)
        stub_path.write_text(content)

        return None

    def __generate_base_init_stub(self) -> None:

        # Get current list of submodules
        submodules: list[str] = []
        for dir in self.stubs_dir.iterdir():

            # Skip if not a directory or if it is a cache directory
            if not dir.is_dir() or dir.name == "__pycache__":
                continue

            # Skip if it is an ignored module
            if dir.name in self.ignored_modules:
                continue

            # Add submodule to list
            submodules.append(dir.name)

        # Initialize stub from contents of actual __init__.py
        self.__generate_single_init_stub(self.stubs_dir)

        # Update with submodules
        init = self.__parse_script(self.stubs_dir / "__init__.pyi")
        init_body = init.body

        # Generate import statement for __init__.pyi
        submodule_import = ast.ImportFrom(
            module="",
            names=[ast.alias(name) for name in submodules],
            level=1,
        )
        init_body.append(submodule_import)

        # Generate __all__ statement
        all_stmt = ast.parse(
            "__all__ = [" + ", ".join([f"'{_name}'" for _name in submodules]) + "]"
        ).body[0]
        init_body.append(all_stmt)

        # Update stub and remove extension.pyi
        init.body = init_body
        content = self.__unparse_script(init)
        (self.stubs_dir / "__init__.pyi").write_text(content)
        (self.stubs_dir / "extension.pyi").unlink()

        return None

    def __generate_init_stubs(self) -> None:

        # Find all the directories in the stubs directory
        for item in self.stubs_dir.rglob("*"):

            # Skip if not a directory or if it is a cache directory
            if not item.is_dir() or item.name == "__pycache__":
                continue

            # Create __init__.pyi stub for the module
            self.__generate_single_init_stub(item)

        # Generate __init__.pyi for the base module
        self.__generate_base_init_stub()

        return None

    def generate_stubs(self) -> None:

        # Create directory for stubs in build
        if self.stubs_dir.exists():
            shutil.rmtree(self.stubs_dir)
        self.stubs_dir.mkdir(exist_ok=False, parents=False)

        # Create directory directory structure for stubs
        self.__create_stubs_directory_structure()

        # Generate default stubs for tudatpy.kernel
        self.__generate_default_kernel_stubs()

        # Generate default stubs for Python source code
        self.__generate_default_python_stubs()

        # Fix autogenerated stubs
        self.__fix_autogenerated_stubs()

        # Missing: Fix __init__.pyi stubs to expand star imports
        self.__generate_init_stubs()

        print("Stubs generated successfully")

        return None


class Builder:

    def __init__(self, args: argparse.Namespace) -> None:

        # Command line arguments to control the build
        self.args = args

        # Resolve build directory
        self.build_dir = Path(self.args.build_dir).resolve()

        # Resolve conda prefix
        self.conda_prefix = Path(os.environ["CONDA_PREFIX"])
        if not self.conda_prefix.exists():
            raise FileNotFoundError(f"Conda prefix {self.conda_prefix} does not exist.")

        # Source directory of tudatpy
        self.python_source_dir = Path(__file__).parent / "src/tudatpy"
        # if not self.python_source_dir.exists():
        #     raise FileNotFoundError(
        #         f"Failed to generate stubs: "
        #         f"Source directory {self.python_source_dir} does not exist."
        #     )

        # Source directory for compiled extensions
        self.extension_source_dir = self.build_dir / "src/tudatpy"
        # if not self.extension_source_dir.exists():
        #     raise FileNotFoundError(
        #         f"Failed to generate stubs: "
        #         f"Source directory {self.extension_source_dir} does not exist."
        #     )

        # Configuration flags
        # self.skip_tudat = "OFF" if self.args.build_tudat else "ON"
        # self.skip_tudatpy = "OFF" if self.args.build_tudatpy else "ON"
        self.build_tests = "ON" if self.args.build_tests else "OFF"
        if self.args.build_github_actions:
            self.build_github_actions = "ON"
        else:
            self.build_github_actions = "OFF"

        return None

    @contextmanager
    def mock_environment(self) -> Generator[Environment, None, None]:

        try:
            # Generate temporary directory
            _tmp = tempfile.TemporaryDirectory()
            tmp = Path(_tmp.name)

            # Mock installation prefix in tmp
            mock_prefix = tmp / "mock_install"
            mock_prefix.mkdir(exist_ok=False, parents=False)

            # Temporary installation of tudatpy
            shutil.copytree(
                self.python_source_dir,
                mock_prefix / "tudatpy",
            )
            shutil.copy(
                self.extension_source_dir / "kernel.so",
                mock_prefix / "tudatpy/kernel.so",
            )

            # Create mock environment with tudatpy in PYTHONPATH
            mock_env = os.environ.copy()
            if "PYTHONPATH" not in mock_env:
                mock_env["PYTHONPATH"] = str(mock_prefix)
            else:
                mock_env["PYTHONPATH"] += f":{mock_prefix}"

            # Pack data and return
            yield Environment(mock_env, tmp, mock_prefix)

        finally:
            # Remove temporary directory
            if _tmp:
                _tmp.cleanup()

    def build_libraries(self) -> None:

        # If clean build is requested, delete build directory
        if self.build_dir.exists() and self.args.clean_build:

            uninstall_required: bool = False

            # Check for new manifests
            manifest_dir = self.build_dir / "manifests"
            if manifest_dir.is_dir():
                content = list(manifest_dir.iterdir())
                if len(content) != 0:
                    uninstall_required = True

            # Check for old manifests
            if (self.build_dir / "custom-manifest.txt").is_file():
                uninstall_required = True

            # If manifest present, ask to uninstall before
            if uninstall_required:
                print(
                    "WARNING\n"
                    "Installation manifests were found in the build "
                    "directory.\nPlease, remove all your current "
                    "installations of tudatpy\nwith the `uninstall.py`"
                    "script before running a clean build."
                )
                exit(1)

            # TODO: Logger
            print(f"Removing pre-existing build directory: {self.build_dir}")
            shutil.rmtree(self.build_dir)

        # If output to file is requested, set up output redirection
        if self.args.output_to_file:
            _output_dest = self.build_dir / "build_output.txt"
        else:
            _output_dest = None

        # If build directory still exists, skip setup
        if not self.build_dir.exists() or self.args.force_setup:

            # Set up build directory
            self.build_dir.mkdir(parents=True, exist_ok=True)
            with chdir(self.build_dir):

                # Base cmake command
                cmake_command = [
                    "cmake",
                    f"-DCMAKE_PREFIX_PATH={self.conda_prefix}",
                    f"-DCMAKE_INSTALL_PREFIX={self.conda_prefix}",
                    f"-DCMAKE_CXX_STANDARD={self.args.cxx_standard}",
                    "-DBoost_NO_BOOST_CMAKE=ON",
                    f"-DCMAKE_BUILD_TYPE={self.args.build_type}",
                    f"-DTUDAT_BUILD_TESTS={self.build_tests}",
                    f"-DTUDAT_BUILD_WITH_MCD={'ON' if self.args.build_with_mcd else 'OFF'}",
                    f"-DTUDAT_BUILD_GITHUB_ACTIONS={self.build_github_actions}",
                ]

                # Add Fortran compiler path if building with MCD
                if self.args.build_with_mcd:
                    gfortran_path = self.conda_prefix / "bin" / "gfortran"
                    if gfortran_path.exists():
                        cmake_command.append(
                            f"-DCMAKE_Fortran_COMPILER={gfortran_path}"
                        )
                    else:
                        print(
                            f"WARNING: gfortran not found at {gfortran_path}. "
                            "CMake will attempt to find it automatically."
                        )

                # Add build and source directories
                cmake_command.extend(["-B", f"{self.build_dir}", "-S", ".."])

                if _output_dest is None:
                    outcome = subprocess.run(cmake_command)
                else:
                    with _output_dest.open("w") as output_dest:
                        outcome = subprocess.run(
                            cmake_command,
                            stdout=output_dest,
                            stderr=output_dest,
                        )
            if outcome.returncode:
                shutil.rmtree(self.build_dir)
                raise RuntimeError("CMake setup failed")

        # Run build command
        with chdir(self.build_dir):

            build_command = ["cmake", "--build", ".", f"-j{self.args.j}"]
            if self.args.verbose:
                build_command.append("--verbose")
            outcome = subprocess.run(build_command)
            if outcome.returncode:
                exit(outcome.returncode)

        return None

    def generate_api_docs(self, mock_env: Environment) -> None:
        """Generate API documentation for the library

        :param mock_env: Mock environment with tudatpy installed
        """

        print("Building API documentation...")

        # Create directory for API docs in build
        api_docs_build_dir = Path(self.build_dir) / "api-docs"
        if api_docs_build_dir.exists():
            shutil.rmtree(api_docs_build_dir)
        api_docs_build_dir.mkdir(parents=True, exist_ok=False)

        # -E is added to force a clean API docs build, as changes in the compiled kernel are not detected otherwise
        api_docs_build_command = [
            "sphinx-build",
            "-b",
            "html",
            "docs/tudatpy/source",
            str(api_docs_build_dir),
            "-E",
        ]
        if self.args.verbose:
            api_docs_build_command.append("--verbose")

        outcome = subprocess.run(api_docs_build_command, env=mock_env.variables)
        if outcome.returncode:
            exit(outcome.returncode)

        return None

    def generate_stubs(self, mock_env: Environment) -> None:

        StubGenerator(self.build_dir, mock_env).generate_stubs()

        return None


if __name__ == "__main__":

    args = BuildParser().parse_args()
    builder = Builder(args)

    # Build libraries
    if not args.skip_build:
        builder.build_libraries()

    with builder.mock_environment() as mock_env:

        # Generate stubs
        if not args.skip_stubs:
            builder.generate_stubs(mock_env)

        # Generate API documentation
        if args.build_api_docs:
            builder.generate_api_docs(mock_env)
