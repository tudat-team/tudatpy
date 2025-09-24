import argparse
from pathlib import Path
import shutil
from contextlib import chdir, contextmanager
import subprocess
import os
import ast
import tempfile
from dataclasses import dataclass
from typing import Generator, Callable
import logging
import enum
import types
import inspect
import re
import importlib
import sys


class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    blue = "\x1b[34;20m"
    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    _format = "%(levelname)-8s :: %(asctime)s :: %(message)s"

    FORMATS = {
        logging.DEBUG: green + _format + reset,
        logging.INFO: blue + _format + reset,
        logging.WARNING: yellow + _format + reset,
        logging.ERROR: red + _format + reset,
        logging.CRITICAL: bold_red + _format + reset,
    }

    def format(self, record) -> str:
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


log = logging.getLogger(__package__)
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
log.addHandler(ch)


class MethodTypes(enum.IntEnum):

    normal_method = enum.auto()
    class_method = enum.auto()
    static_method = enum.auto()


class CppEnumProcessor:

    def __init__(self, enum: type, module: "CppModuleProcessor") -> None:

        self.enum = enum
        self.name = enum.__name__
        self.docstring = enum.__doc__
        self.options = [
            attr
            for attr in dir(self.enum)
            if (not attr.startswith("_")) and (attr not in ("name", "value"))
        ]
        self.module = module

        # Set has_enum flag to true for module
        self.module.required_imports.add("enum")

        return None

    def process(self) -> ast.ClassDef:

        # Initialize empty body
        enum_body: list[ast.stmt] = []

        # Add docstring if present
        if self.docstring is not None:
            enum_body.append(ast.Expr(ast.Constant(self.docstring)))

        # Generate body for ClassDef
        for option in self.options:
            enum_body.append(
                ast.AnnAssign(
                    target=ast.Name(id=option, ctx=ast.Store()),
                    annotation=ast.Name(id=self.name, ctx=ast.Load()),
                    simple=1,
                )
            )

        return ast.ClassDef(
            name=self.name,
            bases=[ast.Name(id="enum.Enum", ctx=ast.Load())],
            keywords=[],
            body=enum_body,
            decorator_list=[],
        )


class CppFunctionProcessor:

    def __init__(
        self,
        name: str,
        function: Callable,
        module: "CppModuleProcessor",
        custom_decorators: list[str] | None = None,
    ) -> None:

        # Function name
        self.name = name
        self.__function = function
        self.processed_signature: str | None = None
        self.decorators: str | None = None
        self.docstring: str | None = None

        # Traceback and invalid flag
        self.invalid: bool = False
        self.warn: bool = False
        self.fatal: bool = False
        self.skip: bool = False
        self.traceback: str = ""

        # The module to which the function belongs
        self.module = module

        # Retrieve documentation for function
        is_python_function: bool = False
        if function.__doc__ is None:

            # Check if it is a Python function (e.g. expections)
            if isinstance(function, types.FunctionType):
                is_python_function = True
            # Otherwise fail due to lack of information about function
            else:
                self.invalid = True
                self.fatal = True
                self.traceback = (
                    self.module.cerror
                    + f"Missing documentation for {name} function ({function})"
                )
                return None

        # If python function, generate documentation
        if is_python_function:

            # Get signature from function object
            signature = inspect.signature(function)
            return_type = signature.return_annotation
            documentation = f"def {name}{str(signature)}"

            # Define return type if not specified
            if isinstance(return_type, inspect._empty):

                # If name is __init__, set return type to None
                if name == "__init__":
                    documentation += " -> None"
                else:
                    documentation += " -> typing.Any"
        else:

            # Sanity
            assert function.__doc__ is not None

            # Special case: constructor without signature
            if (
                self.name == "__init__"
                and "See help(type(self))" in function.__doc__
            ):
                documentation = (
                    "__init__(self) -> None\n\nCONSTRUCTOR NOT EXPOSED"
                )
                self.skip = True
            else:
                documentation = function.__doc__.strip()

        # Separate components of documentation
        self.__separate_documentation_components(
            documentation=documentation,
            custom_decorators=custom_decorators,
        )

        # Update import statements for module
        if not self.invalid:
            self.__update_required_imports()

        return None

    @property
    def documentation(self) -> str:

        common_traceback: str = (
            f"Failed to generate documentation for {self.name}: "
        )
        if self.decorators is None:
            log.fatal(common_traceback + "Missing decorators")
            exit(1)
        if self.signature is None:
            log.fatal(common_traceback + "Missing signature")
            exit(1)
        if self.docstring is None:
            log.fatal(common_traceback + "Missing docstring")
            exit(1)

        return f"{self.decorators}\n{self.signature}\n\t'''{self.docstring}'''".strip()

    def __update_required_imports(self) -> None:

        if self.signature is None:
            return None

        # Extract imported modules from signature with regex
        used_modules = re.findall(
            r"\s[a-zA-Z][a-zA-Z]*(?=\.)",
            self.signature,
        )
        # Update list of imported statements
        if len(used_modules) > 0:
            for item in used_modules:
                _item = str(item).strip()
                self.module.required_imports.add(_item)

        return None

    def __separate_documentation_components(
        self, documentation: str, custom_decorators: list[str] | None = None
    ) -> None:

        # Extract decorators if present
        __decorator_match = re.match(
            pattern=r"(@.*?\n?)*(?=def)",
            string=documentation,
        )
        decorators: str = (
            ""
            if custom_decorators is None
            else ("@" + "\n@".join(custom_decorators))
        )
        if __decorator_match is not None:
            decorators += "\n" + __decorator_match.group(0)
        self.decorators = decorators.strip()

        # Extract signature
        __signature_match = re.search(r"\(.+(?=\n|$)", documentation)
        if __signature_match is None:
            self.invalid = True
            self.fatal = True
            self.traceback = (
                self.module.cerror
                + f"Missing signature for {self.name} function"
            )
            return None
        signature = f"def {self.name}{str(__signature_match[0]).strip()}:"
        self.signature = self.__fix_signature(signature)

        __docstring_match = re.search(r"(?<=\n)(.|\n)+", documentation)
        self.docstring = "No documentation found"
        if __docstring_match is not None:
            self.docstring = __docstring_match.group(0).strip()

    def __fix_signature(self, signature: str) -> str:

        if self.invalid:
            return signature

        # Fail if C++ types in signature
        invalid_types_match = re.search(
            r"\w*(\:\:(\w+<.*[^-]>|\w+))+", signature
        )
        if invalid_types_match is not None:

            invalid_types = invalid_types_match.group(0).rstrip(",")
            self.invalid = True
            self.traceback = (
                self.module.cerror
                + f"Invalid type {invalid_types} in function {self.name} "
            )
            return signature

        # Detect invalid default values
        invalid_defaults_match = re.search(r"(?<=\= )(\W)*<.*?>", signature)
        if invalid_defaults_match is not None:
            invalid_defaults = invalid_defaults_match.group(0)
            self.invalid = True
            self.traceback = (
                self.module.cerror
                + f"Invalid default argument {invalid_defaults} in function {self.name} "
            )
            return signature

        # Replace long array typings
        signature = re.sub(
            r"typing\.Annotated\[numpy\.typing\.ArrayLike,.*?\[.*?\].*?\]",
            "numpy.ndarray",
            signature,
        )
        signature = re.sub(
            r"typing\.Annotated\[numpy\.typing\.NDArray\[.*?\].*?\[.*?\].*?\]",
            "numpy.ndarray",
            signature,
        )

        # Fix: Replace array with numpy.array
        signature = re.sub(r"(?<=[^.])array(?=\()", "numpy.array", signature)

        # Replace collections.abc with typing
        signature = re.sub(r"collections\.abc", "typing", signature)

        # Replace tudatpy.kernel with tudatpy
        signature = signature.replace("tudatpy.kernel", "tudatpy")

        # Replace arg0 with self (properties and setters)
        signature = signature.replace("(arg0", "(self")

        # Remove type annotation for self
        signature = re.sub(r"self:.*?(?=\)|,)", "self", signature)

        return signature

    def handle_errors(self, extra_traceback: str = "") -> None:

        if not self.invalid and not self.warn:
            return None

        traceback = self.traceback + extra_traceback
        if self.warn:
            log.warning(traceback)
            self.invalid = False
            return None

        if not self.fatal:
            log.error(traceback)
            return None

        log.fatal(traceback)
        exit(1)

    def process(self) -> ast.FunctionDef:

        try:
            fun = ast.parse(self.documentation).body[0]
        except SyntaxError as e:
            print(self.__function)
            print(f"Processed: {self.documentation}\n\n")
            raise e
        assert isinstance(fun, ast.FunctionDef)
        return fun


class CppClassProcessor:

    def __init__(self, class_: type, module: "CppModuleProcessor") -> None:

        self.class_ = class_

        # Documentation
        if self.class_.__doc__ is None:
            self.docstring = "No documentation found"
        else:
            self.docstring = class_.__doc__

        # Name and module
        self.name = class_.__name__
        self.module = module

        # Body
        self.body: list[ast.stmt] = [ast.Expr(ast.Constant(self.docstring))]

        return None

    def __get_type_of_method(
        self, method: "CppFunctionProcessor"
    ) -> MethodTypes:

        # Get name of first argument from signature
        first_argument_match = re.search(r"(?<=\()\w*", method.signature)

        # If no arguments, static method
        if first_argument_match is None:
            return MethodTypes.static_method

        match first_argument_match.group(0):
            case "self":
                return MethodTypes.normal_method
            case "cls":
                return MethodTypes.class_method
            case _:
                return MethodTypes.static_method

    def process(self) -> ast.ClassDef:

        # Get parents
        parents: list[str] = []
        for parent in self.class_.__bases__:
            if "pybind11_object" in str(parent):
                continue
            if parent.__module__ != self.class_.__module__:
                parents.append(parent.__module__ + parent.__name__)
            else:
                parents.append(parent.__name__)

        # Process the attributes of the class
        attributes: list[str] = ["__init__"]
        for x in dir(self.class_):
            if x.startswith("_") or x.endswith("_"):
                continue
            attributes.append(x)
        for attribute_name in attributes:

            # Get attribute
            attribute = getattr(self.class_, attribute_name)

            # Process methods
            if inspect.isroutine(attribute):

                # Initialize function processor
                method = CppFunctionProcessor(
                    name=attribute.__name__,
                    function=attribute,
                    module=self.module,
                    custom_decorators=None,
                )
                method.handle_errors(f"{self.name} class")

                # Adjust decorators based on type of method
                assert isinstance(method.decorators, str)
                match self.__get_type_of_method(method):

                    case MethodTypes.normal_method:
                        pass
                    case MethodTypes.class_method:
                        method.decorators += "@classmethod"
                    case MethodTypes.static_method:
                        method.decorators += "@staticmethod"

                # Check if method is static and fix
                if (not method.invalid) and (not method.skip):
                    self.body.append(method.process())

            # Process getters and setters
            elif hasattr(attribute, "fget") or hasattr(attribute, "fset"):

                # Process getter
                if (
                    hasattr(attribute, "fget")
                    and getattr(attribute, "fget") is not None
                ):

                    getter = CppFunctionProcessor(
                        name=attribute_name,
                        function=getattr(attribute, "fget"),
                        module=self.module,
                        custom_decorators=["property"],
                    )
                    getter.handle_errors(f"{self.name} class")
                    if not getter.invalid:
                        self.body.append(getter.process())

                # Process setter
                if (
                    hasattr(attribute, "fset")
                    and getattr(attribute, "fset") is not None
                ):

                    setter = CppFunctionProcessor(
                        name=attribute_name,
                        function=getattr(attribute, "fset"),
                        module=self.module,
                        custom_decorators=[f"{attribute_name}.setter"],
                    )
                    setter.handle_errors(f"{self.name} class")
                    if not setter.invalid:
                        self.body.append(setter.process())

            else:
                print(type(attribute), attribute)

        # Generate ClassDef object
        return ast.ClassDef(
            name=self.name,
            bases=[ast.Name(id=item, ctx=ast.Load()) for item in parents],
            keywords=[],
            body=self.body,
            decorator_list=[],
        )


class CppModuleProcessor:

    def __init__(
        self,
        module: types.ModuleType,
        output_directory: Path,
        ignored_submodules: list[str] | None = None,
    ) -> None:

        self.module = module
        self.name = module.__name__
        self.required_imports = set()
        self.syntax_tree = ast.Module(body=[], type_ignores=[])
        self.submodules: dict[str, "CppModuleProcessor"] = {}
        self.ignored_submodules = (
            ignored_submodules if ignored_submodules is not None else []
        )
        self.outdir: Path = output_directory.resolve()
        self.outdir.mkdir(exist_ok=True, parents=True)

        # Initialize container for items in __all__
        self.all_contents: list[str] = []

        # Common part of error message
        self.cerror = f"In {self.module.__name__}: "

        return None

    def __is_enum(self, item: type) -> bool:

        # Get list of non-special attributes
        non_special_attrs = [
            attr
            for attr in dir(item)
            if (not attr.startswith("_")) and (not attr.endswith("_"))
        ]

        # If name and value not present, not an enum
        if not (
            ("name" in non_special_attrs) and ("value" in non_special_attrs)
        ):
            return False

        # Get list of options
        options = [
            attr for attr in non_special_attrs if attr not in ("name", "value")
        ]

        # Ensure that all options are of the right type
        if not all(
            [isinstance(getattr(item, option), item) for option in options]
        ):
            return False

        # Raise warning if the enumeration exports the options
        if all([option in dir(self.module) for option in options]):
            log.warning(f"Enumeration {item} exports values")

        return True

    def process(self) -> None:

        log.info(f"Processing {self.name}")

        # Load contents of module
        components: list[str] = [
            x
            for x in dir(self.module)
            if (not x.startswith("_"))
            and (not x.endswith("_"))
            and (x not in self.ignored_submodules)
        ]

        # Process components based on type
        for component_name in components:

            # Load component
            component = getattr(self.module, component_name)

            # Modules and submodules
            if isinstance(component, types.ModuleType):
                self.submodules[component_name] = CppModuleProcessor(
                    module=component,
                    output_directory=self.outdir / component_name,
                )
                self.submodules[component_name].process()

            # Classes and enumerations
            elif isinstance(component, type):
                if self.__is_enum(component):

                    # Add item to syntax tree
                    self.syntax_tree.body.append(
                        CppEnumProcessor(component, self).process()
                    )

                    # Add item to __all__
                    self.all_contents.append(component_name)

                else:
                    class_ = CppClassProcessor(component, self).process()
                    if class_ is not None:

                        # Add item to syntax tree
                        self.syntax_tree.body.append(class_)

                        # Add item to __all__
                        self.all_contents.append(component_name)

            # Functions and methods from the kernel
            elif isinstance(component, types.BuiltinFunctionType):

                # Initialize function processor
                function_ = CppFunctionProcessor(
                    name=component.__name__,
                    function=component,
                    module=self,
                    custom_decorators=None,
                )
                function_.handle_errors()
                if not function_.invalid:

                    # Add item to syntax tree
                    self.syntax_tree.body.append(function_.process())

                    # Add item to __all__
                    self.all_contents.append(component_name)

            # Variables
            elif isinstance(
                component, int | float | complex | str | dict | tuple | list
            ):

                node = ast.parse(
                    f"{component_name}: {component.__class__.__name__}"
                ).body[0]
                assert isinstance(node, ast.AnnAssign)
                self.syntax_tree.body.append(node)

                # Add item to __all__
                self.all_contents.append(component_name)

            # Skip exported values of enumerations
            elif hasattr(component, "name") and getattr(
                component, "name"
            ) in dir(type(component)):
                pass

            # Not implemented
            else:
                log.fatal(f"INVALID: {component} :: {type(component)}")
                exit(1)

        # Define import statements
        import_statements: list[ast.stmt] = [
            ast.Import(names=[ast.alias(name)])
            for name in self.required_imports
        ]

        # Define __all__ statement
        all_statement_str = (
            '__all__ = ["' + '", "'.join(self.all_contents) + '"]'
        )
        all_statement = [ast.parse(all_statement_str).body[0]]

        # Update syntax tree with new statements
        updated_body: list[ast.stmt] = (
            import_statements + all_statement + self.syntax_tree.body
        )
        self.syntax_tree.body = updated_body

        # Add import statements to the body
        with open(f"{self.outdir}/extension.pyi", "w") as buffer:
            buffer.write(ast.unparse(self.syntax_tree))

        return None

    # def __parse_function_signature(self, signature: str) -> str:

    #     # If something is in brackets, put surround it by ""
    #     updated_sig = re.sub(r"<([^>]*)>", r'"<\1>"', signature)

    #     updated_sig = updated_sig.replace(
    #         "tudat::Time", "tudatpy.kernel.astro.time_representation.Time"
    #     )

    #     updated_sig = updated_sig.replace("tudatpy.kernel.", "")

    #     return updated_sig

    # def process_function(self, function: types.FunctionType) -> None:

    #     # Get name of the function
    #     name = function.__name__

    #     # Ensure that the function has documentation
    #     if function.__doc__ is None:
    #         log.error(self.cerror + f"Missing documentation for {name}")
    #         exit(1)
    #     description = function.__doc__.split("\n")

    #     # Ensure that signature is present in description
    #     signature: str | None = None
    #     docstring: str = ""
    #     for idx, line in enumerate(description):
    #         if line.split("(")[0] == name:
    #             signature = (
    #                 "\n".join(description[:idx]) + f"\ndef {line}:\n\t..."
    #             )
    #             docstring = "\n".join(description[idx + 1 :])
    #             break

    #     # Fail if signature is missing
    #     if signature is None:
    #         log.error(self.cerror + f"Signature missing for {name}")
    #         exit(1)

    #     x = ast.parse(signature).body[0]
    #     assert isinstance(x, ast.FunctionDef)

    #     x.body[0] = ast.Expr(ast.Constant(docstring))
    #     self.syntax_tree.body.append(x)

    #     # # Check if function signature is present in description
    #     # if description.split("(")[0] != name:
    #     #     log.error(self.cerror + f"Signature missing from doc for {name}")

    #     return None


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
            default="14",
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
        misc_group.add_argument(
            "--debug", action="store_true", help="Display debug information"
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


class PythonStatementProcessor:

    @staticmethod
    def process_statement(statement: ast.stmt) -> ast.stmt:

        if isinstance(statement, ast.Expr):
            return PythonStatementProcessor.process_expression_statement(
                statement
            )
        elif isinstance(statement, ast.Import):
            return PythonStatementProcessor.process_import(statement)
        elif isinstance(statement, ast.ImportFrom):
            return PythonStatementProcessor.process_import_from(statement)
        elif isinstance(statement, ast.ClassDef):
            return PythonStatementProcessor.process_classdef(statement)
        elif isinstance(statement, ast.FunctionDef):
            return PythonStatementProcessor.process_functiondef(statement)
        elif isinstance(statement, ast.AsyncFunctionDef):
            return PythonStatementProcessor.process_asyncfunctiondef(statement)
        elif isinstance(statement, ast.Assign):
            return PythonStatementProcessor.process_assign(statement)
        # elif isinstance(statement, ast.AnnAssign):
        #     return PythonStatementProcessor.process_annassign(statement)

        return statement

    @staticmethod
    def process_docstring(statement: ast.stmt) -> ast.Expr:

        # Fail if statement is not an expression
        if not isinstance(statement, ast.Expr):
            raise ValueError(
                f"Attempted to process {type(statement)} as docstring"
            )

        # Fail if value is not Constant
        if not isinstance(statement.value, ast.Constant):
            raise ValueError(
                f"Attempted to process {type(statement.value)} as docstring"
            )

        # Process docstring
        statement.value = PythonExpressionProcessor.process_constant(
            statement.value, docstring=True
        )

        return statement

    @staticmethod
    def process_expression_statement(expr: ast.Expr) -> ast.Expr:
        expr.value = PythonExpressionProcessor.process_expression(expr.value)
        return expr

    @staticmethod
    def process_import(statement: ast.Import) -> ast.Import:
        return statement

    @staticmethod
    def process_import_from(statement: ast.ImportFrom) -> ast.ImportFrom:
        return statement

    @staticmethod
    def process_classdef(statement: ast.ClassDef) -> ast.ClassDef:

        # Get class docstring
        docstring = ast.get_docstring(statement)

        # If docstring present, process it first
        if docstring is not None:
            statement.body[0] = PythonStatementProcessor.process_docstring(
                statement.body[0]
            )

        # Process the rest of the class normally
        for idx, stmt in enumerate(statement.body):
            statement.body[idx] = PythonStatementProcessor.process_statement(
                stmt
            )

        return statement

    @staticmethod
    def process_functiondef(statement: ast.FunctionDef) -> ast.FunctionDef:

        # Get function docstring
        docstring = ast.get_docstring(statement)

        # Remove body and process docstring if present
        if docstring is None:
            statement.body = [ast.Expr(ast.Constant("Missing docstring"))]
        else:
            # Fail if first statement in body is not docstring
            stmt = statement.body[0]
            if not isinstance(stmt, ast.Expr):
                raise ValueError(
                    f"Attempted to process {type(stmt)} as docstring"
                )

            statement.body = [PythonStatementProcessor.process_docstring(stmt)]

        return statement

    @staticmethod
    def process_asyncfunctiondef(
        statement: ast.AsyncFunctionDef,
    ) -> ast.AsyncFunctionDef:

        # Get function docstring
        docstring = ast.get_docstring(statement)

        # Remove body and process docstring if present
        if docstring is None:
            statement.body = []
        else:
            statement.body = [
                PythonStatementProcessor.process_docstring(statement.body[0])
            ]

        return statement

    @staticmethod
    def process_assign(statement: ast.Assign) -> ast.Assign:

        log.warning(
            f"Missing type annotation in assignment statement: {statement.lineno}"
        )

        return statement

    @staticmethod
    def process_annassign(statement: ast.AnnAssign) -> ast.AnnAssign:

        # Remove value from statement
        statement.value = None

        return statement


class PythonExpressionProcessor:

    @staticmethod
    def process_expression(expr: ast.expr) -> ast.expr:

        if isinstance(expr, ast.Constant):
            return PythonExpressionProcessor.process_constant(expr)

        return expr

    @staticmethod
    def process_constant(
        const: ast.Constant, docstring: bool = False
    ) -> ast.Constant:

        # Skip processing if not a docstring
        if not docstring:
            return const

        # Fail if not a string (should not happen)
        if not isinstance(const.value, str):
            raise NotImplementedError(
                f"Attempted to process non-string constant: {const.lineno}"
            )

        # Extract components of docstring
        components = [x.strip() for x in const.value.splitlines()]

        # Remove empty lines in the beginning
        first_line = components[0]
        while first_line == "":
            components.pop(0)
            first_line = components[0]

        # Remove empty lines in the end
        last_line = components[-1]
        while last_line == "":
            components.pop(-1)
            last_line = components[-1]

        # Fix indentation
        indentation = " " * const.col_offset
        const.value = (f"\n{indentation}").join(components) + f"\n{indentation}"

        return const


class StubGenerator:

    # Default indentation length in pybind11-stubgen
    indentation: str = " " * 4

    # Ignored modules and methods
    ignored_modules: list[str] = [
        "temp",
        "io",
        "numerical_simulation",
        "exceptions",
    ]
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

    def __fix_tudatpy_imports(
        self, module: ast.Module, stub: Path
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
                internal_imports: list[str] = []
                external_imports: list[ast.alias] = []
                for alias in statement.names:

                    if "tudatpy.kernel." not in alias.name:
                        external_imports.append(alias)
                    else:

                        # Remove tudatpy.kernel. prefix from import
                        relative_alias = alias.name.replace(
                            "tudatpy.kernel.", ""
                        )
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
                if (
                    ast.unparse(statement)
                    == "from __future__ import annotations"
                ):
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

        # Process module-level docstring if present
        module_docstring = ast.get_docstring(module)
        if module_docstring is not None:
            module.body[0] = PythonStatementProcessor.process_docstring(
                module.body[0]
            )

        # Process the whole module with normal rules
        for idx, statement in enumerate(module.body):
            module.body[idx] = PythonStatementProcessor.process_statement(
                statement
            )

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

    def __generate_kernel_stubs(self) -> None:

        # Generate stubs for kernel
        sys.path.append(str(self.mock_env.prefix))
        CppModuleProcessor(
            module=importlib.import_module("tudatpy.kernel"),
            output_directory=self.stubs_dir,
            ignored_submodules=self.ignored_modules,
        ).process()
        sys.path.pop(-1)

        return None

    def __generate_default_kernel_stubs(self) -> None:

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
            capture_output=True,
        )

        # Interrupt if stub-generation failed
        if outcome.returncode:

            # Show stderr as error
            for stderr_line in str(outcome.stderr)[2:-1].split(r"\n"):
                log.error(stderr_line)

            log.error("Failed to generate stubs for tudatpy kernel")
            exit(1)

        # Show stub generation stderr as debug
        for stderr_line in str(outcome.stderr)[2:-1].split(r"\n"):
            log.debug(stderr_line)

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
                    self.stubs_dir
                    / relative_path.with_suffix("")
                    / "extension.pyi",
                )

        return None

    def __generate_default_python_stubs(self) -> None:

        log.info("Generating stubs for Python source code")

        # Process all the python scripts in the library
        mock_installation_dir = self.mock_env.prefix / "tudatpy"
        for script in mock_installation_dir.rglob("*.py"):

            # Skip if __init__.py or part of __pycache__
            if (script.name == "__init__.py") or ("__pycache__" in str(script)):
                continue

            # Skip if part of ignored module
            if (
                script.relative_to(mock_installation_dir).parts[0]
                in self.ignored_modules
            ):
                continue

            # Generate syntax tree from script contents
            script_content = self.__parse_script(script)

            # Initialize syntax tree for stub
            stub_content = ast.Module(body=[], type_ignores=[])

            # Process module-level docstring if present
            module_docstring = ast.get_docstring(script_content)
            if module_docstring is not None:
                script_content.body[0] = (
                    PythonStatementProcessor.process_docstring(
                        script_content.body[0]
                    )
                )

            # Fill stub body with script contents
            for statement in script_content.body:
                stub_content.body.append(
                    PythonStatementProcessor.process_statement(statement)
                )

            # Define path to stub and generate parent directory if missing
            stub_path = self.stubs_dir / script.relative_to(
                mock_installation_dir
            ).with_suffix(".pyi")
            # stub_path.parent.mkdir(exist_ok=True, parents=True)

            # Generate stub file from syntax tree
            with stub_path.open("w") as buffer:
                buffer.write(self.__unparse_script(stub_content))

        #     # Generate stub with stubgen
        #     outcome = subprocess.run(
        #         [
        #             "stubgen",
        #             script,
        #             "-o",
        #             str(python_stubs_dir),
        #             "--include-docstrings",
        #             "--parse-only",
        #             "--quiet",
        #         ]
        #     )
        #     if outcome.returncode:
        #         raise RuntimeError(f"Failed to generate stub for {script}")

        # # Move stubs to build directory
        # for stub in (python_stubs_dir / "tudatpy").rglob("*.pyi"):
        #     shutil.copy(
        #         stub,
        #         self.stubs_dir / stub.relative_to(python_stubs_dir / "tudatpy"),
        #     )

        # # Move stubs to build directory
        # print("The normal ones")
        # for stub in (python_stubs_dir / "tudatpy").rglob("*.pyi"):

        #     print(stub)
        #     shutil.copy(
        #         stub,
        #         self.stubs_dir / stub.relative_to(python_stubs_dir / "tudatpy"),
        #     )

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
            if not (
                isinstance(_name, ast.Constant) and isinstance(_name.value, str)
            ):
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
            raise ValueError(
                f"Script {script} does not contain an __all__ statement."
            )

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
                        submodule_stubs_path = (
                            module_path / statement.module.replace(".", "/")
                        )

                        # Generate __init__.pyi if it does not exist
                        if not module_has_init_stub(submodule_stubs_path):
                            self.__generate_single_init_stub(
                                submodule_stubs_path
                            )

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
                        __aliases = [
                            ast.alias(item) for item in submodule_contents
                        ]
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

                    # Non-star imports from kernel are not allowed
                    if not (
                        len(statement.names) == 1
                        and statement.names[0].name == "*"
                    ):
                        raise ValueError(
                            f"Failed to generate {stub_path}: "
                            f"Only star imports from kernel are supported. "
                            f"Requested: {ast.unparse(statement)}"
                        )

                    # Process star import
                    statement, stub_all = self.__expand_kernel_star_import(
                        statement, stub_all
                    )

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
                "__all__ = ["
                + ", ".join([f"'{_name}'" for _name in stub_all])
                + "]"
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
            "__all__ = ["
            + ", ".join([f"'{_name}'" for _name in submodules])
            + "]"
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
        # self.__generate_default_kernel_stubs()
        self.__generate_kernel_stubs()

        # Fix autogenerated stubs
        self.__fix_autogenerated_stubs()

        # Generate default stubs for Python source code
        self.__generate_default_python_stubs()

        # Missing: Fix __init__.pyi stubs to expand star imports
        self.__generate_init_stubs()

        # Display success message
        log.info("Stubs generated successfully")

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
            raise FileNotFoundError(
                f"Conda prefix {self.conda_prefix} does not exist."
            )

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

                if _output_dest is None:
                    outcome = subprocess.run(
                        [
                            "cmake",
                            # f"-DSKIP_TUDAT={self.skip_tudat}",
                            # f"-DSKIP_TUDATPY={self.skip_tudatpy}",
                            f"-DCMAKE_PREFIX_PATH={self.conda_prefix}",
                            f"-DCMAKE_INSTALL_PREFIX={self.conda_prefix}",
                            f"-DCMAKE_CXX_STANDARD={self.args.cxx_standard}",
                            "-DBoost_NO_BOOST_CMAKE=ON",
                            f"-DCMAKE_BUILD_TYPE={self.args.build_type}",
                            f"-DTUDAT_BUILD_TESTS={self.build_tests}",
                            f"-DTUDAT_BUILD_GITHUB_ACTIONS={self.build_github_actions}",
                            "-B",
                            f"{self.build_dir}",
                            "-S",
                            "..",
                        ]
                    )
                else:
                    with _output_dest.open("w") as output_dest:
                        outcome = subprocess.run(
                            [
                                "cmake",
                                # f"-DSKIP_TUDAT={self.skip_tudat}",
                                # f"-DSKIP_TUDATPY={self.skip_tudatpy}",
                                f"-DCMAKE_PREFIX_PATH={self.conda_prefix}",
                                f"-DCMAKE_INSTALL_PREFIX={self.conda_prefix}",
                                f"-DCMAKE_CXX_STANDARD={self.args.cxx_standard}",
                                "-DBoost_NO_BOOST_CMAKE=ON",
                                f"-DCMAKE_BUILD_TYPE={self.args.build_type}",
                                f"-DTUDAT_BUILD_TESTS={self.build_tests}",
                                f"-DTUDAT_BUILD_GITHUB_ACTIONS={self.build_github_actions}",
                                "-B",
                                f"{self.build_dir}",
                                "-S",
                                "..",
                            ],
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

    # Parse command-line arguments
    args = BuildParser().parse_args()

    # Set logging level
    if args.debug:
        log.setLevel(logging.DEBUG)

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
