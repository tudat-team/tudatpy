"""Check the source dependency boundary without importing or building TudatPy.

Run with ``python -m unittest discover -s tests/architecture -v``.
These checks traverse local Python imports and C++ includes (including the
implementation corresponding to an included header).  Compiled binding smoke
tests additionally check registration and linking; this is not a C++ parser.
"""

import ast
from collections import deque
import importlib.util
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "src"
DATA_INPUT = PYTHON_ROOT / "tudatpy/data_input"
FORBIDDEN_MODULES = (
    "tudatpy.dynamics",
    "tudatpy.estimation",
    "tudatpy.numerical_simulation",
    "tudatpy.kernel.dynamics",
    "tudatpy.kernel.estimation",
    "tudatpy.kernel.numerical_simulation",
)
FORBIDDEN_CPP = (
    "tudat/simulation/",
    "tudat/astro/orbit_determination/",
    "tudat/astro/propagators/",
    "tudat/astro/ephemerides/",
    "tudatpy/dynamics/",
    "tudatpy/estimation/",
    "tudatpy/numerical_simulation/",
)
INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]([^">]+)[">]', re.MULTILINE)


def module_name(path):
    parts = list(path.relative_to(PYTHON_ROOT).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def imported_modules(path):
    name = module_name(path)
    package = name if path.name == "__init__.py" else name.rpartition(".")[0]
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"), filename=str(path))):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported = "." * node.level + (node.module or "")
            base = importlib.util.resolve_name(imported, package) if node.level else imported
            yield base
            for alias in node.names:
                if alias.name != "*":
                    yield base + "." + alias.name
        elif isinstance(node, ast.Call) and node.args:
            function = node.func
            is_dynamic_import = (
                isinstance(function, ast.Name) and function.id == "__import__"
            ) or (isinstance(function, ast.Attribute) and function.attr == "import_module")
            if is_dynamic_import and isinstance(node.args[0], ast.Constant):
                value = node.args[0].value
                if isinstance(value, str) and not value.startswith("."):
                    yield value


def local_python_files(name):
    # Package initialisers also execute when a submodule is imported.
    parts = name.split(".")
    for depth in range(1, len(parts) + 1):
        base = PYTHON_ROOT.joinpath(*parts[:depth])
        initialiser = base / "__init__.py"
        if initialiser.is_file():
            yield initialiser
    module = PYTHON_ROOT.joinpath(*parts).with_suffix(".py")
    if module.is_file():
        yield module


def python_dependency_violations():
    queue = deque((path, [module_name(path)]) for path in sorted(DATA_INPUT.rglob("*.py")))
    visited = set()
    violations = []
    while queue:
        path, chain = queue.popleft()
        if path in visited:
            continue
        visited.add(path)
        for name in imported_modules(path):
            if any(name == prefix or name.startswith(prefix + ".") for prefix in FORBIDDEN_MODULES):
                violations.append(" -> ".join(chain + [name]))
            elif name == "tudatpy" or name.startswith("tudatpy."):
                for imported_path in local_python_files(name):
                    queue.append((imported_path, chain + [module_name(imported_path)]))
    return violations


def cpp_dependency_violations():
    roots = sorted(DATA_INPUT.rglob("*.cpp")) + sorted(DATA_INPUT.rglob("*.h"))
    # Include the standalone parser directly, even if its binding registration changes.
    roots += [
        ROOT / "include/tudat/io/readAtmosphericTideData.h",
        ROOT / "src/tudat/io/readAtmosphericTideData.cpp",
        ROOT / "include/tudat/io/readStationLoadingData.h",
        ROOT / "src/tudat/io/readStationLoadingData.cpp",
        ROOT / "include/tudat/io/readTerrestrialFrameTransformationData.h",
        ROOT / "src/tudat/io/readTerrestrialFrameTransformationData.cpp",
    ]
    queue = deque((path, [str(path.relative_to(ROOT))]) for path in roots)
    visited = set()
    violations = []
    while queue:
        path, chain = queue.popleft()
        if path in visited:
            continue
        visited.add(path)
        relative = path.relative_to(ROOT)
        if relative.parts[0] in ("include", "src"):
            logical_path = "/".join(relative.parts[1:])
            if logical_path.startswith(FORBIDDEN_CPP):
                violations.append(" -> ".join(chain))
                continue
        if relative.parts[0] == "include":
            implementation = ROOT / "src" / Path(*relative.parts[1:]).with_suffix(".cpp")
            if implementation.is_file():
                queue.append((implementation, chain + [str(implementation.relative_to(ROOT))]))
        for included in INCLUDE.findall(path.read_text(encoding="utf-8")):
            candidates = (
                path.parent / included,
                ROOT / "include" / included,
                ROOT / "src/tudatpy" / included,
            )
            for candidate in candidates:
                if candidate.is_file():
                    candidate = candidate.resolve()
                    queue.append((candidate, chain + [str(candidate.relative_to(ROOT))]))
                    break
    return violations


class DataInputDependencies(unittest.TestCase):
    def test_python_has_no_dynamics_or_estimation_imports(self):
        self.assertTrue(DATA_INPUT.is_dir(), "Run against a complete source checkout")
        self.assertEqual(python_dependency_violations(), [])

    def test_cpp_has_no_dynamics_or_estimation_includes(self):
        self.assertTrue(DATA_INPUT.is_dir(), "Run against a complete source checkout")
        self.assertEqual(cpp_dependency_violations(), [])


if __name__ == "__main__":
    unittest.main()
