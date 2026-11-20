#!/usr/bin/env python3
"""Build a standalone experiment against an existing Tudat CMake build."""

import argparse
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--build-dir", type=Path, default=Path("build-python-arcs"))
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
build = args.build_dir.resolve()
target = (
    build
    / "tests/test_tudat/src/astro/propagators/CMakeFiles/test_propagators_MultiArcVariationalEquations.dir"
)
flags = {}
for line in (target / "flags.make").read_text().splitlines():
    if " = " in line:
        key, value = line.split(" = ", 1)
        flags[key] = shlex.split(value)
link = shlex.split((target / "link.txt").read_text())
out = build / "benchmarks/buffer_reuse"
out.mkdir(parents=True, exist_ok=True)
sources = [
    root / "benchmarks/inner_loop/buffer_trial.cpp",
    root / "src/tudat/math/basic/legendrePolynomials.cpp",
]
objects = [out / (src.stem + ".o") for src in sources]
processes = []
for src, obj in zip(sources, objects):
    command = [
        link[0],
        *flags["CXX_DEFINES"],
        *flags["CXX_INCLUDES"],
        *flags["CXX_FLAGS"],
        "-fopenmp",
        "-DTUDAT_EXPERIMENT_PARALLEL_ACCELERATION",
        "-DTUDAT_EXPERIMENT_PARALLEL_LEGENDRE",
        "-g",
        "-c",
        str(src),
        "-o",
        str(obj),
    ]
    processes.append(subprocess.Popen(command))
statuses = [p.wait() for p in processes]
if any(statuses):
    raise SystemExit("Compilation failed: " + str(statuses))
rewritten = []
i = 0
while i < len(link):
    token = link[i]
    if token == "-o":
        rewritten += ["-o", str(out / "buffer_trial")]
        i += 2
        continue
    if token.startswith("-Wl,--dependency-file="):
        i += 1
        continue
    if token.endswith(".cpp.o"):
        rewritten += list(map(str, objects))
    else:
        rewritten.append(token)
    i += 1
subprocess.run([*rewritten, "-fopenmp", "-pthread"], cwd=target.parent.parent, check=True)
print(out / "buffer_trial")
