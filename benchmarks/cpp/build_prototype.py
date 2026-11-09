#!/usr/bin/env python3
"""Build against an existing local Tudat CMake build, using its test link configuration."""

import argparse
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--build-dir", type=Path, default=Path("build-parallel-arcs"))
parser.add_argument("--cereal-include", type=Path, help="directory containing cereal/access.hpp")
parser.add_argument("--source", type=Path, default=Path("benchmarks/cpp/parallelMultiArc.cpp"))
parser.add_argument("--name", default="parallel_multi_arc", help="output executable basename")
stages = parser.add_mutually_exclusive_group()
stages.add_argument("--compile-only", action="store_true")
stages.add_argument(
    "--link-only",
    action="store_true",
    help="reuse an object compiled with this build configuration",
)
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
out = build / "benchmarks"
out.mkdir(exist_ok=True)
source = args.source.resolve()
if Path(args.name).name != args.name:
    parser.error("--name must be a basename")
obj = out / (source.stem + ".o")
binary = out / args.name
command = [
    link[0],
    *flags["CXX_DEFINES"],
    *flags["CXX_INCLUDES"],
    *flags["CXX_FLAGS"],
    "-pthread",
    "-g",
    "-c",
    str(source),
    "-o",
    str(obj),
]
if args.cereal_include:
    command.extend(["-I", str(args.cereal_include.resolve())])
if not args.link_only:
    subprocess.run(command, check=True)
if args.compile_only:
    print(obj)
    raise SystemExit(0)
rewritten = []
i = 0
while i < len(link):
    token = link[i]
    if token == "-o":
        rewritten.extend(["-o", str(binary)])
        i += 2
        continue
    if token.startswith("-Wl,--dependency-file="):
        i += 1
        continue
    rewritten.append(str(obj) if token.endswith(".cpp.o") else token)
    i += 1
subprocess.run([*rewritten, "-pthread"], cwd=target.parent.parent, check=True)
print(binary)
