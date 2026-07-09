#!/usr/bin/env python3
"""
Serialization Compliance Tracker

Scans the Tudat codebase and produces an Excel report showing:
  1. What "Settings" classes exist
  2. What classes implement load/save for serialization
  3. What classes implement operator==
  4. What classes have C++ roundtrip testing
  5. What classes have Python exposure for: Equality, Pickle, save_to_*
  6. What classes have Python roundtrip testing

Usage:
    python track_serialization.py [--csv] [--output PATH]

Output: track_serialization.xlsx (default) or CSV if --csv is given.
"""

import re
import argparse
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field

# ── workspace root ──────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent
INCLUDE_DIR = ROOT / "include" / "tudat"
SRC_TUDATPY_DIR = ROOT / "src" / "tudatpy"
TESTS_CPP_DIR = ROOT / "tests" / "test_tudat" / "src" / "io"
TESTS_PY_DIR = ROOT / "tests" / "test_tudatpy"

DEFAULT_EXCEL = ROOT / "track_serialization.xlsx"


# ── data structures ─────────────────────────────────────────────────────────


@dataclass
class ClassInfo:
    name: str
    is_settings: bool = False
    header_file: str = ""
    has_save_load: bool = False
    has_operator_eq: bool = False
    has_equals_method: bool = False
    base_classes: list = field(default_factory=list)  # direct base class names
    cpp_test_files: list = field(default_factory=list)
    py_expose_file: str = ""
    py_has_equals: bool = False
    py_has_pickle: bool = False
    py_has_save_to: bool = False
    py_test_files: list = field(default_factory=list)


# ── helpers ──────────────────────────────────────────────────────────────────


def strip_template(name: str) -> str:
    """Strip template params and namespace: 'Foo<T>' or 'ns::Foo' → 'Foo'"""
    base = name.split("<")[0].strip()
    if "::" in base:
        base = base.split("::")[-1]
    return base


def is_settings_class(name: str) -> bool:
    return name.lower().endswith("settings") and name != "Settings"


# ── Phase 1: Scan headers for class inventory ───────────────────────────────


def scan_headers() -> dict[str, ClassInfo]:
    """Build inventory of all interesting classes from headers.

    For each header file, find all class declarations, then check
    per-class whether the class body contains save/load pairs and
    operator==, and which classes are settings classes.
    """
    results: dict[str, ClassInfo] = {}

    # Patterns for class/struct bodies: capture class name and position of '{'
    class_start_re = re.compile(
        r"(?:^|\n)\s*(?:class|struct)\s+(\w+)\s+[^{;]{0,2000}(\{)", re.MULTILINE
    )
    save_re = re.compile(r"void\s+save\s*\(\s*[^)]*Archive[^)]*\)")
    load_re = re.compile(r"void\s+load\s*\(\s*[^)]*Archive[^)]*\)")
    operator_eq_re = re.compile(r"(?:friend\s+)?(?:bool|auto)\s+operator\s*==\s*\([^)]*\)")
    equals_method_re = re.compile(r"bool\s+equals\s*\([^)]*\)\s*const(?:\s+override)?")

    for hdr in sorted(INCLUDE_DIR.rglob("*.h")):
        rel = str(hdr.relative_to(ROOT))
        try:
            content = hdr.read_text(errors="ignore")
        except Exception:
            continue

        # Find each class/struct with its opening brace position
        # Build list of (cls_name, body_start, body_end) where body
        # is the text between the opening '{' and matching '}'.
        classes = []
        for m in class_start_re.finditer(content):
            cls_name = m.group(1)
            # Capture inheritance from the declaration line
            decl_line = content[m.start() : m.start(2)]
            base_re = re.compile(r":\s*(?:public|protected|private)\s+(\w+)")
            base_m = base_re.search(decl_line)
            base_name = base_m.group(1) if base_m else ""
            brace_start = m.start(2)  # position of '{'
            # Walk forward to find matching '}' (tracking brace depth)
            depth = 0
            pos = brace_start
            while pos < len(content):
                ch = content[pos]
                if ch == "{":
                    depth += 1
                elif ch == "}":
                    depth -= 1
                    if depth == 0:
                        break
                pos += 1
            if depth != 0:
                # Unmatched braces — skip
                continue
            classes.append((cls_name, base_name, brace_start, pos + 1))

        for cls_name, base_name, body_start, body_end in classes:
            # Extract class body text
            body = content[body_start:body_end]

            is_set = is_settings_class(cls_name)
            has_save = bool(save_re.search(body))
            has_load = bool(load_re.search(body))
            has_op_eq = bool(operator_eq_re.search(body))
            has_eq_method = bool(equals_method_re.search(body))

            if not (is_set or has_save or has_load or has_op_eq or has_eq_method):
                continue

            if cls_name not in results:
                results[cls_name] = ClassInfo(
                    name=cls_name,
                    is_settings=is_set,
                    header_file=rel,
                )
            info = results[cls_name]
            if is_set:
                info.is_settings = True
            if has_save and has_load:
                info.has_save_load = True
            if has_op_eq:
                info.has_operator_eq = True
            if has_eq_method:
                info.has_equals_method = True
            if base_name:
                info.base_classes.append(base_name)

    # Propagate has_operator_eq down the inheritance chain.
    # operator== is inherited — if a base class has it, all derived classes do too.
    changed = True
    while changed:
        changed = False
        for cls_name, info in results.items():
            if info.has_operator_eq:
                continue
            for base in info.base_classes:
                if base in results and results[base].has_operator_eq:
                    info.has_operator_eq = True
                    changed = True
                    break

    return results


# ── Phase 2: C++ tests ────────────────────────────────────────────────────


def scan_cpp_tests(class_infos: dict[str, ClassInfo]) -> None:
    """Map serialization test files to the classes they test."""
    if not TESTS_CPP_DIR.exists():
        return

    for test_file in sorted(TESTS_CPP_DIR.glob("*.cpp")):
        if "serial" not in test_file.name.lower():
            continue
        try:
            content = test_file.read_text(errors="ignore")
        except Exception:
            continue
        rel = str(test_file.relative_to(ROOT))

        for cls_name, info in class_infos.items():
            if not info.has_save_load:
                continue
            base = strip_template(cls_name)
            if re.search(rf"\b{re.escape(base)}\b", content):
                info.cpp_test_files.append(rel)


# ── Phase 3: Python exposure ──────────────────────────────────────────────


def _find_py_bindings(content: str) -> list[tuple[str, int, int]]:
    """Return [(short_class_name, start, template_end), ...] for py::class_<...>."""
    results = []
    for m in re.finditer(r"py::class_\s*<", content):
        pos = m.end()
        depth = 0
        end = pos
        while end < len(content):
            c = content[end]
            if c == "<":
                depth += 1
            elif c == ">":
                if depth == 0:
                    break
                depth -= 1
            end += 1
        template = content[pos:end]

        # Get first argument (split on commas at depth 0)
        first_arg = ""
        d = 0
        for c in template:
            if c == "," and d == 0:
                break
            first_arg += c
            if c == "<":
                d += 1
            elif c == ">":
                d -= 1
        first_arg = first_arg.strip()

        # class name: last identifier before '<' or end
        lt = first_arg.find("<")
        name_part = first_arg[:lt] if lt >= 0 else first_arg
        nm = re.search(r"(\w+)\s*$", name_part)
        if nm:
            results.append((nm.group(1), m.start(), end))
    return results


def scan_python_exposure(class_infos: dict[str, ClassInfo]) -> None:
    """Detect __eq__, pickle, save_to in expose*.cpp files."""
    if not SRC_TUDATPY_DIR.exists():
        return

    # Build lookup: short_name → [cls_name, ...]
    name_map: dict[str, list[str]] = defaultdict(list)
    for cls_name in class_infos:
        name_map[strip_template(cls_name)].append(cls_name)

    for exp_file in sorted(SRC_TUDATPY_DIR.rglob("expose*.cpp")):
        try:
            content = exp_file.read_text(errors="ignore")
        except Exception:
            continue
        rel = str(exp_file.relative_to(ROOT))

        bindings = _find_py_bindings(content)

        for i, (short_name, _start, template_end) in enumerate(bindings):
            if short_name not in name_map:
                continue

            # Segment from template end to next binding start (or EOF)
            if i + 1 < len(bindings):
                seg_end = bindings[i + 1][1]
            else:
                seg_end = len(content)
            segment = content[template_end:seg_end]

            has_eq = bool(
                re.search(r'"__eq__"', segment)
                or re.search(r"py::self\s*==", segment)
                or re.search(r"TUDATPY_DEF_EQ_NE\s*\(", segment)
            )
            has_pickle = (
                "py::pickle" in segment
                or "make_pickle" in segment
                or bool(re.search(r"TUDATPY_DEF_PICKLE", segment))
            )
            has_save = bool(
                re.search(r'"save_(?:to_)?(?:json|binary|bin)"', segment)
                or re.search(r'"load_(?:from_)?(?:json|binary|bin)"', segment)
            )

            for cls_name in name_map[short_name]:
                info = class_infos[cls_name]
                info.py_expose_file = rel
                if has_eq:
                    info.py_has_equals = True
                if has_pickle:
                    info.py_has_pickle = True
                if has_save:
                    info.py_has_save_to = True


# ── Phase 4: Python tests ─────────────────────────────────────────────────


def scan_python_tests(class_infos: dict[str, ClassInfo]) -> None:
    """Detect Python roundtrip tests."""
    if not TESTS_PY_DIR.exists():
        return

    for py_test in sorted(TESTS_PY_DIR.rglob("*.py")):
        try:
            content = py_test.read_text(errors="ignore")
        except Exception:
            continue
        lower = content.lower()
        if "pickle" not in lower and "serial" not in lower:
            continue

        rel = str(py_test.relative_to(ROOT))

        for cls_name, info in class_infos.items():
            base = strip_template(cls_name)
            if re.search(rf"\b{re.escape(base)}\b", content):
                info.py_test_files.append(rel)


# ── Excel / CSV output ────────────────────────────────────────────────────


def _yesno(b: bool) -> str:
    return "✓" if b else "✗"


def write_excel(class_infos: dict[str, ClassInfo], output_path: Path) -> Path:
    try:
        import openpyxl
    except ImportError:
        print("openpyxl not installed. Install with: pip install openpyxl")
        return write_csv(class_infos, output_path.with_suffix(".csv"))

    from openpyxl.styles import Font, PatternFill, Alignment

    wb = openpyxl.Workbook()
    hdr_font = Font(bold=True, size=11, color="FFFFFF")
    hdr_fill = PatternFill("solid", fgColor="4472C4")
    yes_fill = PatternFill("solid", fgColor="C6EFCE")
    no_fill = PatternFill("solid", fgColor="FFC7CE")
    settings_fill = PatternFill("solid", fgColor="FFEB9C")
    center = Alignment(horizontal="center", vertical="top")
    wrap = Alignment(wrap_text=True, vertical="top")

    # ── Sheet 1: All Classes ──
    ws = wb.active
    ws.title = "Serialization Status"

    headers = [
        "Class",
        "Is Settings",
        "Save/Load",
        "operator==",
        "equals() method",
        "C++ Roundtrip Test",
        "Python Exposure File",
        "Py Binding: operator==",
        "Py Binding: pickle",
        "Py Binding: file I/O (save/load)",
        "Python Roundtrip Test",
        "Header File",
        "C++ Test Files",
        "Python Test Files",
    ]

    for c, h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font = hdr_font
        cell.fill = hdr_fill
        cell.alignment = wrap

    row = 2
    for cls_name in sorted(class_infos.keys()):
        info = class_infos[cls_name]
        vals = [
            cls_name,
            _yesno(info.is_settings),
            _yesno(info.has_save_load),
            _yesno(info.has_operator_eq),
            _yesno(info.has_equals_method),
            _yesno(bool(info.cpp_test_files)),
            info.py_expose_file or "",
            _yesno(info.py_has_equals),
            _yesno(info.py_has_pickle),
            _yesno(info.py_has_save_to),
            _yesno(bool(info.py_test_files)),
            info.header_file,
            "; ".join(info.cpp_test_files),
            "; ".join(info.py_test_files),
        ]
        for c, v in enumerate(vals, 1):
            cell = ws.cell(row=row, column=c, value=v)
            cell.alignment = wrap
            if c in (2, 3, 4, 5, 6, 8, 9, 10, 11):
                cell.alignment = center
                if v == "✓":
                    cell.fill = yes_fill
                elif v == "✗":
                    cell.fill = no_fill
        if info.is_settings:
            ws.cell(row=row, column=1).fill = settings_fill
        row += 1

    col_widths = [40, 12, 12, 12, 16, 20, 50, 24, 22, 30, 22, 55, 55, 55]
    for i, w in enumerate(col_widths, 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{openpyxl.utils.get_column_letter(len(headers))}{row - 1}"

    # ── Sheet 2: Summary ──
    ws2 = wb.create_sheet("Summary")
    total = len(class_infos)
    settings = sum(1 for i in class_infos.values() if i.is_settings)

    def _c(attr, so=False):
        if so:
            return sum(1 for i in class_infos.values() if i.is_settings and getattr(i, attr))
        return sum(1 for i in class_infos.values() if getattr(i, attr))

    def _cb(attr, so=False):
        if so:
            return sum(1 for i in class_infos.values() if i.is_settings and bool(getattr(i, attr)))
        return sum(1 for i in class_infos.values() if bool(getattr(i, attr)))

    for r, row_data in enumerate(
        [
            ("Metric", "All Classes", "Settings Only"),
            ("Total", total, settings),
            ("Has save/load", _c("has_save_load"), _c("has_save_load", True)),
            ("Has operator==", _c("has_operator_eq"), _c("has_operator_eq", True)),
            ("Has equals() method", _c("has_equals_method"), _c("has_equals_method", True)),
            ("C++ roundtrip test", _cb("cpp_test_files"), _cb("cpp_test_files", True)),
            ("Py Binding: operator==", _c("py_has_equals"), _c("py_has_equals", True)),
            ("Py Binding: pickle", _c("py_has_pickle"), _c("py_has_pickle", True)),
            ("Py Binding: file I/O (save/load)", _c("py_has_save_to"), _c("py_has_save_to", True)),
            ("Python roundtrip test", _cb("py_test_files"), _cb("py_test_files", True)),
        ],
        1,
    ):
        for c, val in enumerate(row_data, 1):
            cell = ws2.cell(row=r, column=c, value=val)
            if r == 1:
                cell.font = hdr_font
                cell.fill = hdr_fill
    ws2.column_dimensions["A"].width = 30
    ws2.column_dimensions["B"].width = 16
    ws2.column_dimensions["C"].width = 16

    # ── Sheet 3: Settings classes ──
    ws3 = wb.create_sheet("Settings Classes")
    s_headers = [
        "Class",
        "Save/Load",
        "operator==",
        "equals() method",
        "C++ Test",
        "Py Binding: operator==",
        "Py Binding: pickle",
        "Py Binding: file I/O (save/load)",
        "Py Test",
        "Header File",
    ]
    for c, h in enumerate(s_headers, 1):
        cell = ws3.cell(row=1, column=c, value=h)
        cell.font = hdr_font
        cell.fill = hdr_fill

    row = 2
    for cls_name in sorted(class_infos.keys()):
        info = class_infos[cls_name]
        if not info.is_settings:
            continue
        vals = [
            cls_name,
            _yesno(info.has_save_load),
            _yesno(info.has_operator_eq),
            _yesno(info.has_equals_method),
            _yesno(bool(info.cpp_test_files)),
            _yesno(info.py_has_equals),
            _yesno(info.py_has_pickle),
            _yesno(info.py_has_save_to),
            _yesno(bool(info.py_test_files)),
            info.header_file,
        ]
        for c, v in enumerate(vals, 1):
            cell = ws3.cell(row=row, column=c, value=v)
            if c in (2, 3, 4, 5, 6, 7, 8, 9):
                cell.alignment = center
                if v == "✓":
                    cell.fill = yes_fill
                elif v == "✗":
                    cell.fill = no_fill
        row += 1

    for i, w in enumerate([45, 12, 14, 16, 12, 24, 22, 30, 12, 60], 1):
        ws3.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w
    ws3.freeze_panes = "A2"
    ws3.auto_filter.ref = f"A1:{openpyxl.utils.get_column_letter(len(s_headers))}{row - 1}"

    wb.save(output_path)
    return output_path


def write_csv(class_infos: dict[str, ClassInfo], output_path: Path) -> Path:
    import csv

    rows = []
    for cls_name in sorted(class_infos.keys()):
        info = class_infos[cls_name]
        rows.append(
            {
                "Class": cls_name,
                "Is Settings": "yes" if info.is_settings else "no",
                "Save/Load": "yes" if info.has_save_load else "no",
                "operator==": "yes" if info.has_operator_eq else "no",
                "equals() method": "yes" if info.has_equals_method else "no",
                "C++ Roundtrip Test": "yes" if info.cpp_test_files else "no",
                "Python Exposure File": info.py_expose_file,
                "Py Binding: operator==": "yes" if info.py_has_equals else "no",
                "Py Binding: pickle": "yes" if info.py_has_pickle else "no",
                "Py Binding: file I/O (save/load)": "yes" if info.py_has_save_to else "no",
                "Python Roundtrip Test": "yes" if info.py_test_files else "no",
                "Header File": info.header_file,
                "C++ Test Files": "; ".join(info.cpp_test_files),
                "Python Test Files": "; ".join(info.py_test_files),
            }
        )

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        if not rows:
            f.write("No classes found\n")
            return output_path
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    return output_path


# ── main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Track serialization compliance across the Tudat codebase."
    )
    parser.add_argument("--csv", action="store_true", help="Output CSV instead of Excel")
    parser.add_argument(
        "--output", "-o", type=Path, default=None, help=f"Output path (default: {DEFAULT_EXCEL})"
    )
    args = parser.parse_args()

    print("Scanning headers for serializable classes...")
    class_infos = scan_headers()
    print(f"  Found {len(class_infos)} classes")

    print("Scanning C++ test files...")
    scan_cpp_tests(class_infos)

    print("Scanning Python exposure files...")
    scan_python_exposure(class_infos)

    print("Scanning Python test files...")
    scan_python_tests(class_infos)

    # Output
    if args.csv:
        path = write_csv(class_infos, args.output or DEFAULT_EXCEL.with_suffix(".csv"))
    else:
        path = write_excel(class_infos, args.output or DEFAULT_EXCEL)

    print(f"\nDone. Report written to: {path}")

    # Quick summary
    total = len(class_infos)
    sl = sum(1 for i in class_infos.values() if i.has_save_load)
    op_eq = sum(1 for i in class_infos.values() if i.has_operator_eq)
    eq_method = sum(1 for i in class_infos.values() if i.has_equals_method)
    cpp = sum(1 for i in class_infos.values() if i.cpp_test_files)
    py_eq = sum(1 for i in class_infos.values() if i.py_has_equals)
    py_pk = sum(1 for i in class_infos.values() if i.py_has_pickle)
    py_st = sum(1 for i in class_infos.values() if i.py_has_save_to)
    py_t = sum(1 for i in class_infos.values() if i.py_test_files)
    settings = sum(1 for i in class_infos.values() if i.is_settings)
    print(f"\n{'='*60}")
    print(f"  Total: {total}  |  Settings: {settings}")
    print(
        f"  save/load: {sl}  |  operator==: {op_eq}  |  equals(): {eq_method}  |  C++ tests: {cpp}"
    )
    print(
        f"  Py Binding: operator==: {py_eq}  |  Py Binding: pickle: {py_pk}  |  Py Binding: file I/O: {py_st}"
    )
    print(f"  Py tests:  {py_t}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
