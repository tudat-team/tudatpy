import subprocess
import re
import csv
import argparse
from pathlib import Path
from collections import defaultdict

EXTENSIONS = (".h", ".hpp", ".hh", ".cpp", ".cc", ".cxx")
TEST_PATTERNS = re.compile(r"(test|spec|fixture)", re.IGNORECASE)
CLASS_RE = re.compile(r"\b(class|struct)\s+([A-Za-z_]\w*)")
SERIALIZE_RE = re.compile(r"\b(?:serialize|save|load)\s*\(")
CEREAL_ARCHIVE_RE = re.compile(r"cereal::(Binary|JSON|XML|Portable)(Input|Output)Archive")
PICKLE_RE = re.compile(r"(?:__getstate__|__setstate__|pickle|enable_pickling)", re.IGNORECASE)

DEFAULT_CSV_OUTPUT = Path(__file__).with_name("find_tested_serialized.csv")

def get_files():
    p = subprocess.run(["rg", "--files"], stdout=subprocess.PIPE, text=True)
    return [f for f in p.stdout.splitlines() if f.endswith(EXTENSIONS)]

def is_test_file(path: str) -> bool:
    return bool(TEST_PATTERNS.search(path))

def _extract_serializable_classes(file):
    try:
        lines = Path(file).read_text(errors="ignore").splitlines()
    except:
        return set()

    results = set()
    current_class = None
    class_brace_depth = None
    pending_class = None
    pending_depth = None
    brace_depth = 0

    for line in lines:
        m = CLASS_RE.search(line)
        if m:
            pending_class = m.group(2)
            pending_depth = brace_depth

        brace_depth += line.count("{")
        brace_depth -= line.count("}")

        if pending_class is not None and brace_depth > pending_depth:
            current_class = pending_class
            class_brace_depth = brace_depth
            pending_class = None
            pending_depth = None

        if current_class is not None and brace_depth < class_brace_depth:
            current_class = None
            class_brace_depth = None

        if SERIALIZE_RE.search(line):
            if current_class and current_class != "Archive":
                results.add(current_class)

    return results

def find_classes_with_serialize(files):
    """Returns dict of {class_name: [source_files]}"""
    found = defaultdict(list)
    for file in files:
        if is_test_file(file):
            continue
        classes = _extract_serializable_classes(file)
        for cls in classes:
            found[cls].append(file)
    return found

def _makes_archive_call_with(content, cls):
    archive_vars = re.findall(
        r"cereal::(?:Binary|JSON|XML|Portable)(?:Input|Output)Archive\s+(\w+)",
        content
    )
    if not archive_vars:
        return False
    cls_re = re.compile(rf"\b{re.escape(cls)}\b")
    for var in archive_vars:
        archive_call_re = re.compile(rf"\b{re.escape(var)}\s*\(")
        call_positions = [m.start() for m in archive_call_re.finditer(content)]
        cls_positions  = [m.start() for m in cls_re.finditer(content)]
        for cp in call_positions:
            for clp in cls_positions:
                if abs(cp - clp) < 300:
                    return True
    return False

def find_serialize_tests(test_files, serializable_classes):
    """Returns dict of {class_name: [test_files_that_reference_it]}"""
    test_file_contents = {}
    for f in test_files:
        try:
            text = Path(f).read_text(errors="ignore")
            if CEREAL_ARCHIVE_RE.search(text):
                test_file_contents[f] = text
        except:
            pass

    coverage = defaultdict(list)
    for cls in serializable_classes:
        for tf, content in test_file_contents.items():
            if _makes_archive_call_with(content, cls):
                coverage[cls].append(tf)

    return coverage

def find_pickle_support(serializable_classes):
    """Returns dict of {class_name: [expose_files_with_pickle_support]}"""
    try:
        result = subprocess.run(
            ["find", ".", "-name", "expose*.cpp", "-o", "-name", "expose*.hpp"],
            stdout=subprocess.PIPE,
            text=True,
            cwd="."
        )
        expose_files = [f for f in result.stdout.splitlines() if f]
    except:
        expose_files = []

    pickle_coverage = defaultdict(list)
    for cls in serializable_classes:
        for expose_file in expose_files:
            try:
                content = Path(expose_file).read_text(errors="ignore")
                # Look for class name near py::pickle patterns
                cls_re = re.compile(rf"\b{re.escape(cls)}\b")
                pickle_re = re.compile(r"py::pickle", re.IGNORECASE)
                
                if cls_re.search(content) and pickle_re.search(content):
                    pickle_coverage[cls].append(expose_file)
            except:
                pass

    return pickle_coverage

def write_csv_report(output_path, serializable, coverage, pickle_coverage):
    rows = []
    for cls in sorted(serializable):
        rows.append(
            {
                "class_name": cls,
                "tested": "yes" if coverage[cls] else "no",
                "pickle_support": "yes" if pickle_coverage[cls] else "no",
                "source_files": "; ".join(serializable[cls]),
                "test_files": "; ".join(sorted(coverage[cls])),
            }
        )

    with Path(output_path).open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=["class_name", "tested", "pickle_support", "source_files", "test_files"],
        )
        writer.writeheader()
        writer.writerows(rows)

    return output_path

def main():
    parser = argparse.ArgumentParser(description="Find serializable classes and matching tests.")
    parser.add_argument(
        "--csv",
        nargs="?",
        const=str(DEFAULT_CSV_OUTPUT),
        default=None,
        help=f"Write a CSV report for Excel. Defaults to {DEFAULT_CSV_OUTPUT} when used without a path.",
    )
    args = parser.parse_args()

    all_files = get_files()
    source_files = [f for f in all_files if not is_test_file(f)]
    test_files = [f for f in all_files if is_test_file(f)]

    serializable = find_classes_with_serialize(source_files)
    coverage = find_serialize_tests(test_files, serializable)
    pickle_coverage = find_pickle_support(serializable)

    tested = {cls for cls in serializable if coverage[cls]}
    untested = {cls for cls in serializable if not coverage[cls]}
    pickle_enabled = {cls for cls in serializable if pickle_coverage[cls]}
    pickle_disabled = {cls for cls in serializable if not pickle_coverage[cls]}

    print(f"\n{'='*50}")
    print(f"  {len(tested)} tested / {len(untested)} untested / {len(serializable)} total")
    print(f"  {len(pickle_enabled)} with pickle / {len(pickle_disabled)} without pickle")
    print(f"{'='*50}\n")

    if untested:
        print("❌ NO SERIALIZE TEST FOUND:")
        for cls in sorted(untested):
            sources = ", ".join(serializable[cls])
            pickle_status = "✓ pickle" if cls in pickle_enabled else "✗ no pickle"
            print(f"   {cls:<40} {pickle_status:<15} (defined in: {sources})")

    if tested:
        print("\n✅ SERIALIZE TEST EXISTS:")
        for cls in sorted(tested):
            pickle_status = "✓ pickle" if cls in pickle_enabled else "✗ no pickle"
            for tf in sorted(coverage[cls]):
                print(f"   {cls:<40} {pickle_status:<15} → {tf}")

    if args.csv is not None:
        output_path = write_csv_report(args.csv, serializable, coverage, pickle_coverage)
        print(f"\nCSV written to: {output_path}")

if __name__ == "__main__":
    main()