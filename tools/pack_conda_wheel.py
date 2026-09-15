#!/usr/bin/env python3
"""Repack an installed conda tudatpy into a pip/uv wheel.

Run inside the conda/pixi/micromamba env that has tudatpy. Does not compile
anything: it copies kernel.so/.pyd and vendors the native libs it links.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import os
import platform
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import zipfile
from pathlib import Path

SKIP_SUFFIXES = {".cpp", ".cc", ".cxx", ".h", ".hpp", ".c", ".cmake", ".so", ".pyd", ".dylib"}
SKIP_NAMES = {"CMakeLists.txt"}
SKIP_DIR_NAMES = {"__pycache__", "temp"}

LINUX_SKIP_SONAMES = {
    "linux-vdso.so.1",
    "ld-linux-x86-64.so.2",
    "ld-linux-aarch64.so.1",
    "libc.so.6",
    "libm.so.6",
    "libdl.so.2",
    "librt.so.1",
    "libpthread.so.0",
    "libresolv.so.2",
}


def conda_prefix() -> Path:
    prefix = os.environ.get("CONDA_PREFIX") or sys.prefix
    return Path(prefix).resolve()


def site_packages() -> Path:
    return Path(sysconfig.get_path("purelib")).resolve()


def copy_python_package(src: Path, dst: Path) -> None:
    for path in src.rglob("*"):
        rel = path.relative_to(src)
        if any(part in SKIP_DIR_NAMES for part in rel.parts):
            continue
        if path.is_dir():
            (dst / rel).mkdir(exist_ok=True)
            continue
        if path.suffix.lower() in SKIP_SUFFIXES or path.name in SKIP_NAMES:
            continue
        (dst / rel.parent).mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dst / rel)


def kernel_path(pkg: Path) -> Path:
    for name in ("kernel.so", "kernel.pyd", "kernel.abi3.so"):
        candidate = pkg / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No kernel extension in {pkg}")


def wheel_tag() -> str:
    impl = f"cp{sys.version_info.major}{sys.version_info.minor}"
    machine = platform.machine().lower()
    if sys.platform == "darwin":
        arch = "arm64" if machine in {"arm64", "aarch64"} else "x86_64"
        return f"{impl}-{impl}-macosx_13_0_{arch}"
    if sys.platform.startswith("linux"):
        arch = "aarch64" if machine in {"arm64", "aarch64"} else "x86_64"
        return f"{impl}-{impl}-linux_{arch}"
    if sys.platform == "win32":
        arch = "win_arm64" if machine in {"arm64", "aarch64"} else "win_amd64"
        return f"{impl}-{impl}-{arch}"
    raise SystemExit(f"Unsupported platform: {sys.platform} {machine}")


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kwargs)


def macos_deps(binary: Path) -> list[str]:
    out = run(["otool", "-L", str(binary)]).stdout
    deps = []
    for line in out.splitlines()[1:]:
        name = line.strip().split(" ", 1)[0]
        if name.startswith("/usr/lib/") or name.startswith("/System/"):
            continue
        deps.append(name)
    return deps


def linux_deps(binary: Path) -> list[tuple[str, str]]:
    out = run(["ldd", str(binary)]).stdout
    found = []
    for line in out.splitlines():
        line = line.strip()
        if " => " in line:
            soname, rest = line.split(" => ", 1)
            soname = soname.strip()
            path = rest.split(" ", 1)[0].strip()
            if path in {"not", "not found"} or rest.startswith("not found"):
                if soname not in LINUX_SKIP_SONAMES:
                    raise RuntimeError(f"{binary} has unresolved dependency {soname}")
                continue
            if soname in LINUX_SKIP_SONAMES:
                continue
            found.append((soname, path))
        elif line.startswith("/") and "ld-linux" not in line:
            found.append((Path(line.split()[0]).name, line.split()[0]))
    return found


def resolve_in_prefix(name: str, prefix: Path) -> Path | None:
    base = Path(name.split("/")[-1]).name
    candidates = [
        prefix / "lib" / Path(name).name,
        prefix / "lib" / name,
        prefix / "lib" / base,
        prefix / "Library" / "bin" / Path(name).name,
        prefix / "Library" / "lib" / Path(name).name,
        prefix / "Library" / "bin" / base,
        prefix / "DLLs" / Path(name).name,
        prefix / "bin" / Path(name).name,
    ]
    for path in candidates:
        if path.exists():
            return path
    for folder in (prefix / "lib", prefix / "Library" / "bin"):
        if folder.is_dir():
            matches = list(folder.glob(base)) + list(folder.glob(base + "*"))
            if matches:
                return matches[0]
    return None


def vendor_macos(kernel: Path, prefix: Path, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    queue = [kernel]
    seen: set[str] = set()
    while queue:
        current = queue.pop()
        for dep in macos_deps(current):
            base = Path(dep).name
            if base == "libc++.1.dylib":
                run(
                    [
                        "install_name_tool",
                        "-change",
                        dep,
                        "/usr/lib/libc++.1.dylib",
                        str(current),
                    ]
                )
                continue
            if base in seen:
                rel = (
                    f"@loader_path/.dylibs/{base}" if current == kernel else f"@loader_path/{base}"
                )
                run(["install_name_tool", "-change", dep, rel, str(current)])
                continue
            src = resolve_in_prefix(dep, prefix)
            if src is None:
                raise RuntimeError(
                    f"{current.name} links {dep} which is not in {prefix} and was not vendored"
                )
            seen.add(base)
            shutil.copy2(src, dest_dir / base)
            copied = dest_dir / base
            rel = f"@loader_path/.dylibs/{base}" if current == kernel else f"@loader_path/{base}"
            run(["install_name_tool", "-change", dep, rel, str(current)])
            run(["install_name_tool", "-id", f"@loader_path/.dylibs/{base}", str(copied)])
            queue.append(copied)
    for path in [kernel, *dest_dir.glob("*")]:
        run(["codesign", "--force", "--sign", "-", str(path)])


def vendor_linux(kernel: Path, prefix: Path, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    queue = [kernel]
    seen: set[str] = set()
    while queue:
        current = queue.pop()
        for soname, path in linux_deps(current):
            if soname in LINUX_SKIP_SONAMES or soname in seen:
                continue
            src = (
                Path(path)
                if path.startswith("/") and Path(path).exists()
                else resolve_in_prefix(soname, prefix)
            )
            if src is None:
                raise RuntimeError(f"{current.name} links {soname} which could not be resolved")
            resolved = src.resolve()
            if resolved.is_relative_to(Path("/usr")) and not str(resolved).startswith(str(prefix)):
                continue
            if not str(resolved).startswith(str(prefix)):
                if soname.startswith(("libstdc++", "libgcc_s", "libgomp")):
                    continue
                raise RuntimeError(
                    f"{current.name} links {soname} -> {resolved}, outside conda prefix {prefix}"
                )
            seen.add(soname)
            shutil.copy2(src, dest_dir / src.name)
            queue.append(dest_dir / src.name)
    run(["patchelf", "--set-rpath", "$ORIGIN/.libs", str(kernel)])
    for path in dest_dir.glob("*.so*"):
        run(["patchelf", "--set-rpath", "$ORIGIN", str(path)])


def windows_imports(dll: Path) -> list[str]:
    try:
        import pefile  # type: ignore
    except ImportError as exc:
        raise RuntimeError("pefile is required to vendor Windows DLLs") from exc
    # Parse from bytes so Windows does not keep a file lock on the copy
    # (TemporaryDirectory cleanup otherwise hits WinError 32).
    pe = pefile.PE(data=dll.read_bytes())
    names = []
    try:
        if not hasattr(pe, "DIRECTORY_ENTRY_IMPORT"):
            return names
        for entry in pe.DIRECTORY_ENTRY_IMPORT:
            names.append(entry.dll.decode("ascii", "ignore"))
    finally:
        pe.close()
    return names


def vendor_windows(kernel: Path, prefix: Path, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    search_dirs = [
        prefix / "Library" / "bin",
        prefix / "Library" / "lib",
        prefix / "DLLs",
        prefix / "bin",
        prefix / "Library" / "mingw-w64" / "bin",
    ]
    queue = [kernel]
    seen: set[str] = set()
    systemish = {
        "kernel32.dll",
        "user32.dll",
        "advapi32.dll",
        "shell32.dll",
        "ole32.dll",
        "oleaut32.dll",
        "ws2_32.dll",
        "ntdll.dll",
        "msvcrt.dll",
        "vcruntime140.dll",
        "vcruntime140_1.dll",
        "msvcp140.dll",
        "python3.dll",
        f"python{sys.version_info.major}{sys.version_info.minor}.dll",
    }
    while queue:
        current = queue.pop()
        for name in windows_imports(current):
            key = name.lower()
            if key in seen or key in systemish or key.startswith(("api-ms-win-", "ext-ms-")):
                continue
            src = None
            for folder in search_dirs:
                if not folder.is_dir():
                    continue
                direct = folder / name
                if direct.exists():
                    src = direct
                    break
                for candidate in folder.iterdir():
                    if candidate.name.lower() == key:
                        src = candidate
                        break
                if src is not None:
                    break
            if src is None:
                if any(token in key for token in ("boost", "cspice", "nrlmsise", "sofa")):
                    raise RuntimeError(
                        f"{current.name} imports {name} which was not found under {prefix}"
                    )
                continue
            seen.add(key)
            shutil.copy2(src, dest_dir / src.name)
            queue.append(dest_dir / src.name)


def write_dist_info(stage: Path, version: str, tag: str) -> None:
    dist = stage / f"tudatpy-{version}.dist-info"
    dist.mkdir()
    (dist / "METADATA").write_text(
        "\n".join(
            [
                "Metadata-Version: 2.1",
                "Name: tudatpy",
                f"Version: {version}",
                "Summary: TU Delft Astrodynamics Toolbox (repacked official conda build)",
                f"Requires-Python: =={sys.version_info.major}.{sys.version_info.minor}.*",
                "Requires-Dist: numpy>=1.26.4,<2",
                "Requires-Dist: scipy",
                "Requires-Dist: pandas",
                "Requires-Dist: matplotlib",
                "Requires-Dist: astropy",
                "Requires-Dist: astropy-healpix",
                "Requires-Dist: astroquery>=0.4.8",
                "Requires-Dist: spiceypy",
                "Requires-Dist: tabulate",
                "Requires-Dist: colorama",
                "Requires-Dist: tqdm",
                "Requires-Dist: requests",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (dist / "WHEEL").write_text(
        "\n".join(
            [
                "Wheel-Version: 1.0",
                "Generator: pack_conda_wheel.py",
                "Root-Is-Purelib: false",
                f"Tag: {tag}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (dist / "top_level.txt").write_text("tudatpy\n", encoding="utf-8")
    license_src = Path(__file__).resolve().parents[1] / "LICENSE"
    if license_src.is_file():
        shutil.copy2(license_src, dist / "LICENSE")


def write_record(stage: Path, version: str) -> None:
    dist = stage / f"tudatpy-{version}.dist-info"
    lines = []
    for path in sorted(stage.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(stage).as_posix()
        data = path.read_bytes()
        digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
        lines.append(f"{rel},sha256={digest},{len(data)}")
    lines.append(f"tudatpy-{version}.dist-info/RECORD,,")
    (dist / "RECORD").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_index(directory: Path) -> None:
    files = sorted(
        p.name for p in directory.iterdir() if p.suffix == ".whl" or p.name.endswith(".tar.gz")
    )
    links = "\n".join(f'    <a href="{name}">{name}</a><br/>' for name in files)
    (directory / "index.html").write_text(
        "<!DOCTYPE html>\n<html><head><title>tudatpy wheels</title></head>\n"
        f"<body>\n<h1>tudatpy</h1>\n{links}\n</body></html>\n",
        encoding="utf-8",
    )


def pack(output_dir: Path) -> Path:
    import tudatpy

    version = getattr(tudatpy, "__version__", "0.0.0")
    src = site_packages() / "tudatpy"
    if not src.is_dir():
        raise SystemExit(f"tudatpy is not installed in {site_packages()}")
    prefix = conda_prefix()
    tag = wheel_tag()
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="tudatpy-wheel-", ignore_cleanup_errors=True) as raw:
        stage = Path(raw)
        pkg = stage / "tudatpy"
        pkg.mkdir()
        copy_python_package(src, pkg)
        src_kernel = kernel_path(src)
        dst_kernel = pkg / src_kernel.name
        shutil.copy2(src_kernel, dst_kernel)

        if sys.platform == "darwin":
            vendor_macos(dst_kernel, prefix, pkg / ".dylibs")
        elif sys.platform.startswith("linux"):
            vendor_linux(dst_kernel, prefix, pkg / ".libs")
        elif sys.platform == "win32":
            vendor_windows(dst_kernel, prefix, pkg)

        write_dist_info(stage, version, tag)
        write_record(stage, version)
        wheel_path = output_dir / f"tudatpy-{version}-{tag}.whl"
        if wheel_path.exists():
            wheel_path.unlink()
        with zipfile.ZipFile(wheel_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for path in stage.rglob("*"):
                if path.is_file():
                    zf.write(path, path.relative_to(stage).as_posix())
        return wheel_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("wheels"))
    parser.add_argument("--index", action="store_true", help="Write index.html in the output dir")
    args = parser.parse_args()
    wheel = pack(args.output_dir)
    print(f"wrote {wheel} ({wheel.stat().st_size} bytes)")
    if args.index:
        write_index(args.output_dir)
        print(f"wrote {args.output_dir / 'index.html'}")


if __name__ == "__main__":
    main()
