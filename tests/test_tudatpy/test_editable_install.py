"""Refresh editable installations without compiling or replacing the kernel."""

import importlib.util
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


def load_script(name):
    path = Path(__file__).resolve().parents[2] / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"tudatpy_{name}_script", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def installation(monkeypatch, tmp_path):
    source = tmp_path / "checkout/src/tudatpy"
    source.mkdir(parents=True)
    (source / "__init__.py").write_text("")
    prefix = tmp_path / "environment"
    pylib = prefix / "site-packages"
    pylib.mkdir(parents=True)
    build = tmp_path / "build"
    kernel_dir = build / "src/tudatpy"
    kernel_dir.mkdir(parents=True)
    kernel = kernel_dir / ("kernel.pyd" if sys.platform == "win32" else "kernel.so")
    kernel.write_bytes(b"existing kernel")
    (build / "tudatpy-stubs").mkdir()

    probe = tmp_path / "symlink-probe"
    try:
        probe.symlink_to(source / "__init__.py")
    except OSError as error:
        if sys.platform == "win32" and error.winerror == 1314:
            pytest.skip("Creating editable installs requires symlink privileges")
        raise
    probe.unlink()

    installer = load_script("install")
    monkeypatch.setattr(installer, "__file__", str(tmp_path / "checkout/install.py"))
    monkeypatch.setenv("CONDA_PREFIX", str(prefix))
    monkeypatch.setattr(installer.sysconfig, "get_path", lambda name: str(pylib))
    monkeypatch.setattr(
        sys,
        "argv",
        ["install.py", "-e", "--skip-tudat", "--build-dir", str(build)],
    )
    return SimpleNamespace(
        installer=installer,
        source=source,
        package=pylib / "tudatpy",
        pylib=pylib,
        kernel=kernel,
        build=build,
        manifest=build / "manifests/environment.txt",
    )


def test_refresh_makes_new_modules_and_packages_importable(installation):
    observations = installation.source / "estimation/observations"
    observations.mkdir(parents=True)
    (observations.parent / "__init__.py").write_text("")
    (observations / "__init__.py").write_text("")
    installation.installer.Installer().install()
    original_manifest = installation.manifest.read_text()

    # Reproduce a checkout gaining a module used by an already-linked initializer.
    (observations / "_query.py").write_text("value = 42\n")
    (observations / "__init__.py").write_text("from ._query import value\n")
    new_package = installation.source / "new_package"
    new_package.mkdir()
    (new_package / "__init__.py").write_text("value = 43\n")
    (installation.source / "py.typed").touch()

    installation.installer.Installer().install()
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-c",
            f"import sys; sys.path.insert(0, {str(installation.pylib)!r}); "
            "from tudatpy.estimation.observations import value; "
            "from tudatpy import new_package; assert (value, new_package.value) == (42, 43)",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (installation.package / "py.typed").is_symlink()
    assert (installation.package / installation.kernel.name).resolve() == installation.kernel
    assert installation.kernel.read_bytes() == b"existing kernel"
    # The original tracked package directory already owns all new files.
    assert installation.manifest.read_text() == original_manifest
    installation.installer.Installer().install()
    assert installation.manifest.read_text() == original_manifest


def test_refresh_preserves_uninstall_ownership_in_an_existing_directory(installation, monkeypatch):
    installation.package.mkdir()
    unrelated = installation.package / "unrelated.txt"
    unrelated.write_text("preserve me")
    installation.installer.Installer().install()

    (installation.source / "new_module.py").write_text("value = 42\n")
    installation.installer.Installer().install()
    assert (installation.package / "new_module.py").is_symlink()

    uninstaller = load_script("uninstall")
    monkeypatch.setattr(sys, "argv", ["uninstall.py", "--build-dir", str(installation.build)])
    uninstaller.Remover().remove()
    assert unrelated.read_text() == "preserve me"
    assert not list(installation.package.glob("*.py"))
    assert not (installation.package / installation.kernel.name).is_symlink()
    assert not (installation.pylib / "tudatpy-stubs").is_symlink()
    assert not installation.manifest.exists()
    assert (installation.source / "new_module.py").exists()


@pytest.mark.parametrize("existing_installation", ["different_checkout", "regular_copy"])
def test_refresh_refuses_to_mix_python_sources(installation, existing_installation, tmp_path):
    installation.installer.Installer().install()
    original_manifest = installation.manifest.read_text()
    init_link = installation.package / "__init__.py"
    init_link.unlink()
    if existing_installation == "different_checkout":
        foreign_init = tmp_path / "other_checkout/__init__.py"
        foreign_init.parent.mkdir()
        foreign_init.write_text("")
        init_link.symlink_to(foreign_init)
    else:
        init_link.write_text("")
    (installation.source / "new_module.py").write_text("value = 42\n")

    with pytest.raises(RuntimeError, match="not editable from this checkout"):
        installation.installer.Installer()
    assert not (installation.package / "new_module.py").exists()
    assert installation.manifest.read_text() == original_manifest
    assert installation.kernel.read_bytes() == b"existing kernel"
