# Pip wheels from the official conda builds

Conda remains the supported way to install TudatPy. This repository can also publish **pip wheels** so students who already use `venv` + `pip` (or [uv](https://docs.astral.sh/uv/), a faster drop-in for pip) do not need a conda distribution.

Nothing is compiled here. GitHub Actions installs the existing `tudat-team::tudatpy` conda package and **repacks** `kernel.so` / `kernel.pyd` plus the native libraries it links (Boost, CSPICE, NRLMSISE-00) into a standard `.whl`.

## Install a wheel

Wheels are attached to GitHub Releases (Linux x86_64, macOS arm64, macOS Intel, Windows amd64 × CPython 3.10/3.11/3.12). After a release tagged `v1.0.0`:

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install tudatpy \
  --find-links https://github.com/tudat-team/tudatpy/releases/download/v1.0.0/index.html
```

Optional, same command with uv:

```bash
uv pip install tudatpy \
  --find-links https://github.com/tudat-team/tudatpy/releases/download/v1.0.0/index.html
```

SPICE kernels and other data files are **not** inside the wheel. Unpack [tudat-resources v2.4](https://github.com/tudat-team/tudat-resources/releases/download/v2.4/resource.tar.gz) into `~/.tudat` so that `~/.tudat/resource/quadrature/gaussianNodes.txt` exists.

## Produce wheels

The **Wheels** workflow does **not** run on push or pull request. It runs on:

- **`release` published** — builds wheels and uploads them (plus `index.html`) to that release
- **`workflow_dispatch`** — builds wheels as Actions artifacts only (nothing attached to a release)

### First 1.0.0 wheel release

CI installs conda `tudatpy=1.0.0` (pinned in the workflow). Merging this PR does not publish wheels.

1. Merge the workflow to the default branch (so `release` / `workflow_dispatch` are available).
2. GitHub → **Releases → Draft a new release**.
3. Create tag **`v1.0.0`** on that commit (or any commit that contains `.github/workflows/wheels.yml`).
4. Publish the release.

That is the first trigger that both builds and attaches wheels. `workflow_dispatch` is only for a dry run.

## What this does not replace

- The conda `environment.yaml` path in the docs.
- Building Tudat from this source tree (`install.py` / CMake).
- Optional conda-only extras (FFTW, PaGMO, Mars Climate Database).
