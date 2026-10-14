#!/usr/bin/env python3
"""Execute the unchanged pinned IGN routine; no network or reference-source copy."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_cache", type=Path)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    source = args.reference_cache / "parametric.f"
    digest = "e3817d992f3120ff420563151a87b7705377b3ffaebcbc6e9eae8615cf4dbda3"
    if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
        raise ValueError("Original IGN routine checksum mismatch")
    driver = """program reference
  implicit none
  integer :: model, i, k
  real(8) :: dt, a1, t1, a2, t2, d
  real(8), parameter :: elapsed(7) = [0d0,1d-14,1d-7,0.01d0,1d0,20d0,1d5]
  do k=1,3
    a1 = (-1d0)**k * (12d0 + k)
    a2 = (-1d0)**(k+1) * (31d0 - k)
    t1 = 0.003d0 * k
    t2 = 0.7d0 * k
    do model=0,5
      do i=1,7
        dt=elapsed(i)
        call parametric(model,dt,a1,t1,a2,t2,d)
        write(*,'(I1,7(",",ES25.17E3))') model,dt,a1,t1,a2,t2,d
      end do
    end do
  end do
end program reference
"""
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory)
        (temporary / "driver.f90").write_text(driver)
        executable = temporary / "reference"
        command = [
            args.compiler,
            "-O0",
            str(source),
            str(temporary / "driver.f90"),
            "-o",
            str(executable),
        ]
        subprocess.run(command, check=True)
        values = subprocess.check_output([str(executable)], text=True)
    data = "model,elapsed_julian_years,a1_mm,t1_years,a2_mm,t2_years,correction_mm\n" + values
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "ign_psd_reference.csv").write_text(data)
    metadata = {
        "reference": json.loads((args.reference_cache / "parametric.f.source.json").read_text()),
        "compiler": subprocess.check_output([args.compiler, "--version"], text=True).splitlines()[
            0
        ],
        "flags": ["-O0"],
        "driver_sha256": hashlib.sha256(driver.encode()).hexdigest(),
        "output_sha256": hashlib.sha256(data.encode()).hexdigest(),
        "rows": len(values.splitlines()),
        "scope": "Unmodified IGN equation oracle, selectors 0--5; no complete trajectory reference claimed",
    }
    (args.output / "ign_psd_reference.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
