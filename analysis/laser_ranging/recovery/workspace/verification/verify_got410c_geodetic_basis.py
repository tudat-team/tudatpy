"""Manual source-formula controls at frozen epochs, without running range cases."""

from pathlib import Path
import hashlib
import importlib.util
import json
import math

import numpy as np
import tudatpy.kernel
from tudatpy.astro import time_representation as time

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
HELPER = Path(
    "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-signature/analysis/laser_ranging/ocean_loading.py"
)
spec = importlib.util.spec_from_file_location("geodetic_loading_helper", HELPER)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
product = WORK / "research/ocean-loading/toc_got410c_harmod.hps"
rows = [
    line.split()
    for line in product.read_text().splitlines()
    if line and line[0] in "HSD" and line[1:2] == " "
]
harmonics = {row[1]: [float(x.replace("D", "E")) for x in row[2:]] for row in rows if row[0] == "H"}
xyz = {row[1]: np.array([float(x) for x in row[2:5]]) for row in rows if row[0] == "S"}
historical_path = WORK / "logs/got410c-helper-simplified-manual-check.json"
historical_bytes = historical_path.read_bytes()
historical = json.loads(historical_bytes)


class IdentityConverter:
    """Expose the callback's TT harmonic calculation at exactly selected epochs."""

    def convert_time(self, input_scale, output_scale, value, position):
        """Assert TDB-to-TT geocentre routing and return the prescribed TT value."""
        assert input_scale == time.tdb_scale and output_scale == time.tt_scale
        np.testing.assert_array_equal(position, np.zeros(3))
        return value


results = []
for old in historical["sites"]:
    site = old["site"]
    # Independent closed-form REF_ELL latitude from the pinned source, not its execution.
    x, y, z = xyz[site]
    horizontal = math.hypot(x, y)
    radius = np.linalg.norm(xyz[site])
    a, f = 6378137.0, 1 / 298.257
    e2 = 2 * f - f * f
    mu = math.atan(z / horizontal * ((1 - f) + e2 * a / radius))
    phi = math.atan(
        ((1 - f) * z + e2 * a * math.sin(mu) ** 3)
        / ((1 - f) * (horizontal - e2 * a * math.cos(mu) ** 3))
    )
    lam = math.atan2(y, x)
    cp, sp, cl, sl = math.cos(phi), math.sin(phi), math.cos(lam), math.sin(lam)
    hand_basis = np.array([[cp * cl, -sl, -sp * cl], [cp * sl, cl, -sp * sl], [sp, 0, cp]])
    callback, meta = helper.got410c_ocean_loading(product, site, xyz[site], IdentityConverter())
    basis = np.array(meta["uen_to_itrs"])
    np.testing.assert_allclose(basis.T @ basis, np.eye(3), rtol=0, atol=1e-15)
    np.testing.assert_allclose(np.cross(basis[:, 1], basis[:, 2]), basis[:, 0], rtol=0, atol=1e-15)
    np.testing.assert_allclose(basis, hand_basis, rtol=0, atol=2e-12)
    up = xyz[site] / radius
    east = np.cross([0, 0, 1], up)
    east /= np.linalg.norm(east)
    radial_basis = np.column_stack((up, east, np.cross(up, east)))
    samples = []
    for saved in old["samples"]:
        tt = saved["tt_seconds"]
        uen = np.zeros(3)
        for row in rows:
            if row[0] != "D" or row[2] != site:
                continue
            phase, omega, accel = harmonics[row[1]]
            argument = phase + omega * tt + 0.5 * accel * tt**2
            coeff = np.array([float(x) for x in row[3:]])
            uen += coeff[:3] * math.cos(argument) + coeff[3:] * math.sin(argument)
        u, e, n = uen
        hand = np.array(
            [cl * (cp * u - sp * n) - sl * e, sl * (cp * u - sp * n) + cl * e, sp * u + cp * n]
        )
        actual = callback(tt)
        radial = radial_basis @ uen
        np.testing.assert_allclose(radial, saved["itrs_displacement_m"], rtol=0, atol=1e-17)
        np.testing.assert_allclose(actual, hand, rtol=0, atol=5e-14)
        assert np.isfinite(actual).all()
        samples.append(
            {
                "tt_seconds": tt,
                "uen_harmonic_sum_m": uen.tolist(),
                "geodetic_xyz_m": actual.tolist(),
                "historical_radial_xyz_m": radial.tolist(),
                "geodetic_minus_radial_xyz_m": (actual - radial).tolist(),
                "geodetic_minus_radial_norm_m": float(np.linalg.norm(actual - radial)),
                "independent_hand_rotation_error_m": float(np.linalg.norm(actual - hand)),
            }
        )
    results.append(
        {
            "site": site,
            "provider_xyz_m": xyz[site].tolist(),
            "basis_ellipsoid": meta["basis_ellipsoid"],
            "basis_determinant": float(np.linalg.det(basis)),
            "orthonormal_max_abs_error": float(np.max(abs(basis.T @ basis - np.eye(3)))),
            "basis_vs_source_ref_ell_max_abs_error": float(np.max(abs(basis - hand_basis))),
            "samples": samples,
        }
    )
assert historical_path.read_bytes() == historical_bytes
result = {
    "helper_sha256": hashlib.sha256(HELPER.read_bytes()).hexdigest(),
    "kernel_sha256": hashlib.sha256(Path(tudatpy.kernel.__file__).read_bytes()).hexdigest(),
    "product_sha256": hashlib.sha256(product.read_bytes()).hexdigest(),
    "historical_radial_control": str(historical_path),
    "historical_radial_control_sha256": hashlib.sha256(historical_bytes).hexdigest(),
    "historical_radial_control_unchanged": True,
    "provider_basis_source_revision": helper.BASIS_SOURCE_REVISION,
    "scope": "Three selected sites, frozen TT epochs0/765201600; no range solves or residual-based choice.",
    "sites": results,
}
out = WORK / "logs/got410c-geodetic-basis-verification.json"
out.write_text(json.dumps(result, indent=2) + "\n")
for row in results:
    print(
        row["site"],
        "basis_error",
        row["basis_vs_source_ref_ell_max_abs_error"],
        "delta_um",
        [1e6 * s["geodetic_minus_radial_norm_m"] for s in row["samples"]],
    )
