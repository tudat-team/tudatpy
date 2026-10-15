"""Pinned GOT4.10c loading for analysis station callbacks, without BLQ/HARDISP.

HARPOS arguments use TT (TDT) seconds from J2000 and geocentric radial Up,
East, North amplitudes in metres: https://massloading.net/harpos_format.txt.
The service realization (https://massloading.net/toc/) uses STW105/REAR and
the whole-Earth centre of mass. It differs from ASI's 11-constituent BLQ plus
HARDISP and JPL FES2004; using it does not reproduce either provider exactly.
Only the 17 supplied gravity-tide terms are included, without ocean-pole or
equilibrium terms. Do not duplicate this realization with a BLQ ocean-loading
model. The product narrative calls Up ellipsoid-normal, conflicting with the
HARPOS radial definition; this analysis explicitly follows the latter.
"""

from pathlib import Path
import hashlib

import numpy as np

from tudatpy.astro import time_representation as time

PRODUCT_URL = "https://massloading.net/imsl/load_har_list/toc/got410c/toc_got410c_harmod.hps"
PRODUCT_SHA256 = "67e160364809cec482dc6e122df256d279a87478aeebfc67586cd0fe08311150"
CONSTITUENTS = (
    "K1",
    "K1-",
    "K1+",
    "K2",
    "K2+",
    "M2",
    "M2-",
    "M4",
    "N2",
    "N2-",
    "O1",
    "O1-",
    "P1",
    "Q1",
    "Q1-",
    "S1",
    "S2",
)


def got410c_ocean_loading(file_name, site_name, nominal_station_xyz, time_converter):
    """Return ``(displacement_callback, metadata)`` for custom_station_motion.

    ``file_name`` is the cached, unmodified product, verified by its pinned hash.
    ``site_name`` is case sensitive; use Graz, Herstmon or Apollo for the selected
    laser monuments (GRAZ/HERS identify different nearby monuments). The nominal
    station XYZ is in terrestrial metres. Pass the time converter owned by the
    selected Earth rotation's angles calculator, not an independent default.
    The callback accepts TDB seconds from J2000 and returns terrestrial XYZ metres.
    """
    path = Path(file_name)
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != PRODUCT_SHA256:
        raise ValueError("GOT4.10c product does not match the pinned input")

    harmonics, amplitudes, provider_xyz = {}, {}, None
    radius = 3000.0
    for line in raw.decode("ascii").splitlines():
        if line.startswith("H "):
            name = line[3:11].strip()
            harmonics[name] = [
                float(line[a:b].replace("D", "E")) for a, b in ((13, 26), (28, 47), (49, 59))
            ]
        elif line.startswith("S ") and line[3:11].strip() == site_name:
            provider_xyz = np.array([float(line[a:b]) for a, b in ((13, 26), (27, 40), (41, 54))])
        elif line.startswith("D ") and line[13:21].strip() == site_name:
            amplitudes[line[3:11].strip()] = [
                float(line[a:b])
                for a, b in ((24, 32), (33, 41), (42, 50), (53, 61), (62, 70), (71, 79))
            ]
    if (
        set(harmonics) != set(CONSTITUENTS)
        or set(amplitudes) != set(CONSTITUENTS)
        or provider_xyz is None
    ):
        raise ValueError(f"Expected site {site_name!r} and all 17 defined constituent rows")
    xyz = np.asarray(nominal_station_xyz, dtype=float).reshape(3)
    separation = float(np.linalg.norm(xyz - provider_xyz))
    if not np.all(np.isfinite(xyz)) or separation > radius:
        raise ValueError(f"Station lies outside {site_name}'s {radius:g} m applicability radius")

    up = provider_xyz / np.linalg.norm(provider_xyz)
    east = np.array([-up[1], up[0], 0.0])
    east /= np.linalg.norm(east)
    north = np.cross(up, east)
    basis = np.column_stack((up, east, north))
    phase, omega, acceleration = np.array([harmonics[name] for name in CONSTITUENTS]).T
    coefficients = np.array([amplitudes[name] for name in CONSTITUENTS])

    def displacement(tdb_seconds):
        tt = time_converter.convert_time(
            time.tdb_scale, time.tt_scale, float(tdb_seconds), np.zeros(3)
        )
        argument = phase + omega * tt + 0.5 * acceleration * tt**2
        uen = np.cos(argument) @ coefficients[:, :3] + np.sin(argument) @ coefficients[:, 3:]
        return basis @ uen

    metadata = {
        "model": "GOT4.10c",
        "source_url": PRODUCT_URL,
        "sha256": digest,
        "file": str(path.resolve()),
        "product_release": "GOT4.10c; generated 2017-06-17; HARPOS 2005.03.28",
        "site": site_name,
        "provider_xyz_m": provider_xyz.tolist(),
        "nominal_xyz_m": xyz.tolist(),
        "site_separation_m": separation,
        "applicability_radius_m": radius,
        "constituents": list(CONSTITUENTS),
        "argument_time": "TT seconds since J2000; selected Earth rotation converter at geocentre",
        "amplitude_axes": "geocentric radial Up, East, North (HARPOS specification)",
        "axis_wording_discrepancy": "Product narrative says ellipsoid-normal Up; HARPOS says radial Up. This analysis follows HARPOS.",
        "uen_to_itrs": basis.tolist(),
        "origin": "whole-Earth centre of mass",
        "earth_response": "STW105/REAR",
        "scope": "17 supplied gravity-tide terms including nodal sidebands; ocean-pole and equilibrium terms are not included",
        "comparison": "Different realization from ASI 11-constituent BLQ+HARDISP and JPL FES2004; not exact provider replay",
    }
    return displacement, metadata
