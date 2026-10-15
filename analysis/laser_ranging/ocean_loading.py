"""Pinned GOT4.10c loading for analysis station callbacks, without BLQ/HARDISP.

HARPOS arguments use TT (TDT) seconds from J2000 and geocentric radial Up,
East, North amplitudes in metres: https://massloading.net/harpos_format.txt.
The service realization (https://massloading.net/toc/) uses STW105/REAR and
the whole-Earth centre of mass. It differs from ASI's 11-constituent BLQ plus
HARDISP and JPL FES2004; using it does not reproduce either provider exactly.
Do not also attach a BLQ ocean-loading model, ocean-pole or equilibrium add-ons.
"""

from pathlib import Path
import hashlib
import json

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

    ``file_name`` is the cached, unmodified product with its adjacent SOURCES.json.
    ``site_name`` is case sensitive; use Graz, Herstmon or Apollo for the selected
    laser monuments (GRAZ/HERS identify different nearby monuments). The nominal
    station XYZ is in terrestrial metres. Pass the time converter owned by the
    selected Earth rotation's angles calculator, not an independent default.
    The callback accepts TDB seconds from J2000 and returns terrestrial XYZ metres.
    """
    path = Path(file_name)
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    provenance = json.loads(path.with_name("SOURCES.json").read_text())
    if digest != PRODUCT_SHA256 or any(
        (
            provenance.get("sha256") != PRODUCT_SHA256,
            provenance.get("url") != PRODUCT_URL,
            provenance.get("file") != "toc_got410c_harmod.hps",
            provenance.get("bytes") != len(raw),
        )
    ):
        raise ValueError("GOT4.10c product or source provenance does not match the pinned input")

    lines = raw.decode("ascii").splitlines()
    if lines[0] != "HARPOS  Format version of 2005.03.28" or lines[-1] != lines[0]:
        raise ValueError("Missing pinned HARPOS header/trailer")
    harmonics, amplitudes, sites, radii = {}, {}, {}, []
    for line in lines[1:-1]:
        if not line or line.startswith("#"):
            continue
        kind = line[0]
        if kind == "A":
            radii.append(float(line[3:17]))
        elif kind == "H":
            name = line[3:11].strip()
            if name in harmonics:
                raise ValueError(f"Duplicate harmonic {name}")
            harmonics[name] = [
                float(line[a:b].replace("D", "E")) for a, b in ((13, 26), (28, 47), (49, 59))
            ]
        elif kind == "S":
            name = line[3:11].strip()
            if name in sites:
                raise ValueError(f"Duplicate site {name}")
            sites[name] = np.array([float(line[a:b]) for a, b in ((13, 26), (27, 40), (41, 54))])
        elif kind == "D":
            name, site = line[3:11].strip(), line[13:21].strip()
            if name not in harmonics or site not in sites:
                raise ValueError("Displacement references an undefined harmonic or site")
            if site == site_name:
                if name in amplitudes:
                    raise ValueError(f"Duplicate displacement {name} for {site}")
                amplitudes[name] = [
                    float(line[a:b])
                    for a, b in ((24, 32), (33, 41), (42, 50), (53, 61), (62, 70), (71, 79))
                ]
        else:
            raise ValueError(f"Unexpected HARPOS record {kind}")
    if (
        set(harmonics) != set(CONSTITUENTS)
        or set(amplitudes) != set(CONSTITUENTS)
        or radii != [3000.0]
    ):
        raise ValueError(
            "Expected exactly 17 defined constituents, all selected-site rows, and 3000 m radius"
        )
    xyz = np.asarray(nominal_station_xyz, dtype=float).reshape(3)
    provider_xyz = sites[site_name]
    separation = float(np.linalg.norm(xyz - provider_xyz))
    if not np.all(np.isfinite(xyz)) or separation > radii[0]:
        raise ValueError(f"Station lies outside {site_name}'s {radii[0]:g} m applicability radius")

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
        "source_provenance": provenance,
        "site": site_name,
        "provider_xyz_m": provider_xyz.tolist(),
        "nominal_xyz_m": xyz.tolist(),
        "site_separation_m": separation,
        "applicability_radius_m": radii[0],
        "constituents": list(CONSTITUENTS),
        "argument_time": "TT seconds since J2000; selected Earth rotation converter at geocentre",
        "amplitude_axes": "geocentric radial Up, East, North (HARPOS specification)",
        "uen_to_itrs": basis.tolist(),
        "origin": "whole-Earth centre of mass",
        "earth_response": "STW105/REAR",
        "scope": "17 supplied gravity-tide terms including nodal sidebands; no additional terms",
        "comparison": "Different realization from ASI 11-constituent BLQ+HARDISP and JPL FES2004; not exact provider replay",
    }
    return displacement, metadata
