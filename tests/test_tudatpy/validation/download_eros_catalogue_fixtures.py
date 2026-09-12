"""Explicitly refresh the Eros fixtures from MPC and JPL (requires network access).

    python tests/test_tudatpy/validation/download_eros_catalogue_fixtures.py --output /tmp/eros-fixtures

This acquisition script does not import Tudat. It independently samples the JPL
bias map, allowing the residual test to detect changes in the reader's sign.
"""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.time import Time
from astropy_healpix import HEALPix
from astroquery.jplhorizons import Horizons
from astroquery.mpc import MPC

CATALOGUES = "a b c d e g i j l m n o p q r t u v w L N Q R S U Y".split()
CASES = (("T05", 2016), ("711", 1996))


def sample_bias(table, bias_map):
    """Read positive RA*cos(DEC), DEC biases and proper motions from the JPL map."""
    pixels = HEALPix(nside=64, order="ring").lonlat_to_healpix(
        table.RA.to_numpy() * u.deg, table.DEC.to_numpy() * u.deg
    )
    catalogue_ids = np.array([CATALOGUES.index(code) for code in table.catalog])
    entries = bias_map[pixels, catalogue_ids]
    years = (table.epoch.to_numpy() - 2451545.0) / 365.25
    return entries[:, :2] + years[:, None] * entries[:, 2:] / 1000.0


def download_case(table, observatories, bias_map, bias_hash, station, year, output):
    """Keep all catalogue-coded CCD rows for one station/year and freeze both Horizons products."""
    years = Time(table.epoch, format="jd", scale="utc").datetime
    mask = (
        (table.observatory == station)
        & np.array([date.year == year for date in years])
        & (table.note2 == "C")
        & table.catalog.isin(CATALOGUES)
    )
    selected = table.loc[
        mask, ["number", "epoch", "RA", "DEC", "note2", "catalog", "observatory"]
    ].copy()
    selected = selected.sort_values("epoch", kind="stable").reset_index(drop=True)
    if selected.empty:
        raise ValueError(f"No catalogue-coded CCD observations for {station}/{year}")
    case = output / f"eros_{station.lower()}_{year}"
    case.mkdir(parents=True, exist_ok=True)
    selected[["bias_ra_cosdec_arcsec", "bias_dec_arcsec"]] = sample_bias(selected, bias_map)

    # Observer-table quantity 1 is astrometric ICRF at UTC reception epochs, not apparent of-date.
    columns = {
        "RA": "horizons_ra_deg",
        "DEC": "horizons_dec_deg",
        "RA_3sigma": "RA_3sigma",
        "DEC_3sigma": "DEC_3sigma",
        "delta": "delta",
        "elong": "elong",
    }
    for start in range(0, len(selected), 50):
        query = Horizons(
            id="433;", location=station, epochs=selected.epoch.iloc[start : start + 50].tolist()
        )
        query.TIMEOUT = 60
        reference = query.ephemerides(quantities="1,20,23,36", extra_precision=True)
        for source, destination in columns.items():
            selected.loc[start : start + len(reference) - 1, destination] = np.asarray(
                reference[source]
            )
    selected.to_csv(case / "observations.csv", index=False)

    # Query all four translations with identical TDB grids and conventions; never integrate dynamics.
    start_jd = np.floor(selected.epoch.min()) - 2
    stop_jd = np.ceil(selected.epoch.max()) + 2
    epochs = {
        "start": Time(start_jd, format="jd", scale="tdb").iso,
        "stop": Time(stop_jd, format="jd", scale="tdb").iso,
        "step": "1h",
    }
    vectors = {}
    for name, identifier in (("Eros", "433;"), ("Earth", "399"), ("Sun", "10"), ("Moon", "301")):
        print(f"Downloading {station}/{year}: {name}", flush=True)
        query = Horizons(id=identifier, location="@0", epochs=epochs)
        query.TIMEOUT = 60
        data = query.vectors(refplane="earth", aberrations="geometric")
        times = (np.asarray(data["datetime_jd"]) - 2451545.0) * 86400.0
        if "epochs_tdb" in vectors:
            np.testing.assert_array_equal(times, vectors["epochs_tdb"])
        vectors["epochs_tdb"] = times
        state = np.column_stack(
            [np.asarray(data[column]) for column in ("x", "y", "z", "vx", "vy", "vz")]
        )
        state[:, :3] *= u.au.to(u.m)
        state[:, 3:] *= (u.au / u.day).to(u.m / u.s)
        vectors[name] = state
    np.savez_compressed(case / "horizons_vectors.npz", **vectors)

    # MPC parallax constants are geocentric: form Cartesian ITRS coordinates directly.
    site = observatories[observatories["Code"] == station][0]
    longitude = np.deg2rad(float(site["Longitude"]))
    position = 6378137.0 * np.array(
        [
            float(site["cos"]) * np.cos(longitude),
            float(site["cos"]) * np.sin(longitude),
            float(site["sin"]),
        ]
    )
    metadata = {
        "target": "433 Eros",
        "station": station,
        "station_name": str(site["Name"]),
        "station_position_itrs_m": position.tolist(),
        "station_source": "MPC longitude and geocentric parallax constants; equatorial radius 6378137 m",
        "number_of_observations": len(selected),
        "start_utc": Time(selected.epoch.min(), format="jd", scale="utc").isot,
        "end_utc": Time(selected.epoch.max(), format="jd", scale="utc").isot,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "selection": "All catalogue-coded CCD observations from this observatory/year; no residual-based rejection.",
        "bias_model": "EFCC18 / Eggl et al. 2020, debias_2018, NSIDE 64, RING",
        "bias_file_sha256": bias_hash,
        "horizons_vectors": "Geometric SSB, ICRF/J2000 equator, TDB seconds since J2000, m and m/s; hourly with two-day padding",
        "horizons_observer_reference": "Quantity 1, MPC site, UTC, extra precision; astrometric ICRF/light-time only",
        "sources": [
            "https://minorplanetcenter.net/",
            "https://ssd.jpl.nasa.gov/horizons/manual.html",
            "https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_2018.tgz",
        ],
    }
    (case / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--bias-file",
        type=Path,
        default=Path.home() / ".tudat/resource/star_catalog_biases/debias_2018/bias.dat",
    )
    args = parser.parse_args()
    table = MPC.get_observations(433).to_pandas()
    observatories = MPC.get_observatory_codes()
    bias_map = np.loadtxt(args.bias_file, comments="!").reshape(49152, 26, 4)
    bias_hash = hashlib.sha256(args.bias_file.read_bytes()).hexdigest()
    for station, year in CASES:
        download_case(table, observatories, bias_map, bias_hash, station, year, args.output)


if __name__ == "__main__":
    main()
