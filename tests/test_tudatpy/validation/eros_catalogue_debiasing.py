"""Compare both catalogue-bias signs using Eros astrometry and Horizons ephemerides.

Run with the Python package/kernel built from this PR, for example:
    python tests/test_tudatpy/validation/eros_catalogue_debiasing.py --output /tmp/eros-check

The default run uses frozen MPC observations and geometric Horizons state tables.
It computes observations with Tudat; it never integrates equations of motion.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tudatpy import kernel
from tudatpy.data_input import resource_paths
from tudatpy.data_input.environment_data import spice
from tudatpy.data_input.tracking_data.optical_utilities import read_optical_data
from tudatpy.data_input.tracking_data.optical_utilities.optical_utilities import BIAS_LOWRES_FILE
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import links, model_settings, light_time_corrections
from tudatpy.estimation.observations_setup import observations_simulation_settings

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
CASES = ("eros_t05_2016", "eros_711_1996")
ARCSEC_PER_RADIAN = 180.0 * 3600.0 / np.pi


def file_hash(path):
    """Record the exact resource/kernel used for the numerical experiment."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def create_bodies(case, stride=1, tides=True, station_offset=None):
    """Use Horizons translations, IERS Earth orientation, WGS84, and station deformation."""
    setup = environment_setup
    spice.load_standard_kernels()
    settings = setup.get_default_body_settings(["Earth", "Sun", "Moon"], "SSB", "J2000")
    settings.add_empty_settings("Eros")
    with np.load(case / "horizons_vectors.npz") as vectors:
        epochs = vectors["epochs_tdb"][::stride]
        for name in ("Eros", "Earth", "Sun", "Moon"):
            settings.get(name).ephemeris_settings = setup.ephemeris.tabulated(
                dict(zip(epochs, vectors[name][::stride])), "SSB", "J2000"
            )

    earth = settings.get("Earth")
    earth.shape_settings = setup.shape.oblate_spherical(6378137.0, 1.0 / 298.257223563)
    earth.rotation_model_settings = setup.rotation_model.gcrs_to_itrs(
        setup.rotation_model.iau_2006, "J2000"
    )
    earth.gravity_field_settings.associated_reference_frame = "ITRS"
    if tides:
        earth.shape_deformation_settings = [
            setup.shape_deformation.iers_2010_solid_body_tidal(),
            setup.shape_deformation.pole_tidal(),
        ]
    metadata = json.loads((case / "metadata.json").read_text())
    position = np.array(metadata["station_position_itrs_m"])
    if station_offset is not None:
        position += np.asarray(station_offset)
    motions = [setup.ground_station.BodyDeformationStationMotionSettings()] if tides else []
    earth.ground_station_settings = [
        setup.ground_station.basic_station(
            metadata["station"], position, station_motion_settings=motions
        )
    ]
    return setup.create_system_of_bodies(settings)


def run_case(case, stride=1, tides=True, station_offset=None, relativity=True):
    """Compute raw, PR-subtracted, and reversed-sign residuals through the same Tudat model."""
    case = Path(case)
    table = pd.read_csv(
        case / "observations.csv", dtype={"number": str, "observatory": str, "catalog": str}
    )
    table = table.sort_values("epoch", kind="stable").reset_index(drop=True)
    bodies = create_bodies(case, stride, tides, station_offset)
    tracking, _ = read_optical_data(
        table,
        custom_name="Eros",
        add_weights=False,
        add_star_catalog_corrections=True,
        add_ancillary_data=True,
    )
    if len(tracking) != 1:
        raise ValueError("Each case must contain one target and one ground station.")
    settings = []
    station = table.observatory.iloc[0]
    link = links.link_definition(
        {
            links.transmitter: links.body_origin_link_end_id("Eros"),
            links.receiver: links.body_reference_point_link_end_id("Earth", station),
        }
    )
    correction = (
        [light_time_corrections.first_order_relativistic_light_time_correction(["Sun"])]
        if relativity
        else []
    )
    settings.append(model_settings.angular_position(link, correction))
    simulators = observations_simulation_settings.create_observation_simulators(settings, bodies)

    raw_angles = np.deg2rad(table[["RA", "DEC"]].to_numpy())
    raw_angles[:, 0] = (raw_angles[:, 0] + np.pi) % (2.0 * np.pi) - np.pi
    projection = np.column_stack([np.cos(raw_angles[:, 1]), np.ones(len(table))])
    stored_corrections = np.asarray(tracking[0].get_observation_corrections()).reshape(-1, 2)
    result = table.copy()
    for name in ("raw", "subtract", "add"):
        if name == "add":
            tracking[0].set_observation_corrections((-stored_corrections).tolist())
        collection = observations.create_observation_collection_from_tracking_data(
            tracking, bodies, apply_corrections=name != "raw"
        )
        observations.compute_residuals_and_dependent_variables(collection, simulators, bodies)
        residuals = np.asarray(collection.get_concatenated_residuals()).reshape(-1, 2)
        # Wrap RA at the meridian, then express both components in the tangent plane.
        residuals[:, 0] = (residuals[:, 0] + np.pi) % (2.0 * np.pi) - np.pi
        residuals *= ARCSEC_PER_RADIAN * projection
        result[[f"{name}_ra_arcsec", f"{name}_dec_arcsec"]] = residuals
        if name == "raw":
            # Ensure the converter preserved the input rows and coordinates before comparing signs.
            np.testing.assert_allclose(
                np.asarray(collection.concatenated_observations).reshape(-1, 2),
                raw_angles,
                atol=1e-14,
                rtol=0,
            )

    result[["applied_ra_arcsec", "applied_dec_arcsec"]] = (
        stored_corrections * ARCSEC_PER_RADIAN * projection
    )
    reference = raw_angles - np.deg2rad(table[["horizons_ra_deg", "horizons_dec_deg"]].to_numpy())
    reference[:, 0] = (reference[:, 0] + np.pi) % (2.0 * np.pi) - np.pi
    result[["horizons_raw_ra_arcsec", "horizons_raw_dec_arcsec"]] = (
        reference * ARCSEC_PER_RADIAN * projection
    )
    return result


def statistics(result):
    """Summarize both residual components and their combined RMS without rejecting data."""
    summary = {"count": len(result)}
    for name in ("raw", "subtract", "add"):
        values = result[[f"{name}_ra_arcsec", f"{name}_dec_arcsec"]].to_numpy()
        summary[name] = {
            "rms_arcsec": float(np.sqrt(np.mean(np.sum(values**2, axis=1)))),
            "mean_ra_arcsec": float(values[:, 0].mean()),
            "mean_dec_arcsec": float(values[:, 1].mean()),
        }
    applied = result[["applied_ra_arcsec", "applied_dec_arcsec"]].to_numpy()
    summary["correction_rms_arcsec"] = float(np.sqrt(np.mean(np.sum(applied**2, axis=1))))
    summary["correction_to_raw_rms"] = (
        summary["correction_rms_arcsec"] / summary["raw"]["rms_arcsec"]
    )
    model_difference = (
        result[["raw_ra_arcsec", "raw_dec_arcsec"]].to_numpy()
        - result[["horizons_raw_ra_arcsec", "horizons_raw_dec_arcsec"]].to_numpy()
    )
    summary["horizons_observer_difference_maximum_arcsec"] = float(
        np.max(np.linalg.norm(model_difference, axis=1))
    )
    return summary


def plot_case(result, path):
    """Show every retained measurement for all three signs on identical time and angle axes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.time import Time

    dates = Time(result.epoch, format="jd", scale="utc").to_datetime()
    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 6), constrained_layout=True)
    for axis, component, label in zip(
        axes, ("ra", "dec"), (r"$\Delta\alpha\cos\delta$", r"$\Delta\delta$")
    ):
        for name, color in (("raw", "0.45"), ("subtract", "tab:blue"), ("add", "tab:red")):
            axis.scatter(dates, result[f"{name}_{component}_arcsec"], s=14, label=name, color=color)
        axis.axhline(0, color="black", linewidth=0.5)
        axis.set_ylabel(label + " [arcsec]")
        axis.grid(alpha=0.2)
    axes[0].legend(ncol=3)
    axes[0].set_title(f"Eros — MPC {result.observatory.iloc[0]} — {len(result)} observations")
    axes[1].set_xlabel("Reception epoch (UTC)")
    figure.savefig(path)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    eop_file = Path(resource_paths.get_earth_orientation_path()) / "eopc04_14_IAU2000.62-now.txt"
    summary = {
        "kernel_path": kernel.__file__,
        "kernel_sha256": file_hash(kernel.__file__),
        "eop_path": str(eop_file),
        "eop_sha256": file_hash(eop_file),
        "bias_sha256": file_hash(BIAS_LOWRES_FILE),
        "cases": {},
    }
    for name in args.cases:
        print(f"Computing {name} with {kernel.__file__}", flush=True)
        result = run_case(FIXTURES / name)
        result.to_csv(args.output / f"{name}_residuals.csv", index=False)
        plot_case(result, args.output / f"{name}_residuals.pdf")
        info = statistics(result)
        # Quantify interpolation, Earth deformation, and first-order light-time model sensitivity.
        baseline = result[["raw_ra_arcsec", "raw_dec_arcsec"]].to_numpy()
        for label, options in (
            ("two_hour_grid", {"stride": 2}),
            ("without_tides", {"tides": False}),
            ("without_solar_shapiro", {"relativity": False}),
            ("station_shift_100m", {"station_offset": [100.0, 0, 0]}),
        ):
            alternative = run_case(FIXTURES / name, **options)
            delta = alternative[["raw_ra_arcsec", "raw_dec_arcsec"]].to_numpy() - baseline
            info[label + "_maximum_arcsec"] = float(np.max(np.linalg.norm(delta, axis=1)))
        summary["cases"][name] = info
        print(json.dumps(info, indent=2), flush=True)
        (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
