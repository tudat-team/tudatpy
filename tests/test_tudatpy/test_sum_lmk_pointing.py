"""Python-level tests for loading SPC SUM/LMK files as tracking data, and for the per-image camera
pointing correction.

The SUM/LMK fixtures used by the C++ tests live inside ``tests/test_tudat/data.zip`` and are only
unpacked into the build tree, so they are not reachable from an installed tudatpy. The files are
small and rigidly formatted, so they are written out here instead, which also keeps this test
independent of the C++ test data.
"""

import numpy as np
import pytest

from tudatpy.data_input.tracking_data import sum_lmk
from tudatpy.dynamics import environment_setup, parameters_setup
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import model_settings
from tudatpy.estimation.observations_setup import observations_simulation_settings
from tudatpy.data_input.environment_data import spice


def _fortran_double(value: float) -> str:
    """Format a float the way the SPC SUM/LMK files do (``0.1000000000D+03``)."""
    mantissa, exponent = f"{value:.10E}".split("E")
    mantissa = float(mantissa) / 10.0
    return f"{mantissa:.10f}D{int(exponent) + 1:+03d}"


def _write_sum_file(directory, image_id, utc, landmark_pixels, pointing_sigma):
    """Write one SUM image: spacecraft 10 km behind the target, identity camera axes."""
    lines = [
        image_id,
        utc,
        "  1024  1024   500 65535                         NPX, NLN, THRSH",
        f"{_fortran_double(100.0)} {_fortran_double(512.0)} {_fortran_double(512.0)}   MMFL, CTR",
        f"{_fortran_double(0.0)} {_fortran_double(0.0)} {_fortran_double(10.0)}   SCOBJ",
        f"{_fortran_double(1.0)} {_fortran_double(0.0)} {_fortran_double(0.0)}   CX",
        f"{_fortran_double(0.0)} {_fortran_double(1.0)} {_fortran_double(0.0)}   CY",
        f"{_fortran_double(0.0)} {_fortran_double(0.0)} {_fortran_double(1.0)}   CZ",
        f"{_fortran_double(0.0)} {_fortran_double(0.0)} {_fortran_double(1.0)}   SZ",
        " ".join(_fortran_double(v) for v in (10.0, 0.0, 0.0, 0.0, 10.0, 0.0)) + " K-MATRIX",
        f"{_fortran_double(pointing_sigma)} SIGMA_PTG",
        "LANDMARKS",
    ]
    for landmark_id, (sample, line) in landmark_pixels.items():
        lines.append(f"{landmark_id} {sample:.2f} {line:.2f}")
    lines.append("END FILE")

    path = directory / f"{image_id}.SUM"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _write_lmk_file(directory, landmark_id, body_fixed_position, image_ids):
    lines = [
        f"{landmark_id} T",
        "64 SIZE",
        f"{_fortran_double(0.1)} SCALE",
        "0 0 0 0 HORIZON",
        f"{_fortran_double(1.0e-3)} SIGKM",
        f"{_fortran_double(2.0e-3)} RMSLMK",
        " ".join(_fortran_double(v) for v in body_fixed_position) + " VLM",
        f"{_fortran_double(1.0)} {_fortran_double(0.0)} {_fortran_double(0.0)} UX",
        f"{_fortran_double(0.0)} {_fortran_double(1.0)} {_fortran_double(0.0)} UY",
        f"{_fortran_double(0.0)} {_fortran_double(0.0)} {_fortran_double(1.0)} UZ",
        f"{_fortran_double(1.0e-3)} SIGMA_LMK",
        "PICTURES",
    ]
    for image_id in image_ids:
        lines.append(f"{image_id} 512.00 512.00 *")
    lines.append("END FILE")

    path = directory / f"{landmark_id}.LMK"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# Landmark body-fixed positions [km], spread so they project away from the optical centre.
LANDMARK_POSITIONS = {
    "LMK0001": (0.0, 0.0, 0.0),
    "LMK0002": (0.3, -0.2, 0.05),
}
IMAGE_IDS = ("IMG001", "IMG002")
POINTING_SIGMA = 1.0e-3


@pytest.fixture(scope="module")
def sum_lmk_files(tmp_path_factory):
    directory = tmp_path_factory.mktemp("sum_lmk")
    sum_files = [
        _write_sum_file(
            directory,
            "IMG001",
            "2015 JUN 05 07:24:42.053",
            {"LMK0001": (512.0, 512.0), "LMK0002": (542.0, 492.0)},
            POINTING_SIGMA,
        ),
        _write_sum_file(
            directory,
            "IMG002",
            "2015 JUN 05 07:34:42.053",
            {"LMK0001": (512.0, 512.0), "LMK0002": (542.0, 492.0)},
            POINTING_SIGMA,
        ),
    ]
    lmk_files = [
        _write_lmk_file(directory, landmark_id, position, IMAGE_IDS)
        for landmark_id, position in LANDMARK_POSITIONS.items()
    ]
    return sum_files, lmk_files


def _make_bodies():
    """Target at the origin, spacecraft 10 km behind it, both with identity attitude."""
    spice.load_standard_kernels()
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    states = {
        "Target": np.zeros(6),
        "Spacecraft": np.array([0.0, 0.0, -10000.0, 0.0, 0.0, 0.0]),
    }
    for name, state in states.items():
        body_settings.add_empty_settings(name)
        settings = body_settings.get(name)
        settings.ephemeris_settings = environment_setup.ephemeris.constant(state, "SSB", "J2000")
        settings.rotation_model_settings = environment_setup.rotation_model.constant_rotation_model(
            "J2000", f"{name}_Fixed", np.eye(3)
        )
    return environment_setup.create_system_of_bodies(body_settings)


def _camera_data(supplementary_data):
    """Map camera name to camera supplementary data."""
    return {
        camera.camera_id: camera
        for entry in supplementary_data
        for camera in entry.instrument_supplementary_data
    }


@pytest.fixture
def loaded(sum_lmk_files):
    """Read the files, apply the supplementary data and create the observation collection."""
    sum_files, lmk_files = sum_lmk_files
    tracking_data, supplementary_data = sum_lmk.read_sum_lmk_data(
        sum_files, lmk_files, "Target", "Spacecraft"
    )
    bodies = _make_bodies()
    observations.set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    collection = observations.create_observation_collection_from_tracking_data(
        tracking_data, bodies
    )
    return bodies, tracking_data, supplementary_data, collection


def _residuals(bodies, collection):
    link_definitions = collection.get_link_definitions_for_observables(
        model_settings.pixel_coordinates_type
    )
    simulators = observations_simulation_settings.create_observation_simulators(
        [model_settings.pixel_coordinates(link) for link in link_definitions], bodies
    )
    observations.compute_residuals_and_dependent_variables(collection, simulators, bodies)
    return np.array(collection.get_concatenated_residuals())


def test_reading_needs_no_bodies(sum_lmk_files):
    sum_files, lmk_files = sum_lmk_files
    tracking_data, supplementary_data = sum_lmk.read_sum_lmk_data(
        sum_files, lmk_files, "Target", "Spacecraft"
    )

    # One observation per (image, landmark) pair.
    assert len(tracking_data) == len(IMAGE_IDS) * len(LANDMARK_POSITIONS)
    for entry in tracking_data:
        assert entry.observable_type == "PixelCoordinates"
        assert entry.time_scale == "UTC"

    landmarks = {
        entry.reference_point_name: entry.reference_point_supplementary_data
        for entry in supplementary_data
        if entry.reference_point_supplementary_data is not None
    }
    assert set(landmarks) == set(LANDMARK_POSITIONS)
    np.testing.assert_allclose(
        landmarks["LMK0002"].body_fixed_position,
        1000.0 * np.array(LANDMARK_POSITIONS["LMK0002"]),
        rtol=1e-12,
    )

    cameras = _camera_data(supplementary_data)
    assert set(cameras) == {f"Camera_{image_id}" for image_id in IMAGE_IDS}
    for camera in cameras.values():
        assert camera.pointing_reference_body_name == "Target"
        np.testing.assert_allclose(camera.pointing_uncertainty, POINTING_SIGMA, rtol=1e-12)


def test_raw_contents_convert_like_files(sum_lmk_files):
    sum_files, lmk_files = sum_lmk_files
    images = sum_lmk.read_raw_sum_files(sum_files)
    landmarks = sum_lmk.read_raw_lmk_files(lmk_files)
    assert {image.image_id for image in images} == set(IMAGE_IDS)
    assert set(landmarks) == set(LANDMARK_POSITIONS)

    tracking_data, supplementary_data = sum_lmk.convert_raw_sum_lmk_data(
        images, landmarks, "Target", "Spacecraft"
    )
    assert len(tracking_data) == len(IMAGE_IDS) * len(LANDMARK_POSITIONS)
    assert len(supplementary_data) == len(IMAGE_IDS) + len(LANDMARK_POSITIONS)


def test_supplementary_data_adds_cameras_and_landmarks(loaded):
    bodies, _, supplementary_data, collection = loaded

    # Each image gets its own camera, so each image gets its own pointing parameter.
    camera_dict = bodies.get("Spacecraft").system_models.camera_dict
    for camera_name in _camera_data(supplementary_data):
        assert camera_name in camera_dict

    ground_stations = bodies.get("Target").ground_station_list
    for landmark_id in LANDMARK_POSITIONS:
        assert landmark_id in ground_stations

    assert collection.concatenated_observations.size == 2 * len(IMAGE_IDS) * len(LANDMARK_POSITIONS)


def test_pointing_parameters_and_apriori(loaded):
    bodies, _, supplementary_data, _ = loaded
    cameras = _camera_data(supplementary_data)

    parameter_settings = [
        parameters_setup.camera_pointing_correction("Spacecraft", camera_name)
        for camera_name in sorted(cameras)
    ]
    parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
    assert parameters.parameter_set_size == 3 * len(IMAGE_IDS)

    # The SUM SIGMA_PTG values give a diagonal a-priori for the pointing parameters.
    inverse_apriori = np.zeros((parameters.parameter_set_size,) * 2)
    for camera_name, camera in cameras.items():
        parameter_type = (
            parameters_setup.EstimatableParameterTypes.camera_pointing_correction_type,
            ("Spacecraft", camera_name),
        )
        ((start, size),) = parameters.indices_for_parameter_type(parameter_type)
        block = slice(start, start + size)
        inverse_apriori[block, block] = np.diag(1.0 / camera.pointing_uncertainty**2)

    np.testing.assert_allclose(np.diag(inverse_apriori), 1.0 / POINTING_SIGMA**2, rtol=1e-12)


def test_pointing_correction_changes_the_camera_orientation(loaded):
    bodies, *_ = loaded
    camera = bodies.get("Spacecraft").system_models.get_camera("Camera_IMG001")

    np.testing.assert_allclose(camera.pointing_correction, np.zeros(3), atol=0.0)
    nominal = np.array(camera.nominal_quaternion)
    np.testing.assert_allclose(camera.quaternion, nominal, atol=1e-15)

    correction = np.array([1.0e-4, -2.0e-4, 3.0e-4])
    camera.pointing_correction = correction
    np.testing.assert_allclose(camera.pointing_correction, correction, rtol=1e-12)

    # The nominal orientation is untouched, the effective orientation moves.
    np.testing.assert_allclose(camera.nominal_quaternion, nominal, atol=1e-15)
    assert np.linalg.norm(np.array(camera.quaternion) - nominal) > 1e-6

    camera.pointing_correction = np.zeros(3)
    np.testing.assert_allclose(camera.quaternion, nominal, atol=1e-15)


def test_pointing_correction_moves_pixel_residuals(loaded):
    bodies, _, _, collection = loaded

    baseline = _residuals(bodies, collection)
    assert baseline.size > 0

    camera = bodies.get("Spacecraft").system_models.get_camera("Camera_IMG001")
    camera.pointing_correction = np.array([1.0e-3, 0.0, 0.0])
    perturbed = _residuals(bodies, collection)

    # The correction must actually reach the observation model.
    assert np.linalg.norm(perturbed - baseline) > 1e-3

    camera.pointing_correction = np.zeros(3)
    restored = _residuals(bodies, collection)
    np.testing.assert_allclose(restored, baseline, atol=1e-10)
