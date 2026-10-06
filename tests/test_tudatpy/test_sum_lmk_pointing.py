"""Python-level tests for the SPC SUM/LMK pixel-landmark workflow and the per-image camera
pointing correction.

The SUM/LMK fixtures used by the C++ tests live inside ``tests/test_tudat/data.zip`` and are only
unpacked into the build tree, so they are not reachable from an installed tudatpy. The files are
small and rigidly formatted, so they are written out here instead, which also keeps this test
independent of the C++ test data.
"""

import numpy as np
import pytest

from tudatpy.dynamics import environment_setup, parameters_setup
from tudatpy.estimation.observations_setup import observations_wrapper
from tudatpy.interface import spice


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


@pytest.fixture
def conversion(sum_lmk_files):
    sum_files, lmk_files = sum_lmk_files
    bodies = _make_bodies()
    settings = observations_wrapper.SumLmkObservationConversionSettings("Target", "Spacecraft")
    result = observations_wrapper.create_sum_lmk_observation_collection(
        sum_files, lmk_files, bodies, settings
    )
    return bodies, result


def test_conversion_registers_one_camera_per_image(conversion):
    bodies, result = conversion

    assert result.receiver_body_name == "Spacecraft"
    image_to_camera = dict(result.image_id_to_camera_name)
    assert set(image_to_camera) == set(IMAGE_IDS)
    # Each image gets its own camera, so each image gets its own pointing parameter.
    assert len(set(image_to_camera.values())) == len(IMAGE_IDS)

    camera_dict = bodies.get("Spacecraft").system_models.camera_dict
    for camera_name in image_to_camera.values():
        assert camera_name in camera_dict

    assert result.observation_collection.concatenated_observations.size > 0
    assert len(result.observation_model_settings) > 0


def test_pointing_parameter_settings(conversion):
    _, result = conversion

    all_settings = observations_wrapper.create_sum_lmk_pointing_parameter_settings(result)
    assert len(all_settings) == len(IMAGE_IDS)

    subset = observations_wrapper.create_sum_lmk_pointing_parameter_settings(result, ["IMG002"])
    assert len(subset) == 1

    with pytest.raises(RuntimeError):
        observations_wrapper.create_sum_lmk_pointing_parameter_settings(result, ["NOT_AN_IMAGE"])


def test_inverse_apriori_covariance_from_sigma_ptg(conversion):
    bodies, result = conversion

    parameter_settings = observations_wrapper.create_sum_lmk_pointing_parameter_settings(result)
    parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
    assert parameters.parameter_set_size == 3 * len(IMAGE_IDS)

    covariance = observations_wrapper.create_sum_lmk_inverse_apriori_covariance(result, parameters)
    assert covariance.shape == (3 * len(IMAGE_IDS), 3 * len(IMAGE_IDS))

    expected = 1.0 / POINTING_SIGMA**2
    np.testing.assert_allclose(np.diag(covariance), expected, rtol=1e-12)
    # Purely diagonal: SIGMA_PTG gives no cross-correlation.
    np.testing.assert_allclose(covariance - np.diag(np.diag(covariance)), 0.0, atol=0.0)

    # Composes with an existing a-priori rather than replacing it.
    base = np.eye(covariance.shape[0])
    combined = observations_wrapper.create_sum_lmk_inverse_apriori_covariance(
        result, parameters, base
    )
    np.testing.assert_allclose(combined, covariance + base, rtol=1e-12)


def test_apriori_skips_parameters_that_are_not_estimated(conversion):
    bodies, result = conversion

    # Only one image's pointing is estimated; the other image's a-priori has nowhere to go.
    parameter_settings = observations_wrapper.create_sum_lmk_pointing_parameter_settings(
        result, ["IMG002"]
    )
    parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
    covariance = observations_wrapper.create_sum_lmk_inverse_apriori_covariance(result, parameters)

    assert covariance.shape == (3, 3)
    np.testing.assert_allclose(np.diag(covariance), 1.0 / POINTING_SIGMA**2, rtol=1e-12)


def test_pointing_correction_changes_the_camera_orientation(conversion):
    bodies, result = conversion
    camera_name = dict(result.image_id_to_camera_name)["IMG001"]
    camera = bodies.get("Spacecraft").system_models.get_camera(camera_name)

    np.testing.assert_allclose(camera.pointing_correction, np.zeros(3), atol=0.0)
    nominal = np.array(camera.nominal_quaternion)
    np.testing.assert_allclose(camera.quaternion, nominal, atol=1e-15)

    correction = np.array([1.0e-4, -2.0e-4, 3.0e-4])
    camera.pointing_correction = correction
    np.testing.assert_allclose(camera.pointing_correction, correction, rtol=1e-12)

    # The nominal orientation is untouched; the effective orientation moves.
    np.testing.assert_allclose(camera.nominal_quaternion, nominal, atol=1e-15)
    assert np.linalg.norm(np.array(camera.quaternion) - nominal) > 1e-6

    camera.pointing_correction = np.zeros(3)
    np.testing.assert_allclose(camera.quaternion, nominal, atol=1e-15)


def test_pointing_correction_moves_pixel_residuals(conversion):
    bodies, result = conversion

    baseline = observations_wrapper.compute_sum_lmk_residuals(result.observation_collection, bodies)
    assert baseline.size > 0

    camera_name = dict(result.image_id_to_camera_name)["IMG001"]
    camera = bodies.get("Spacecraft").system_models.get_camera(camera_name)
    camera.pointing_correction = np.array([1.0e-3, 0.0, 0.0])

    perturbed = observations_wrapper.compute_sum_lmk_residuals(
        result.observation_collection, bodies
    )

    # The correction must actually reach the observation model.
    assert np.linalg.norm(perturbed - baseline) > 1e-3

    camera.pointing_correction = np.zeros(3)
    restored = observations_wrapper.compute_sum_lmk_residuals(result.observation_collection, bodies)
    np.testing.assert_allclose(restored, baseline, atol=1e-10)
