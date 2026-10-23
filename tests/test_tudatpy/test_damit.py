from pathlib import Path

import numpy as np
import pytest

import tudatpy.data_input.environment_data.damit.damit as damit_module
from tudatpy.astro.element_conversion import j2000_to_eclipj2000
from tudatpy.data_input.environment_data.damit import DamitModel
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup.shape import TriAxialEllipsoidBodyShapeSettings

DATA_DIR = Path(__file__).parent
EROS_SPIN = DATA_DIR / "eros_spin.txt"
EROS_SHAPE = DATA_DIR / "eros_shape.obj"

TEST_ASTEROID = 433
TEST_ASTEROID_MODEL = 3083


def _sphere_volume(diameter):
    return 4 / 3 * np.pi * (diameter / 2) ** 3


def _box_model(x, y, z, angle=0.0):
    """DamitModel of a box with half-sides x, y, z, rotated by `angle` around the z-axis"""
    vertices = np.array(
        [
            [-x, -y, -z],
            [x, -y, -z],
            [x, y, -z],
            [-x, y, -z],
            [-x, -y, z],
            [x, -y, z],
            [x, y, z],
            [-x, y, z],
        ],
        dtype=float,
    )
    # Outward-facing triangles
    facets = np.array(
        [
            [0, 2, 1],
            [0, 3, 2],  # -z
            [4, 5, 6],
            [4, 6, 7],  # +z
            [0, 1, 5],
            [0, 5, 4],  # -y
            [3, 7, 6],
            [3, 6, 2],  # +y
            [0, 4, 7],
            [0, 7, 3],  # -x
            [1, 6, 5],
            [1, 2, 6],  # +x
        ]
    )
    rotation_z = np.array(
        [[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
    )
    return DamitModel(
        name="box",
        vertices=vertices @ rotation_z.T,
        facets=facets,
        pole_longitude=0.0,
        pole_latitude=np.pi / 2,
        period=3600.0,
        initial_rotation_angle=0.0,
        reference_epoch=0.0,
    )


@pytest.fixture(scope="module")
def eros():
    return DamitModel.from_file(EROS_SPIN, EROS_SHAPE, "Eros")


def test_from_file_parses_spin_and_shape(eros):
    """Shape and spin files are read correctly and converted to SI units / radians"""
    assert eros._vertices.shape == (856, 3)
    assert eros._facets.shape == (1708, 3)
    assert eros._facets.min() == 0
    assert eros._facets.max() == 855

    assert eros._pole_longitude == pytest.approx(np.deg2rad(17))
    assert eros._pole_latitude == pytest.approx(np.deg2rad(11))
    assert eros._period == pytest.approx(5.27025528 * 3600)
    assert eros._initial_rotation_angle == pytest.approx(np.deg2rad(32.64))
    assert eros._reference_epoch == pytest.approx(0.0)  # JD 2451545.0 is J2000


def test_volume_scaling(eros):
    """Without diameter the shape is converted from km to m, with diameter it is scaled to that size"""
    assert eros.volume == pytest.approx(2.4916158e12, rel=1e-6)

    diameter = 16_840.0
    scaled = DamitModel.from_file(EROS_SPIN, EROS_SHAPE, "Eros", diameter=diameter)
    assert scaled.volume == pytest.approx(_sphere_volume(diameter))


def test_non_obj_shape_file_raises(tmp_path):
    """A shape file not in OBJ format is rejected"""
    shape_file = tmp_path / "shape.txt"
    shape_file.write_text("4 4\n0 0 0\n")
    with pytest.raises(ValueError):
        DamitModel.from_file(EROS_SPIN, shape_file, "test")


def test_ellipsoid_semi_axes_box():
    """A box with half-sides (x, y, z) has an equal-inertia ellipsoid with semi-axes sqrt(5/3) * (x, y, z)"""
    half_sides = np.array([3.0, 2.0, 1.0])
    model = _box_model(*half_sides)

    semi_axes, offset = model._approximate_best_fitting_ellipsoid()

    np.testing.assert_allclose(semi_axes, np.sqrt(5 / 3) * half_sides, rtol=1e-10)
    assert offset == pytest.approx(0.0, abs=1e-10)
    assert isinstance(
        model.best_fitting_ellipsoid_shape_model(), TriAxialEllipsoidBodyShapeSettings
    )


@pytest.mark.parametrize("angle", [0.3, -0.7, 1.2])
def test_ellipsoid_offset_rotated_box(angle):
    """Rotating the shape around the z-axis is recovered as the principal axis offset"""
    half_sides = np.array([3.0, 2.0, 1.0])
    model = _box_model(*half_sides, angle=angle)

    semi_axes, offset = model._approximate_best_fitting_ellipsoid()

    np.testing.assert_allclose(semi_axes, np.sqrt(5 / 3) * half_sides, rtol=1e-10)
    assert offset == pytest.approx(angle, abs=1e-10)


@pytest.mark.parametrize("base_frame", ["ECLIPJ2000", "J2000"])
def test_rotation_model_pole_and_period(eros, base_frame):
    """The body z-axis points along the pole, and the orientation repeats after one period"""
    body_settings = environment_setup.BodyListSettings("SSB", base_frame)
    body_settings.add_empty_settings("Eros")
    body_settings.get("Eros").rotation_model_settings = (
        eros.polyhedron_rotational_ephemeris_settings(base_frame)
    )
    rotation = environment_setup.create_system_of_bodies(body_settings).get("Eros").rotation_model

    longitude, latitude = eros._pole_longitude, eros._pole_latitude
    pole = np.array(
        [
            np.cos(latitude) * np.cos(longitude),
            np.cos(latitude) * np.sin(longitude),
            np.sin(latitude),
        ]
    )
    if base_frame == "J2000":
        pole = j2000_to_eclipj2000().T @ pole

    epoch = eros._reference_epoch
    rotation_at_epoch = rotation.body_fixed_to_inertial_rotation(epoch)
    np.testing.assert_allclose(rotation_at_epoch[:, 2], pole, atol=1e-12)
    np.testing.assert_allclose(
        rotation.body_fixed_to_inertial_rotation(epoch + eros._period), rotation_at_epoch, atol=1e-9
    )


def test_rotation_invalid_base_frame(eros):
    """Only J2000 and ECLIPJ2000 are supported base frames"""
    with pytest.raises(ValueError):
        eros.polyhedron_rotational_ephemeris_settings("IAU_EARTH")


@pytest.mark.parametrize(
    "kwargs",
    [{}, {"mpc_number": TEST_ASTEROID, "damit_model_id": TEST_ASTEROID_MODEL}],
    ids=["none", "both"],
)
def test_query_argument_validation(kwargs):
    """Exactly one of MPC number and model ID must be given"""
    with pytest.raises(ValueError):
        DamitModel.query(**kwargs)


@pytest.fixture
def offline_damit(tmp_path, monkeypatch):
    """Replace the DAMIT tables by small local tables, and model downloads by the local Eros files"""
    asteroids = tmp_path / "asteroids.csv"
    asteroids.write_text("id,number,name\n1,433,Eros\n2,1000,Piazzia\n")
    asteroid_models = tmp_path / "asteroid_models.csv"
    asteroid_models.write_text(
        "id,asteroid_id,lsm,calibrated_size,equiv_diameter,quality_flag\n"
        "10,1,LS,0,16.84,2\n"  # low quality
        "11,1,LS,0,20.0,4\n"  # best quality for 433
        "12,1,LS,1,,1\n"  # calibrated size
        "20,2,LS,0,,3\n"  # no size information
    )  # Text column as in the real table, otherwise pandas turns the selected model ID into a float
    monkeypatch.setattr(damit_module, "_DAMIT_ASTEROIDS", asteroids)
    monkeypatch.setattr(damit_module, "_DAMIT_ASTEROID_MODELS", asteroid_models)
    monkeypatch.setattr(
        damit_module, "_download_model", lambda model_id, force=False: (EROS_SHAPE, EROS_SPIN)
    )


def test_query_selection_offline(offline_damit, eros):
    """Model selection and sizing from the DAMIT tables"""
    # Highest quality flag is selected, and scaled to the catalog diameter (km -> m)
    best = DamitModel.query(mpc_number=433)
    assert best._name == "11"
    assert best.volume == pytest.approx(_sphere_volume(20_000.0))

    # Given diameter overrides the catalog diameter
    assert DamitModel.query(damit_model_id=10, diameter=1_000.0).volume == pytest.approx(
        _sphere_volume(1_000.0)
    )

    # Calibrated model is not rescaled
    assert DamitModel.query(damit_model_id=12).volume == pytest.approx(eros.volume)

    # No size information available
    with pytest.warns(UserWarning):
        DamitModel.query(mpc_number=1000)

    # Unknown asteroid or model
    with pytest.raises(LookupError):
        DamitModel.query(mpc_number=1)
    with pytest.raises(LookupError):
        DamitModel.query(damit_model_id=99)


@pytest.mark.remote_data
def test_query_remote():
    """Querying the online database by MPC number and by model ID gives the same model"""
    model_by_mpc = DamitModel.query(mpc_number=TEST_ASTEROID, force_refresh=True)
    model_by_id = DamitModel.query(damit_model_id=TEST_ASTEROID_MODEL)

    assert model_by_mpc._name == str(TEST_ASTEROID_MODEL)
    assert model_by_mpc.volume == pytest.approx(model_by_id.volume)
