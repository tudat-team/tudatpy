from pathlib import Path

import numpy as np
import pytest

import tudatpy.data_input.environment_data.damit.damit as damit_module
from tudatpy.astro.polyhedron_utilities import volume
from tudatpy.data_input.environment_data.damit import DamitModel
from tudatpy.dynamics import environment_setup
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).parent
EROS_SPIN = DATA_DIR / "eros_spin.txt"
EROS_SHAPE = DATA_DIR / "eros_shape.obj"

TEST_ASTEROID = 433
TEST_ASTEROID_MODEL = 3083


def _sphere_volume(diameter):
    return 4 / 3 * np.pi * (diameter / 2) ** 3


def _model_volume(model):
    return volume(model.vertices, model.facets)


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
    """Shape model of (433) Eros with DAMIT Model ID 3083"""
    return DamitModel.from_file(EROS_SPIN, EROS_SHAPE)


def test_from_file_parses_spin_and_shape(eros):
    """Test parsing of shape and spin files and conversion to correct units"""
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
    """Test that shape model is correctly scaled and that diameter input scales volume to that diameter"""
    assert _model_volume(eros) == pytest.approx(2.4916158e12, rel=1e-6)

    diameter = 16_840.0
    scaled = DamitModel.from_file(EROS_SPIN, EROS_SHAPE, diameter=diameter)
    assert _model_volume(scaled) == pytest.approx(_sphere_volume(diameter))


def test_non_obj_shape_file_raises(tmp_path):
    """A shape file not in OBJ format is rejected"""
    shape_file = tmp_path / "shape.txt"
    shape_file.write_text("4 4\n0 0 0\n")
    with pytest.raises(ValueError):
        DamitModel.from_file(EROS_SPIN, shape_file)


def test_ellipsoid_semi_axes_box():
    """A box with half-sides (x, y, z) has an equal-inertia ellipsoid with semi-axes sqrt(5/3) * (x, y, z)"""
    half_sides = np.array([3.0, 2.0, 1.0])
    model = _box_model(*half_sides)

    semi_axes, offset = model.best_fitting_ellipsoid()

    np.testing.assert_allclose(semi_axes, np.sqrt(5 / 3) * half_sides, rtol=1e-10)
    assert offset == pytest.approx(0.0, abs=1e-10)


def test_shape_settings(eros):
    """Check wiring of DamitModel to shape settings creation"""
    polyhedron = environment_setup.shape.damit_polyhedron(eros)
    np.testing.assert_array_equal(polyhedron.vertices_coordinates, eros.vertices)
    np.testing.assert_array_equal(polyhedron.vertices_defining_each_facet, eros.facets)

    ellipsoid = environment_setup.shape.damit_best_fitting_ellipsoid(eros)
    semi_axes, _ = eros.best_fitting_ellipsoid()
    np.testing.assert_allclose(
        [ellipsoid.semi_axis_a, ellipsoid.semi_axis_b, ellipsoid.semi_axis_c], semi_axes
    )


@pytest.mark.parametrize("angle", [0.3, -0.7, 1.2])
def test_ellipsoid_offset_rotated_box(angle):
    """Rotating the shape around the z-axis is recovered as the principal axis offset"""
    half_sides = np.array([3.0, 2.0, 1.0])
    model = _box_model(*half_sides, angle=angle)

    semi_axes, offset = model.best_fitting_ellipsoid()

    np.testing.assert_allclose(semi_axes, np.sqrt(5 / 3) * half_sides, rtol=1e-10)
    assert offset == pytest.approx(angle, abs=1e-10)


def _passive_rotation_z(angle):
    """Rotate about z-axis"""
    return np.array(
        [[np.cos(angle), np.sin(angle), 0], [-np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
    )


def test_initial_orientation(eros):
    """Tests that initial orientation has the correct pole and that extra rotation rotates the initial orientation"""
    orientation = eros.initial_orientation("ECLIPJ2000")

    longitude, latitude = eros.pole_longitude, eros.pole_latitude
    pole = [
        np.cos(latitude) * np.cos(longitude),
        np.cos(latitude) * np.sin(longitude),
        np.sin(latitude),
    ]
    np.testing.assert_allclose(orientation[2], pole, atol=1e-12)

    np.testing.assert_allclose(
        eros.initial_orientation("ECLIPJ2000", extra_rotation=0.3),
        _passive_rotation_z(0.3) @ orientation,
        atol=1e-12,
    )


@pytest.mark.parametrize("ellipsoid", [False, True], ids=["polyhedron", "ellipsoid"])
def test_rotation_model_settings(ellipsoid):
    """Tests wiring of DamitModel rotation parameters to rotation model settings for the polyhedron and ellipsoid"""
    angle = 0.3
    model = _box_model(3.0, 2.0, 1.0, angle=angle)
    if ellipsoid:
        create_rotation_settings = environment_setup.rotation_model.damit_best_fitting_ellipsoid
        extra_rotation = angle
    else:
        create_rotation_settings = environment_setup.rotation_model.damit_polyhedron
        extra_rotation = 0.0

    body_settings = environment_setup.BodyListSettings("SSB", "ECLIPJ2000")
    body_settings.add_empty_settings("Box")
    body_settings.get("Box").rotation_model_settings = create_rotation_settings(
        model, "ECLIPJ2000", "Box_fixed"
    )
    rotation = environment_setup.create_system_of_bodies(body_settings).get("Box").rotation_model

    initial = model.initial_orientation("ECLIPJ2000", extra_rotation)
    epoch = model.reference_epoch
    np.testing.assert_allclose(rotation.inertial_to_body_fixed_rotation(epoch), initial, atol=1e-12)

    dt = 600.0
    turned = _passive_rotation_z(2 * np.pi * dt / model.period) @ initial
    np.testing.assert_allclose(
        rotation.inertial_to_body_fixed_rotation(epoch + dt), turned, atol=1e-9
    )


def test_rotation_invalid_base_frame(eros):
    """Only J2000 and ECLIPJ2000 are supported base frames"""
    with pytest.raises(ValueError):
        environment_setup.rotation_model.damit_polyhedron(eros, "IAU_EARTH", "Eros_fixed")


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
    """Replace the DAMIT tables by small local tables with several models of Eros, and model downloads by the local
    Eros files. All models load the same Eros files; only the table data (quality flag, size) differs between them.

    Returns the list of model IDs that were downloaded."""
    asteroids = tmp_path / "asteroids.csv"
    asteroids.write_text("id,number,name\n1,433,Eros\n")
    asteroid_models = tmp_path / "asteroid_models.csv"
    asteroid_models.write_text(
        "id,asteroid_id,lsm,calibrated_size,equiv_diameter,quality_flag\n"
        "10,1,LS,0,16.84,2\n"  # low quality
        "11,1,LS,0,20.0,4\n"  # best quality
        "12,1,LS,1,,1\n"  # calibrated size
        "13,1,LS,0,,3\n"  # no size information
    )  # Text column as in the real table, otherwise pandas turns the selected model ID into a float
    monkeypatch.setattr(damit_module, "_DAMIT_ASTEROIDS", asteroids)
    monkeypatch.setattr(damit_module, "_DAMIT_ASTEROID_MODELS", asteroid_models)

    downloaded = []

    def download_model(model_id, force=False):
        downloaded.append(model_id)
        return EROS_SHAPE, EROS_SPIN

    monkeypatch.setattr(damit_module, "_download_model", download_model)
    return downloaded


def test_query_selection_offline(offline_damit, eros):
    """Model selection and sizing from the DAMIT tables"""
    # Highest quality flag is selected, and scaled to its catalog diameter (km -> m)
    best = DamitModel.query(mpc_number=433)
    assert offline_damit == [11]
    assert _model_volume(best) == pytest.approx(_sphere_volume(20_000.0))

    # Given diameter overrides the catalog diameter
    assert _model_volume(DamitModel.query(damit_model_id=10, diameter=1_000.0)) == pytest.approx(
        _sphere_volume(1_000.0)
    )

    # Calibrated model is not rescaled
    assert _model_volume(DamitModel.query(damit_model_id=12)) == pytest.approx(_model_volume(eros))

    # No size information available
    with pytest.warns(UserWarning):
        DamitModel.query(damit_model_id=13)

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

    assert _model_volume(model_by_mpc) == pytest.approx(_model_volume(model_by_id))
    assert model_by_mpc.pole_longitude == pytest.approx(model_by_id.pole_longitude)
    assert model_by_mpc.pole_latitude == pytest.approx(model_by_id.pole_latitude)
    assert model_by_mpc.period == pytest.approx(model_by_id.period)


def test_model_plotting(eros):
    """Tests that plotting the model does not raise and returns Axes"""
    assert isinstance(eros.plot_shape_model(show=False), plt.Axes)
