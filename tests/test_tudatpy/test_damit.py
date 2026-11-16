from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

import tudatpy.data_input.environment_data.damit.damit as damit_module
from tudatpy.astro.polyhedron_utilities import volume
from tudatpy.data_input.environment_data.damit import DamitModel
from tudatpy.dynamics import environment_setup
from tudatpy.interface import spice

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
def damit_dir(data_dir) -> Path:
    """Directory holding the DAMIT spin and shape files of (433) Eros"""
    return data_dir / "damit"


@pytest.fixture(scope="module")
def eros_spin(damit_dir) -> Path:
    return damit_dir / "eros_spin.txt"


@pytest.fixture(scope="module")
def eros_shape(damit_dir) -> Path:
    return damit_dir / "eros_shape.obj"


@pytest.fixture(scope="module")
def eros(eros_spin, eros_shape):
    """Shape model of (433) Eros with DAMIT Model ID 3083"""
    return DamitModel.from_file(eros_spin, eros_shape)


def test_from_file_parses_spin_and_shape(eros):
    """Test parsing of shape and spin files and conversion to correct units"""
    assert eros._vertices.shape == (856, 3)
    assert eros._facets.shape == (1708, 3)
    assert eros._facets.min() == 0
    assert eros._facets.max() == 855

    assert eros._pole_longitude == pytest.approx(np.deg2rad(17))
    assert eros._pole_latitude == pytest.approx(np.deg2rad(11))
    assert eros._period == pytest.approx(5.27025528 * 3600)
    assert eros._reference_epoch == pytest.approx(0.0)  # JD 2451545.0 is J2000


def test_volume_scaling(eros, eros_spin, eros_shape):
    """Test that shape model is correctly scaled and that diameter input scales volume to that diameter"""
    # Volume of Eros of 2535 +/- 20 km^3 from the NEAR MSI shape model (Thomas et al., 2002, Icarus 155(1), 18-37,
    # doi:10.1006/icar.2001.6755, abstract). The coarse 1708-plate version of this shape model used here (DAMIT model
    # 3083, PDS plate model eros_2001012_001708) has a 1.7 % smaller volume.
    assert _model_volume(eros) == pytest.approx(2535e9, rel=0.02)

    diameter = 16_840.0
    scaled = DamitModel.from_file(eros_spin, eros_shape, diameter=diameter)
    assert _model_volume(scaled) == pytest.approx(_sphere_volume(diameter))


def test_non_obj_shape_file_raises(tmp_path, eros_spin):
    """A shape file not in OBJ format is rejected"""
    shape_file = tmp_path / "shape.txt"
    shape_file.write_text("4 4\n0 0 0\n")
    with pytest.raises(ValueError):
        DamitModel.from_file(eros_spin, shape_file)


def test_ellipsoid_semi_axes_box():
    """A box with half-sides (x, y, z) has an equal-inertia ellipsoid with semi-axes sqrt(5/3) * (x, y, z)"""
    half_sides = np.array([3.0, 2.0, 1.0])
    model = _box_model(*half_sides)

    semi_axes = model.best_fitting_ellipsoid()

    np.testing.assert_allclose(semi_axes, np.sqrt(5 / 3) * half_sides, rtol=1e-10)


def test_shape_settings(eros):
    """Check wiring of DamitModel to shape settings creation"""
    polyhedron = environment_setup.shape.damit_polyhedron(eros)
    np.testing.assert_array_equal(polyhedron.vertices_coordinates, eros.vertices)
    np.testing.assert_array_equal(polyhedron.vertices_defining_each_facet, eros.facets)

    ellipsoid = environment_setup.shape.damit_best_fitting_ellipsoid(eros)
    np.testing.assert_allclose(
        [ellipsoid.semi_axis_a, ellipsoid.semi_axis_b, ellipsoid.semi_axis_c],
        eros.best_fitting_ellipsoid(),
    )


@pytest.mark.parametrize("angle", [0.3, -0.7, 1.2])
def test_principal_axis_frame_rotated_box(angle):
    """A shape rotated around the z-axis is rotated back to its principal axes, and the rotation is added to the
    rotation angle"""
    model = _box_model(3.0, 2.0, 1.0, angle=angle)

    np.testing.assert_allclose(model.vertices, _box_model(3.0, 2.0, 1.0).vertices, atol=1e-12)
    assert model.initial_rotation_angle == pytest.approx(angle, abs=1e-12)


def _passive_rotation_z(angle):
    """Rotate about z-axis"""
    return np.array(
        [[np.cos(angle), np.sin(angle), 0], [-np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
    )


def test_initial_orientation(eros, eros_shape):
    """The body-fixed vertices, rotated to ECLIPJ2000 at the reference epoch, are the vertices of the DAMIT shape model
    rotated with the DAMIT convention r_ecl = R_z(lambda) R_y(90 deg - beta) R_z(phi_0) r_ast"""
    # Pole longitude lambda, pole latitude beta and initial rotation angle phi_0 (deg) from eros_spin.txt
    pole_longitude, pole_latitude, initial_rotation_angle = np.deg2rad([17.0, 11.0, 32.64])
    colatitude = np.pi / 2 - pole_latitude
    rotation_y = np.array(
        [
            [np.cos(colatitude), 0, np.sin(colatitude)],
            [0, 1, 0],
            [-np.sin(colatitude), 0, np.cos(colatitude)],
        ]
    )
    damit_to_ecliptic = (
        _passive_rotation_z(-pole_longitude)
        @ rotation_y
        @ _passive_rotation_z(-initial_rotation_angle)
    )
    damit_vertices, _ = DamitModel._load_shape_model(eros_shape, None)

    np.testing.assert_allclose(
        eros.vertices @ eros.initial_orientation("ECLIPJ2000"),
        damit_vertices @ damit_to_ecliptic.T,
        atol=1e-6,
    )


def test_initial_orientation_j2000(eros):
    """The J2000 initial orientation is the ECLIPJ2000 one combined with the SPICE rotation from J2000 to ECLIPJ2000"""
    np.testing.assert_allclose(
        eros.initial_orientation("J2000"),
        eros.initial_orientation("ECLIPJ2000")
        @ spice.compute_rotation_matrix_between_frames("J2000", "ECLIPJ2000", 0.0),
        atol=1e-12,
    )


def test_rotation_model_settings():
    """Tests wiring of DamitModel rotation parameters to rotation model settings"""
    model = _box_model(3.0, 2.0, 1.0, angle=0.3)

    body_settings = environment_setup.BodyListSettings("SSB", "ECLIPJ2000")
    body_settings.add_empty_settings("Box")
    body_settings.get("Box").rotation_model_settings = environment_setup.rotation_model.damit(
        model, "ECLIPJ2000", "Box_fixed"
    )
    rotation = environment_setup.create_system_of_bodies(body_settings).get("Box").rotation_model

    initial = model.initial_orientation("ECLIPJ2000")
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
        environment_setup.rotation_model.damit(eros, "IAU_EARTH", "Eros_fixed")


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
def offline_damit(tmp_path, eros_spin, eros_shape):
    """Replace the DAMIT tables by small local tables with several models of Eros, and model downloads by the local
    Eros files. All models load the same Eros files; only the table data (quality flag, size) differs between them.

    Returns the list of model IDs that were downloaded."""
    asteroids = tmp_path / "asteroids.csv"
    asteroids.write_text("id,number,name\n1,433,Eros\n4,4179,Toutatis\n")
    asteroid_models = tmp_path / "asteroid_models.csv"
    asteroid_models.write_text(
        "id,asteroid_id,lsm,calibrated_size,equiv_diameter,quality_flag,yorp\n"
        "10,1,LS,0,16.84,2,\n"  # low quality
        "11,1,LS,0,20.0,4,\n"  # best quality
        "12,1,LS,1,,1,\n"  # calibrated size
        "18,4,LS,1,,3,\n"  # two equally good models, highest ID first
        "17,4,LS,1,,3,\n"
    )

    downloaded = []

    def download_model(model_id, force=False):
        downloaded.append(model_id)
        return eros_shape, eros_spin

    with (
        patch.object(damit_module, "_DAMIT_ASTEROIDS", asteroids),
        patch.object(damit_module, "_DAMIT_ASTEROID_MODELS", asteroid_models),
        patch.object(damit_module, "_download_model", download_model),
    ):
        yield downloaded


def test_query_selection_offline(offline_damit, eros):
    """Model selection and sizing from the DAMIT tables"""
    # Highest quality flag is selected, and scaled to its catalog diameter (km -> m)
    best = DamitModel.query(mpc_number=433)
    assert offline_damit == [11]
    assert _model_volume(best) == pytest.approx(_sphere_volume(20_000.0))
    assert best.number_of_equally_good_models == 1

    # Equally good models: lowest model ID is selected, whatever the table order
    equally_good = DamitModel.query(mpc_number=4179)
    assert offline_damit[-1] == 17
    assert equally_good.number_of_equally_good_models == 2
    assert DamitModel.query(damit_model_id=18).number_of_equally_good_models is None

    # Given diameter overrides the catalog diameter
    assert _model_volume(DamitModel.query(damit_model_id=10, diameter=1_000.0)) == pytest.approx(
        _sphere_volume(1_000.0)
    )

    # Calibrated model is not rescaled
    assert _model_volume(DamitModel.query(damit_model_id=12)) == pytest.approx(_model_volume(eros))

    # Unknown asteroid or model
    with pytest.raises(LookupError):
        DamitModel.query(mpc_number=1)
    with pytest.raises(LookupError):
        DamitModel.query(damit_model_id=99)


@pytest.mark.remote_data
def test_query_remote(tmp_path):
    """Querying the online database (into a temporary folder) by MPC number, and by the ID of the model that this
    should select, gives the same model"""
    with (
        patch.object(damit_module, "_DAMIT_ASTEROIDS", tmp_path / "tables" / "asteroids.csv"),
        patch.object(
            damit_module, "_DAMIT_ASTEROID_MODELS", tmp_path / "tables" / "asteroid_models.csv"
        ),
        patch.object(damit_module, "_DAMIT_MODELS_DIR", tmp_path / "models"),
    ):
        models = damit_module.list_damit_models(mpc_number=TEST_ASTEROID)
        best_model_id = int(
            models[models["quality_flag"] == models["quality_flag"].max()].index.min()
        )

        model_by_mpc = DamitModel.query(mpc_number=TEST_ASTEROID)
        model_by_id = DamitModel.query(damit_model_id=best_model_id, force_refresh=True)

    assert _model_volume(model_by_mpc) == pytest.approx(_model_volume(model_by_id))
    assert model_by_mpc.pole_longitude == pytest.approx(model_by_id.pole_longitude)
    assert model_by_mpc.pole_latitude == pytest.approx(model_by_id.pole_latitude)
    assert model_by_mpc.period == pytest.approx(model_by_id.period)
