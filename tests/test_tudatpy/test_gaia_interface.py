"""Tests for the Gaia tracking-data and environment-data interfaces."""

from pathlib import Path
from unittest import mock
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from scipy.linalg import block_diag
from tudatpy.data_input.tracking_data.gaia import (
    GaiaAstrometry,
    generate_astrometry_parquet,
)
from tudatpy.data_input.environment_data.gaia import (
    gaia_object_catalog,
    get_kepler_covariance_from_gaia_archive,
    get_state_covariance_from_gaia_archive,
    get_state_from_gaia_archive,
)
from tudatpy.astro.time_representation import (
    date_time_components_to_epoch,
    julian_day_to_seconds_since_epoch,
)
from tudatpy.astro.element_conversion import j2000_to_eclipj2000
from tudatpy.dynamics.environment_setup import get_default_body_settings, create_system_of_bodies
from tudatpy.estimation import observations
from tudatpy.interface import spice
from tudatpy.math.interpolators import (
    lagrange_interpolation,
    create_one_dimensional_vector_interpolator,
    BoundaryInterpolationType,
    LagrangeInterpolatorBoundaryHandling,
)

TEST_ASTEROID_MPC = (673, 779)  # Asteroids used in subsequent tests (673 Edda and 779 Nina)
TEST_DIR = Path(__file__).parent
ASTROMETRY_ARCHIVE_PATH = (
    TEST_DIR / "gaia_astro_archive_for_tests.parquet"
)  # Astrometry for 673 and 779 only
SOURCE_ARCHIVE_PATH = (
    TEST_DIR / "gaia_source_archive_for_tests.parquet"
)  # Orbit data for 673 and 779 only


##############################
# GaiaAstrometry tests
##############################


@pytest.fixture(scope="module")
def spice_kernels():
    """Load spice kernels, necessary for setting up a SystemOfBodies"""
    spice.load_standard_kernels()


@pytest.fixture(scope="module")
def _gaia_astrometry():
    """Create a GaiaAstrometry object for asteroids 673 and 779 from the local test archive.
    Runs once per module"""
    gaia_astrometry = GaiaAstrometry.load_from_local_archive(
        ASTROMETRY_ARCHIVE_PATH, TEST_ASTEROID_MPC
    )
    gaia_astrometry.apply_filters()
    return gaia_astrometry


@pytest.fixture
def gaia_astrometry(_gaia_astrometry):
    """Copy of the GaiaAstrometry object to avoid mutation in subsequent tests"""
    return _gaia_astrometry.copy()


@pytest.fixture
def astrometry_table(gaia_astrometry):
    """Shortcut for gaia_astrometry.table"""
    return gaia_astrometry.table


@pytest.fixture
def observation_dataset(gaia_astrometry, spice_kernels):
    """ObservationDataset of gaia_astrometry with observations for 673 and 779"""
    body_settings = get_default_body_settings(["Sun"], "SSB", "J2000")
    body_settings.add_empty_settings("Gaia")
    body_settings.get("Gaia").ephemeris_settings = gaia_astrometry.get_gaia_ephemeris_settings(
        geocentric=False
    )
    bodies = create_system_of_bodies(body_settings)
    return gaia_astrometry.to_observation_dataset(bodies)


@pytest.mark.remote_data
def test_local_astroquery_consistency(astrometry_table):
    """Test whether data obtained from the local archive and astroquery are consistent"""
    astrometry_from_astroquery = GaiaAstrometry.load_from_astroquery(TEST_ASTEROID_MPC)
    astrometry_from_astroquery.apply_filters()
    astrometry_table_from_astroquery = astrometry_from_astroquery.table

    pdt.assert_frame_equal(astrometry_table, astrometry_table_from_astroquery, check_dtype=False)


@pytest.mark.remote_data
def test_astrometry_retrieval_no_data():
    """Test that no data available raises error"""
    asteroid_no_data = 250000

    # Online variant
    with pytest.raises(RuntimeError):
        GaiaAstrometry.load_from_astroquery([asteroid_no_data])
    # Local variant
    with pytest.raises(RuntimeError):
        GaiaAstrometry.load_from_local_archive(ASTROMETRY_ARCHIVE_PATH, [asteroid_no_data])


def test_mpc_number_zero_is_rejected():
    """Natural satellites cannot be distinguished using their shared MPC number zero."""
    with pytest.raises(ValueError, match="MPC number 0"):
        GaiaAstrometry(pd.DataFrame({"number_mp": [0]}))
    with pytest.raises(ValueError, match="MPC number 0"):
        get_state_from_gaia_archive(0, SOURCE_ARCHIVE_PATH)


def test_generate_astrometry_parquet_accepts_compressed_csv(tmp_path):
    """The parquet generator accepts the compressed files distributed by the Gaia archive."""
    archive_dir = tmp_path / "archive"
    archive_dir.mkdir()

    chunks = []
    for file_number in range(20):
        chunk = pd.DataFrame(
            {
                "number_mp": [20 - file_number],
                "value": [file_number],
            }
        )
        chunk.to_csv(
            archive_dir / f"SsoObservation_{file_number:02d}.csv.gz",
            index=False,
            compression="gzip",
        )
        chunks.append(chunk)

    generate_astrometry_parquet(archive_dir, tmp_path)

    expected = pd.concat(chunks, ignore_index=True).sort_values("number_mp", ignore_index=True)
    result = pd.read_parquet(tmp_path / "gaia_astrometry_archive.parquet")
    pdt.assert_frame_equal(result, expected)


def test_observation_table_angles(astrometry_table):
    """Test units and range of angles in observation table"""
    # Check if angles are in radians and in correct range
    ra, dec, pa = astrometry_table[["ra", "dec", "position_angle_scan"]].to_numpy().T

    assert np.all((ra < np.pi) & (ra > -np.pi))
    assert np.all((dec < np.pi / 2) & (dec > -np.pi / 2))
    assert np.all((pa < 2 * np.pi) & (pa > 0))


def test_table_ordering(astrometry_table):
    """Test ordering of table entries (first by MPC, then by epoch)"""
    # MPC ordering
    assert astrometry_table["number_mp"].is_monotonic_increasing

    # Epoch ordering
    for mpc in TEST_ASTEROID_MPC:
        assert astrometry_table.loc[
            astrometry_table["number_mp"] == mpc, "epoch"
        ].is_monotonic_increasing


def test_observation_table_uncertainty(astrometry_table):
    """Check that uncertainty and correlation are in the correct ranges"""
    unc = astrometry_table[
        ["ra_error_random", "dec_error_random", "ra_error_systematic", "dec_error_systematic"]
    ].to_numpy()
    corr = astrometry_table[
        ["ra_dec_correlation_random", "ra_dec_correlation_systematic"]
    ].to_numpy()

    assert np.all(unc >= 0)
    assert np.all((corr >= -1) & (corr <= 1))


def test_observation_table_epochs(astrometry_table):
    """Test if the epochs are in the range of the epochs for Gaia FPR"""
    earliest_epoch = date_time_components_to_epoch(
        year=2014, month=7, day=26, hour=0, minute=0, seconds=0
    )
    final_epoch = date_time_components_to_epoch(
        year=2020, month=1, day=20, hour=0, minute=0, seconds=0
    )
    epochs = astrometry_table["epoch"].to_numpy()

    assert np.all((epochs >= earliest_epoch) & (epochs <= final_epoch))


# NOTE: Gaia ephemeris tests are loose because possibly different planetary ephemerides models were used for
# the JPL/Gaia orbit estimations
def test_gaia_barycentric_state_with_jpl_ephemeris(astrometry_table):
    """Test that the barycentric states of Gaia in the table are consistent with those from JPL Horizons"""
    columns = ["x_gaia", "y_gaia", "z_gaia", "vx_gaia", "vy_gaia", "vz_gaia"]
    states_from_archive = astrometry_table.sort_values(by="epoch")[columns].to_numpy()
    states_from_jpl = np.loadtxt(TEST_DIR / "gaia_ephemeris_barycentric.txt")[:, 1:]

    np.testing.assert_allclose(states_from_archive, states_from_jpl, rtol=1e-4, atol=1e-4)


def test_gaia_geocentric_state_with_jpl_ephemeris(astrometry_table):
    """Test that the geocentric states of Gaia in the table are consistent with those from JPL Horizons"""
    columns = [
        "x_gaia_geocentric",
        "y_gaia_geocentric",
        "z_gaia_geocentric",
        "vx_gaia_geocentric",
        "vy_gaia_geocentric",
        "vz_gaia_geocentric",
    ]
    states_from_archive = astrometry_table.sort_values(by="epoch")[columns].to_numpy()
    states_from_jpl = np.loadtxt(TEST_DIR / "gaia_ephemeris_geocentric.txt")[:, 1:]

    np.testing.assert_allclose(states_from_archive, states_from_jpl, rtol=1e-4, atol=1e-4)


def test_observation_table_epoch_filter(gaia_astrometry):
    """Test filtering of epochs with float bounds"""
    filter_start = date_time_components_to_epoch(
        year=2017, month=1, day=1, hour=0, minute=0, seconds=0
    )
    filter_end = date_time_components_to_epoch(
        year=2018, month=1, day=1, hour=0, minute=0, seconds=0
    )

    # Check first that observations outside of this time range exist
    assert gaia_astrometry.table["epoch"].min() < filter_start
    assert gaia_astrometry.table["epoch"].max() > filter_end

    gaia_astrometry.apply_filters(epoch_start=filter_start, epoch_end=filter_end)

    assert gaia_astrometry.table["epoch"].min() >= filter_start
    assert gaia_astrometry.table["epoch"].max() <= filter_end


def test_apply_corrections_photocenter(gaia_astrometry):
    """Tests that apply_corrections correctly modifies observation by calculated photocenter offset"""
    get_obs = lambda table, mpc: table.loc[table["number_mp"] == mpc, ["ra", "dec"]].to_numpy()

    astrometry_table = gaia_astrometry.table
    observations_uncorr = {mpc: get_obs(astrometry_table, mpc) for mpc in TEST_ASTEROID_MPC}

    # Make photocenter offset function return a fixed offset
    offset = 1e-9  # 1e-9 radians in RA and DEC
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), offset)

    # Corrections are only applied to asteroids included in the dimensions mapping
    body_dimensions = {TEST_ASTEROID_MPC[0]: 1e3}
    with mock.patch(
        "tudatpy.data_input.tracking_data.gaia.gaia.photocenter_correction_angular_observations",
        side_effect=fake_correction,
    ):
        gaia_astrometry.apply_corrections(
            bodies=None, photocenter_body_dimensions=body_dimensions, light_deflection_bodies=None
        )
    astrometry_table = gaia_astrometry.table

    for mpc in TEST_ASTEROID_MPC:
        observations_corr = get_obs(astrometry_table, mpc)
        expected_offset = offset if mpc in body_dimensions else 0.0
        np.testing.assert_allclose(
            observations_corr - observations_uncorr[mpc],
            np.full_like(observations_corr, expected_offset),
            rtol=1e-7,
            atol=0,
        )


def test_apply_corrections_light_deflection(gaia_astrometry):
    """Test that apply_corrections correctly applies a light deflection offset to the observations"""
    get_obs = lambda table, mpc: table.loc[table["number_mp"] == mpc, ["ra", "dec"]].to_numpy()

    astrometry_table = gaia_astrometry.table
    observations_uncorr = {mpc: get_obs(astrometry_table, mpc) for mpc in TEST_ASTEROID_MPC}

    # Make light deflection function return a fixed offset
    offset = 1e-9  # 1e-9 radians in RA and DEC
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), offset)

    # Corrections are applied to all loaded asteroids
    with mock.patch(
        "tudatpy.data_input.tracking_data.gaia.gaia.light_deflection_correction_angular_observations",
        side_effect=fake_correction,
    ):
        gaia_astrometry.apply_corrections(bodies=None, light_deflection_bodies=["Sun"])
    astrometry_table = gaia_astrometry.table

    for mpc in TEST_ASTEROID_MPC:
        observations_corr = get_obs(astrometry_table, mpc)
        np.testing.assert_allclose(
            observations_corr - observations_uncorr[mpc],
            np.full_like(observations_corr, offset),
            rtol=1e-7,
            atol=0,
        )


def test_apply_corrections_is_atomic(gaia_astrometry):
    """A failed correction must leave the astrometry table unchanged."""
    table_before_correction = gaia_astrometry.table

    def fail_for_second_asteroid(observations, body_name, **kwargs):
        if body_name == str(TEST_ASTEROID_MPC[1]):
            raise ValueError("Missing ephemeris")
        return np.full((len(observations), 2), 1e-9)

    with mock.patch(
        "tudatpy.data_input.tracking_data.gaia.gaia.light_deflection_correction_angular_observations",
        side_effect=fail_for_second_asteroid,
    ):
        with pytest.raises(ValueError, match="Missing ephemeris"):
            gaia_astrometry.apply_corrections(bodies=None, light_deflection_bodies=["Sun"])

    pdt.assert_frame_equal(gaia_astrometry.table, table_before_correction)


def test_apply_corrections_twice_raises_error(gaia_astrometry):
    """Applying corrections twice on the same instance must raise an error"""
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), 1e-9)

    with mock.patch(
        "tudatpy.data_input.tracking_data.gaia.gaia.light_deflection_correction_angular_observations",
        side_effect=fake_correction,
    ):
        gaia_astrometry.apply_corrections(bodies=None, light_deflection_bodies=["Sun"])

        with pytest.raises(RuntimeError):
            gaia_astrometry.apply_corrections(bodies=None, light_deflection_bodies=["Sun"])


def test_get_gaia_ephemeris_settings_geocentric(gaia_astrometry, spice_kernels):
    """Test if states in catalog and those retrieved from ephemeris match (geocentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_astrometry.get_gaia_ephemeris_settings(geocentric=True)
    body_settings = get_default_body_settings(["Sun", "Earth"], "SSB", "J2000")
    body_settings.add_empty_settings("Gaia")
    body_settings.get("Gaia").ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get("Gaia").ephemeris

    # Compare state vectors from Tudat and table
    astrometry_table = gaia_astrometry.table
    states_from_tudat = np.array(
        [gaia_ephemeris.cartesian_state(epoch) for epoch in astrometry_table.epoch]
    )
    states_from_table = astrometry_table[
        [
            "x_gaia_geocentric",
            "y_gaia_geocentric",
            "z_gaia_geocentric",
            "vx_gaia_geocentric",
            "vy_gaia_geocentric",
            "vz_gaia_geocentric",
        ]
    ].to_numpy()

    np.testing.assert_array_equal(states_from_table, states_from_tudat)


def test_get_gaia_ephemeris_settings_barycentric(gaia_astrometry, spice_kernels):
    """Test if states in catalog and those retrieved from ephemeris match (barycentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_astrometry.get_gaia_ephemeris_settings(geocentric=False)
    body_settings = get_default_body_settings(["Sun", "Earth"], "SSB", "J2000")
    body_settings.add_empty_settings("Gaia")
    body_settings.get("Gaia").ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get("Gaia").ephemeris

    # Compare state vectors from Tudat and table
    astrometry_table = gaia_astrometry.table
    states_from_tudat = np.array(
        [gaia_ephemeris.cartesian_state(epoch) for epoch in astrometry_table.epoch]
    )
    states_from_table = astrometry_table[
        ["x_gaia", "y_gaia", "z_gaia", "vx_gaia", "vy_gaia", "vz_gaia"]
    ].to_numpy()

    np.testing.assert_array_equal(states_from_table, states_from_tudat)


def test_summary_does_not_raise(gaia_astrometry, capsys):
    """print_summary() should run without error and print a summary header"""
    gaia_astrometry.print_summary()

    assert "SUMMARY" in capsys.readouterr().out


def test_to_observation_dataset_without_gaia_raises_error(gaia_astrometry, spice_kernels):
    """to_observation_dataset must raise an error if Gaia is not loaded in bodies"""
    body_settings = get_default_body_settings(["Sun"], "SSB", "J2000")
    bodies = create_system_of_bodies(body_settings)

    with pytest.raises(ValueError):
        gaia_astrometry.to_observation_dataset(bodies)


def test_weight_matrix_symmetry(observation_dataset):
    """Test that the weight matrix is symmetric"""
    weight_matrix = observation_dataset.get_weight_matrix().toarray()

    np.testing.assert_allclose(weight_matrix, weight_matrix.T, rtol=1e-10, atol=0)


def test_covariance_matrix_variance_consistency(gaia_astrometry):
    """Test that the covariance matrix diagonal matches the per-observation variances in the
    astrometry table, ordered by epoch (RA then DEC per observation)."""
    # Concatenate the per-transit covariance blocks into one matrix over all observations
    covariance = block_diag(*gaia_astrometry._get_observation_covariance(TEST_ASTEROID_MPC[0]))
    variance_from_matrix = np.diag(covariance)

    table = gaia_astrometry.table
    table = table[table["number_mp"] == TEST_ASTEROID_MPC[0]].reset_index(drop=True)

    # By design, the systematic error is constant per transit: use the first entry of each transit_id
    systematic = table.groupby("transit_id", sort=False)[
        ["ra_error_systematic", "dec_error_systematic"]
    ].transform("first")

    var_ra = table["ra_error_random"] ** 2 + systematic["ra_error_systematic"] ** 2
    var_dec = table["dec_error_random"] ** 2 + systematic["dec_error_systematic"] ** 2
    variance_from_table = np.column_stack([var_ra, var_dec]).ravel()

    np.testing.assert_array_equal(variance_from_matrix, variance_from_table)


def test_covariance_weight_matrix_consistency(observation_dataset, gaia_astrometry):
    """Test that the covariance matrix computed from the astrometry table is consistent with weight matrix
    passed to observation_dataset"""
    condition = observations.observation_query.transmitter == observations.LinkEndId(
        str(TEST_ASTEROID_MPC[0]), ""
    )
    weight_matrix = observation_dataset.get_weight_matrix(condition).toarray()
    covariance_from_obscol = np.linalg.inv(weight_matrix)
    covariance_from_table = block_diag(
        *gaia_astrometry._get_observation_covariance(TEST_ASTEROID_MPC[0])
    )

    # NOTE: loose tolerance because of high condition number so much precision is lost during 2 inversions
    np.testing.assert_allclose(covariance_from_obscol, covariance_from_table, rtol=1e-5, atol=0)


def test_covariance_matrix_nonzero_entries(observation_dataset, astrometry_table):
    """Test block diagonal structure by checking number of nonzero entries in covariance matrix"""
    # Count number of nonzero entries in covariance matrix
    condition = observations.observation_query.transmitter == observations.LinkEndId(
        str(TEST_ASTEROID_MPC[0]), ""
    )
    weight_matrix = observation_dataset.get_weight_matrix(condition).toarray()
    covariance = np.linalg.inv(weight_matrix)
    nonzeros_covariance_matrix = np.count_nonzero(covariance)

    # Calculate expected number of nonzero entries for a block-diagonal structure
    table = astrometry_table[astrometry_table["number_mp"] == TEST_ASTEROID_MPC[0]]
    transit_ids = pd.unique(table["transit_id"])
    transit_lengths = [
        2 * np.count_nonzero(table["transit_id"] == id) for id in transit_ids
    ]  # times 2 because of RA, DEC
    nonzeros_expected = sum([n**2 for n in transit_lengths])

    assert nonzeros_covariance_matrix == nonzeros_expected


def test_observation_consistency(observation_dataset, astrometry_table):
    """Test if observations in the table are consistent with those in ObservationDataset"""
    condition = observations.observation_query.transmitter == observations.LinkEndId(
        str(TEST_ASTEROID_MPC[0]), ""
    )
    observations_from_dataset = np.ravel(observation_dataset.get_observations(condition))
    epochs_from_dataset = np.array(
        [epoch.to_float() for epoch in observation_dataset.get_times(condition)]
    )

    # Get observations and epochs from table
    table = astrometry_table[astrometry_table["number_mp"] == TEST_ASTEROID_MPC[0]]
    epochs_from_table = table.epoch.to_numpy()
    observations_from_table = np.ravel(table[["ra", "dec"]])

    np.testing.assert_array_equal(epochs_from_dataset, epochs_from_table)
    np.testing.assert_array_equal(observations_from_dataset, observations_from_table)


##############################
# Gaia asteroid archive tests
##############################

# Reference J2000 Heliocentric state vectors retrieved from JPL Horizons at the same epoch as Gaia catalog
STATE_EDDA_JPL = [
    4.128085144236671e08,
    7.098887570553194e07,
    4.413226270395131e07,
    -3.659135280911876e00,
    1.621253557540478e01,
    6.231297081940196e00,
]  # KM, KM/s
STATE_NINA_JPL = [
    2.143720589223085e08,
    2.310721906203461e08,
    1.782520574512776e08,
    -1.254695764130771e01,
    1.502661150744586e01,
    4.097474113261074e00,
]


@pytest.fixture(scope="module")
def gaia_asteroid_catalog():
    """Load the test asteroid catalog."""
    return gaia_object_catalog(SOURCE_ARCHIVE_PATH)


@pytest.mark.remote_data
def test_asteroid_local_astroquery_consistency():
    """Test whether asteroid data obtained from the local archive and astroquery are consistent"""
    mpc_number = TEST_ASTEROID_MPC[0]

    epoch_local, state_local = get_state_from_gaia_archive(mpc_number, SOURCE_ARCHIVE_PATH)
    epoch_astroquery, state_astroquery = get_state_from_gaia_archive(mpc_number)
    assert epoch_local == epoch_astroquery
    np.testing.assert_allclose(state_local, state_astroquery)

    for get_covariance in [
        get_state_covariance_from_gaia_archive,
        get_kepler_covariance_from_gaia_archive,
    ]:
        epoch_local, covariance_local = get_covariance(mpc_number, SOURCE_ARCHIVE_PATH)
        epoch_astroquery, covariance_astroquery = get_covariance(mpc_number)
        assert epoch_local == epoch_astroquery
        np.testing.assert_allclose(covariance_local, covariance_astroquery)


def test_asteroid_data_retrieval_no_data():
    """Test that no local data returned raises an error."""
    asteroid_no_data = 250000

    for get_data in [
        get_state_from_gaia_archive,
        get_state_covariance_from_gaia_archive,
        get_kepler_covariance_from_gaia_archive,
    ]:
        with pytest.raises(LookupError):
            get_data(asteroid_no_data, SOURCE_ARCHIVE_PATH)


@pytest.mark.remote_data
def test_asteroid_astroquery_no_data():
    """Test that no online data returned raises an error."""
    with pytest.raises(LookupError):
        get_state_from_gaia_archive(250000)


def test_asteroid_epochs(gaia_asteroid_catalog):
    """Test that the state vector epochs fall within the Gaia FPR observation period"""
    earliest_epoch = date_time_components_to_epoch(
        year=2014, month=7, day=26, hour=0, minute=0, seconds=0
    )
    final_epoch = date_time_components_to_epoch(
        year=2020, month=1, day=20, hour=0, minute=0, seconds=0
    )
    epochs = gaia_asteroid_catalog["epoch_state_vector"].to_numpy()

    assert np.all((epochs >= earliest_epoch) & (epochs <= final_epoch))


def test_asteroid_state_vector_with_jpl_horizons():
    """Test match between Gaia state vectors and that of JPL Horizons"""
    ref_states = [
        np.array(STATE_EDDA_JPL) * 1e3,
        np.array(STATE_NINA_JPL) * 1e3,
    ]  # JPL Horizons vectors
    for mpc_number, ref_state in zip(TEST_ASTEROID_MPC, ref_states):
        # Heliocentric state vectors:
        epoch, state = get_state_from_gaia_archive(mpc_number, SOURCE_ARCHIVE_PATH)

        np.testing.assert_allclose(state, ref_state, rtol=1e-7, atol=1e-7)


def test_asteroid_state_frame_conversion(spice_kernels):
    """Test state origin and orientation conversion."""
    mpc_number = TEST_ASTEROID_MPC[0]
    epoch, state = get_state_from_gaia_archive(mpc_number, SOURCE_ARCHIVE_PATH)

    ecliptic_epoch, ecliptic_state = get_state_from_gaia_archive(
        mpc_number, SOURCE_ARCHIVE_PATH, frame_orientation="ECLIPJ2000"
    )
    rotation = j2000_to_eclipj2000()
    expected_ecliptic_state = np.concatenate((rotation @ state[:3], rotation @ state[3:]))
    assert ecliptic_epoch == epoch
    np.testing.assert_array_equal(ecliptic_state, expected_ecliptic_state)

    earth_epoch, earth_state = get_state_from_gaia_archive(
        mpc_number, SOURCE_ARCHIVE_PATH, frame_origin="Earth"
    )
    sun_from_earth = spice.get_body_cartesian_state_at_epoch(
        target_body_name="Sun",
        observer_body_name="Earth",
        reference_frame_name="J2000",
        aberration_corrections="NONE",
        ephemeris_time=epoch,
    )
    assert earth_epoch == epoch
    np.testing.assert_array_equal(earth_state, state + sun_from_earth)


def test_asteroid_state_invalid_orientation():
    """Test that unsupported state orientations raise an error."""
    with pytest.raises(ValueError, match="frame_orientation"):
        get_state_from_gaia_archive(
            TEST_ASTEROID_MPC[0], SOURCE_ARCHIVE_PATH, frame_orientation="invalid"
        )


def test_asteroid_catalog_columns(gaia_asteroid_catalog):
    """Derived orbit elements and orbit class should not be added to the catalog."""
    assert "orbital_elements" not in gaia_asteroid_catalog
    assert "orbit_class" not in gaia_asteroid_catalog
    assert gaia_asteroid_catalog.index.name == "number_mp"
    assert gaia_asteroid_catalog["number_mp"].is_monotonic_increasing


def test_asteroid_covariance_shape_and_symmetry(gaia_asteroid_catalog):
    """Covariance matrices must be reconstructed from the raw upper triangle as symmetric 6x6 matrices"""
    for _, asteroid_data in gaia_asteroid_catalog.iterrows():
        for column in ["orbital_elements_var_covar_matrix", "h_state_vector_var_covar_matrix"]:
            covariance = asteroid_data[column]

            assert covariance.shape == (6, 6)
            np.testing.assert_allclose(covariance, covariance.T, rtol=1e-12, atol=0)


def test_asteroid_covariance_values(gaia_asteroid_catalog):
    """Sanity checks on the covariance matrices: variances must be positive, and the semi-major axis
    uncertainty must be in a range plausible for Gaia orbit solutions (order of meters, not AU)"""
    for _, asteroid_data in gaia_asteroid_catalog.iterrows():
        for column in ["orbital_elements_var_covar_matrix", "h_state_vector_var_covar_matrix"]:
            assert np.all(np.diag(asteroid_data[column]) > 0)

        sma_uncertainty = np.sqrt(asteroid_data["orbital_elements_var_covar_matrix"][0, 0])
        assert 1e-2 < sma_uncertainty < 1e6  # In meters; catches missing or double unit scaling


def test_asteroid_covariance_functions(gaia_asteroid_catalog):
    """Covariance functions return the state epoch and the requested covariance."""
    mpc_number = TEST_ASTEROID_MPC[0]
    object_row = gaia_asteroid_catalog.loc[mpc_number]

    epoch, covariance = get_state_covariance_from_gaia_archive(mpc_number, SOURCE_ARCHIVE_PATH)
    assert epoch == object_row["epoch_state_vector"]
    np.testing.assert_array_equal(covariance, object_row["h_state_vector_var_covar_matrix"])

    epoch, covariance = get_kepler_covariance_from_gaia_archive(mpc_number, SOURCE_ARCHIVE_PATH)
    assert epoch == object_row["epoch_state_vector"]
    np.testing.assert_array_equal(covariance, object_row["orbital_elements_var_covar_matrix"])
