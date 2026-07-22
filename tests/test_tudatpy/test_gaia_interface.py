"""
Tests for the GaiaAstrometry class in tudatpy.data.gaia.gaia.

The GaiaAsteroids class is not covered by these tests.
"""
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tudatpy.data.gaia.gaia import GaiaAstrometry
from tudatpy.astro.time_representation import (
    date_time_components_to_epoch,
    julian_day_to_seconds_since_epoch,
)
from tudatpy.dynamics.environment_setup import get_default_body_settings, create_system_of_bodies
from tudatpy.estimation.observations.observations_processing import observation_parser
from tudatpy.interface import spice
from tudatpy.math.interpolators import (
    lagrange_interpolation,
    create_one_dimensional_vector_interpolator,
    BoundaryInterpolationType,
    LagrangeInterpolatorBoundaryHandling,
)

# Constants
TEST_ASTEROID_MPC = (433, 779)  # Asteroids used in subsequent tests
TEST_DIR = Path(__file__).parent
LOCAL_ARCHIVE_PATH = TEST_DIR / 'gaia_archive_for_tests.parquet'  # Contains data for 433 and 779 only


# ── HELPERS ───────────────────────────────────────────────────────────────────

def read_jpl_table(path):
    """
    Read and parse JPL Horizons vector table results file.
    """
    with open(path) as f:
        lines = f.readlines()

    # Locate the data block between the $$SOE and $$EOE markers.
    start = next(i for i, line in enumerate(lines) if line.strip() == "$$SOE")
    end = next(i for i, line in enumerate(lines) if line.strip() == "$$EOE")

    rows = []
    for line in lines[start + 1:end]:
        line = line.strip()
        if not line:
            continue
        # Split on commas; columns are JDTDB, Calendar Date, then the 12
        # numeric fields. A trailing comma yields an empty final element.
        fields = [field.strip() for field in line.split(",")]
        jd = float(fields[0])
        # fields[1] is the calendar date string, which we drop.
        values = [float(field) for field in fields[2:14] if field != 'n.a.']
        rows.append([jd] + values)

    data = np.array(rows)

    # Convert time (JDTDB) to TDB seconds since J2000, and states/uncertainty to SI
    data[:, 0] = np.array([julian_day_to_seconds_since_epoch(jd) for jd in data[:, 0]])
    data[:, 1:] *= 1000

    state_history = data[:, :7]
    if data.shape[1] > 7:  # Data contains uncertainty values
        uncertainty = np.column_stack((data[:, 0], data[:, 7:]))

    else:
        uncertainty = None

    return state_history, uncertainty


def interpolate_state_history(state_history: np.ndarray,
                              epochs_output: list | np.ndarray,
                              order: int = 8) -> np.ndarray:
    """Interpolate a state history array (epoch, x, y, z, vx, vy, vz) to the requested epochs
    using Lagrange interpolation."""
    state_history_dict = {epoch: state for epoch, state in zip(state_history[:, 0], state_history[:, 1:])}

    interpolation_settings = lagrange_interpolation(
        order,
        boundary_interpolation=BoundaryInterpolationType.throw_exception_at_boundary,
        lagrange_boundary_handling=LagrangeInterpolatorBoundaryHandling.lagrange_no_boundary_interpolation
    )
    interpolator = create_one_dimensional_vector_interpolator(state_history_dict, interpolation_settings)

    return np.array([interpolator.interpolate(epoch) for epoch in epochs_output])


# ── FIXTURES ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope='session')
def spice_kernels():
    """Load spice kernels, necessary for setting up a SystemOfBodies"""
    spice.load_standard_kernels()


@pytest.fixture(scope='module')
def _gaia_astrometry():
    """Create a GaiaAstrometry object for asteroids 433 and 779 from the local test archive.
    Runs once per module"""
    return GaiaAstrometry.load_from_local_archive(LOCAL_ARCHIVE_PATH, TEST_ASTEROID_MPC)


@pytest.fixture
def gaia_astrometry(_gaia_astrometry):
    """Copy of the GaiaAstrometry object to avoid mutation in subsequent tests"""
    return _gaia_astrometry.copy()


@pytest.fixture
def observation_table(gaia_astrometry):
    """Shortcut for gaia_astrometry.table"""
    return gaia_astrometry.table


@pytest.fixture
def observation_collection(gaia_astrometry, spice_kernels):
    """ObservationCollection of gaia_astrometry constructed with to_observation_collection()"""
    # Gaia must be loaded in bodies with its ephemeris to use to_observation_collection()
    body_settings = get_default_body_settings(['Sun'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = gaia_astrometry.get_gaia_ephemeris(geocentric=False)
    bodies = create_system_of_bodies(body_settings)
    return gaia_astrometry.to_observation_collection(bodies)


# ── TESTS ──────────────────────────────────────────────────────────────────

@pytest.mark.remote_data
def test_local_astroquery_consistency(observation_table):
    """Test whether data obtained from the local archive and astroquery are consistent"""
    astrometry_from_astroquery = GaiaAstrometry.load_from_astroquery(TEST_ASTEROID_MPC)
    observation_table_from_astroquery = astrometry_from_astroquery.table

    pdt.assert_frame_equal(observation_table, observation_table_from_astroquery, check_dtype=False)


@pytest.mark.remote_data
def test_data_retrieval_no_data_astroquery():
    """Test that no data returned raises an error for load_from_astroquery"""
    asteroid_no_data = 250000

    with pytest.raises(LookupError):
        GaiaAstrometry.load_from_astroquery([asteroid_no_data])


def test_data_retrieval_no_data_locally():
    """Test that load_from_local_archive raises an error if no data is found"""
    asteroid_no_data = 250000

    with pytest.raises(LookupError):
        GaiaAstrometry.load_from_local_archive(LOCAL_ARCHIVE_PATH, [asteroid_no_data])


def test_table_for_missing_mpc_raises_error(gaia_astrometry):
    """Error must be raised when requesting a table for an asteroid that is not loaded"""
    asteroid_not_loaded = 250000

    with pytest.raises(ValueError):
        gaia_astrometry.table_for_mpc(asteroid_not_loaded)


def test_observation_table_angles(observation_table):
    """Test units and range of angles in observation table"""
    # Check if angles are in radians and in correct range
    ra, dec, pa = observation_table[['ra', 'dec', 'position_angle_scan']].to_numpy().T

    assert np.all((ra < np.pi) & (ra > -np.pi))
    assert np.all((dec < np.pi / 2) & (dec > -np.pi / 2))
    assert np.all((pa < 2 * np.pi) & (pa > 0))


def test_observation_table_uncertainty(observation_table):
    """Check that uncertainty and correlation are in the correct ranges"""
    unc = observation_table[['ra_error_random', 'dec_error_random', 'ra_error_systematic',
                             'dec_error_systematic']].to_numpy()
    corr = observation_table[['ra_dec_correlation_random', 'ra_dec_correlation_systematic']].to_numpy()

    assert np.all(unc >= 0)
    assert np.all((corr >= -1) & (corr <= 1))


def test_observation_table_epochs(observation_table):
    """Test if the epochs are in the range of the epochs for Gaia FPR"""
    earliest_epoch = date_time_components_to_epoch(year=2014, month=7, day=26, hour=0, minute=0, seconds=0)
    final_epoch = date_time_components_to_epoch(year=2020, month=1, day=20, hour=0, minute=0, seconds=0)
    epochs = observation_table['epoch'].to_numpy()

    assert np.all((epochs >= earliest_epoch) & (epochs <= final_epoch))


def test_gaia_barycentric_state_with_jpl_ephemeris(gaia_astrometry):
    """Test that the barycentric states of Gaia in the table are consistent with those from JPL Horizons"""
    # Retrieve observation table for one asteroid only
    table_test_asteroid = gaia_astrometry.table_for_mpc(TEST_ASTEROID_MPC[0])
    epoch_list = table_test_asteroid['epoch'].to_list()
    states_from_horizons, _ = read_jpl_table(TEST_DIR / 'gaia_ephemeris_jpl_barycentric_12h.txt')
    states_from_horizons = interpolate_state_history(states_from_horizons, epoch_list)  # Interpolate to same times
    states_from_table = table_test_asteroid[['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']]
    states_from_table = states_from_table.to_numpy()

    np.testing.assert_allclose(states_from_horizons, states_from_table, rtol=1e-05, atol=10)
    # We use looser tolerances here because a large deviation can exist due to different adopted planetary ephemerides


def test_gaia_geocentric_state_with_jpl_ephemeris(gaia_astrometry):
    """Test that the geocentric states of Gaia in the table are consistent with those from JPL Horizons"""
    # Retrieve observation table for one asteroid only
    table_test_asteroid = gaia_astrometry.table_for_mpc(TEST_ASTEROID_MPC[0])
    epoch_list = table_test_asteroid['epoch'].to_list()
    states_from_horizons, _ = read_jpl_table(TEST_DIR / 'gaia_ephemeris_jpl_geocentric_12h.txt')
    states_from_horizons = interpolate_state_history(states_from_horizons, epoch_list)  # Interpolate to same times
    states_from_table = table_test_asteroid[['x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                                             'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']]
    states_from_table = states_from_table.to_numpy()

    np.testing.assert_allclose(states_from_horizons, states_from_table, rtol=1e-05, atol=10)
    # We use looser tolerances here because a large deviation can exist due to different adopted planetary ephemerides


def test_observation_table_epoch_filter(gaia_astrometry):
    """Test filtering of epochs with float bounds"""
    filter_start = date_time_components_to_epoch(year=2017, month=1, day=1, hour=0, minute=0, seconds=0)
    filter_end = date_time_components_to_epoch(year=2018, month=1, day=1, hour=0, minute=0, seconds=0)

    # Check first that observations outside of this time range exist
    assert gaia_astrometry.table['epoch'].min() < filter_start
    assert gaia_astrometry.table['epoch'].max() > filter_end

    gaia_astrometry.filter(epoch_start=filter_start, epoch_end=filter_end)

    assert gaia_astrometry.table['epoch'].min() >= filter_start
    assert gaia_astrometry.table['epoch'].max() <= filter_end


def test_correct_observations_photocenter(gaia_astrometry):
    """Tests that correct_observations correctly applies a photocenter offset to the observations"""
    get_obs = lambda table, mpc: table.loc[table['number_mp'] == mpc, ['ra', 'dec']].to_numpy()

    observation_table = gaia_astrometry.table
    observations_uncorr = {mpc: get_obs(observation_table, mpc) for mpc in TEST_ASTEROID_MPC}

    # Make photocenter offset function return a fixed offset
    offset = 1e-9  # 1e-9 radians in RA and DEC
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), offset)

    # Corrections are applied to all loaded asteroids
    diameters = {mpc: 1e3 for mpc in TEST_ASTEROID_MPC}
    with mock.patch('tudatpy.data.gaia.gaia.photocenter_corrections_from_observations',
                    side_effect=fake_correction):
        gaia_astrometry.correct_observations(bodies=None, diameters=diameters,
                                             light_deflection_bodies=[], correct_photocenter=True)
    observation_table = gaia_astrometry.table

    for mpc in TEST_ASTEROID_MPC:
        observations_corr = get_obs(observation_table, mpc)
        np.testing.assert_allclose(observations_corr - observations_uncorr[mpc],
                                   np.full_like(observations_corr, offset), rtol=1e-7, atol=0)


def test_correct_observations_light_deflection(gaia_astrometry):
    """Test that correct_observations correctly applies a light deflection offset to the observations"""
    get_obs = lambda table, mpc: table.loc[table['number_mp'] == mpc, ['ra', 'dec']].to_numpy()

    observation_table = gaia_astrometry.table
    observations_uncorr = {mpc: get_obs(observation_table, mpc) for mpc in TEST_ASTEROID_MPC}

    # Make light deflection function return a fixed offset
    offset = 1e-9  # 1e-9 radians in RA and DEC
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), offset)

    # Corrections are applied to all loaded asteroids
    with mock.patch('tudatpy.data.gaia.gaia.relativistic_light_deflection_from_observations',
                    side_effect=fake_correction):
        gaia_astrometry.correct_observations(bodies=None, light_deflection_bodies=['Sun'],
                                             correct_photocenter=False)
    observation_table = gaia_astrometry.table

    for mpc in TEST_ASTEROID_MPC:
        observations_corr = get_obs(observation_table, mpc)
        np.testing.assert_allclose(observations_corr - observations_uncorr[mpc],
                                   np.full_like(observations_corr, offset), rtol=1e-7, atol=0)


def test_correct_observations_twice_raises_error(gaia_astrometry):
    """Applying corrections twice on the same instance must raise an error"""
    fake_correction = lambda observations, **kwargs: np.zeros((len(observations), 2))

    with mock.patch('tudatpy.data.gaia.gaia.relativistic_light_deflection_from_observations',
                    side_effect=fake_correction):
        gaia_astrometry.correct_observations(bodies=None, light_deflection_bodies=['Sun'],
                                             correct_photocenter=False)

        with pytest.raises(RuntimeError):
            gaia_astrometry.correct_observations(bodies=None, light_deflection_bodies=['Sun'],
                                                 correct_photocenter=False)


def test_get_gaia_ephemeris_geocentric(gaia_astrometry, spice_kernels):
    """Test if states in catalog and ephemeris match (geocentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_astrometry.get_gaia_ephemeris(geocentric=True)
    body_settings = get_default_body_settings(['Sun', 'Earth'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get('Gaia').ephemeris

    # Compare state vectors from Tudat and table
    observation_table = gaia_astrometry.table
    states_from_tudat = [gaia_ephemeris.cartesian_state(epoch) for epoch in observation_table.epoch]
    states_from_tudat = np.array(states_from_tudat)
    states_from_table = observation_table[['x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                                           'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']]
    states_from_table = np.array(states_from_table)

    np.testing.assert_array_equal(states_from_table, states_from_tudat)


def test_get_gaia_ephemeris_barycentric(gaia_astrometry, spice_kernels):
    """Test if states in catalog and ephemeris match (barycentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_astrometry.get_gaia_ephemeris(geocentric=False)
    body_settings = get_default_body_settings(['Sun', 'Earth'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get('Gaia').ephemeris

    # Compare state vectors from Tudat and table
    observation_table = gaia_astrometry.table
    states_from_tudat = [gaia_ephemeris.cartesian_state(epoch) for epoch in observation_table.epoch]
    states_from_tudat = np.array(states_from_tudat)
    states_from_table = observation_table[['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']]
    states_from_table = np.array(states_from_table)

    np.testing.assert_array_equal(states_from_table, states_from_tudat)


def test_summary_does_not_raise(gaia_astrometry, capsys):
    """print_summary() should run without error and print a summary header"""
    gaia_astrometry.print_summary()

    assert 'SUMMARY' in capsys.readouterr().out


def test_to_observation_collection_without_gaia_raises_error(gaia_astrometry, spice_kernels):
    """to_observation_collection must raise an error if Gaia is not loaded in bodies"""
    body_settings = get_default_body_settings(['Sun'], 'SSB', 'J2000')
    bodies = create_system_of_bodies(body_settings)

    with pytest.raises(ValueError):
        gaia_astrometry.to_observation_collection(bodies)


def test_weight_matrix_symmetry(observation_collection):
    """Test that the weight matrix is symmetric"""
    weight_matrix = observation_collection.get_concatenated_weight_matrix().toarray()

    np.testing.assert_allclose(weight_matrix, weight_matrix.T, rtol=1e-10, atol=0)


def test_covariance_matrix_variance_consistency(observation_collection, observation_table):
    """Test that the covariance matrix diagonal is consistent with the observation table"""
    # Get variances from observation collection
    parser = observation_parser(str(TEST_ASTEROID_MPC[0]))  # We check for one asteroid
    weight_matrix = observation_collection.get_concatenated_weight_matrix(parser).toarray()
    covariance = np.linalg.inv(weight_matrix)
    variance_from_obscol = np.diag(covariance)

    # Get variances from table
    table = observation_table[observation_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    random_variance = np.ravel(table[['ra_error_random', 'dec_error_random']])
    sys_variance = np.ravel(table[['ra_error_systematic', 'dec_error_systematic']])
    variance_from_table = random_variance ** 2 + sys_variance ** 2

    np.testing.assert_allclose(variance_from_obscol, variance_from_table, rtol=1e-4, atol=0)


def test_covariance_matrix_nonzero_entries(observation_collection, observation_table):
    """Test block diagonal structure by checking number of nonzero entries in covariance matrix"""
    # Count number of nonzero entries in covariance matrix
    parser = observation_parser(str(TEST_ASTEROID_MPC[0]))  # We test for one asteroid
    weight_matrix = observation_collection.get_concatenated_weight_matrix(parser).toarray()
    covariance = np.linalg.inv(weight_matrix)
    nonzeros_covariance_matrix = np.count_nonzero(covariance)

    # Calculate expected number of nonzero entries for a block-diagonal structure
    table = observation_table[observation_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    transit_ids = pd.unique(table['transit_id'])
    transit_lengths = [2 * np.count_nonzero(table['transit_id'] == id) for id in transit_ids]  # times 2 because of RA, DEC
    nonzeros_expected = sum([n ** 2 for n in transit_lengths])

    assert nonzeros_covariance_matrix == nonzeros_expected


def test_observation_consistency(observation_collection, observation_table):
    """Test if observations in the table are consistent with those in observation_collection"""
    # Get observations and epochs from ObservationCollection
    parser = observation_parser(str(TEST_ASTEROID_MPC[0]))  # We check for one asteroid
    observations_from_collection = observation_collection.get_concatenated_observations(parser)
    epochs_from_collection = observation_collection.get_concatenated_observation_times(parser)

    # Get observations and epochs from table
    table = observation_table[observation_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    epochs_from_table = table.epoch.to_numpy()
    observations_from_table = np.ravel(table[['ra', 'dec']])

    np.testing.assert_array_equal(epochs_from_collection, epochs_from_table)
    np.testing.assert_array_equal(observations_from_collection, observations_from_table)
