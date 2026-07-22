from data.gaia_astrometry import GaiaQuery
import pandas.testing as pdt
import pytest
import numpy as np
from validation.utils import read_jpl_table
from postprocessing.helpers import interpolate_state_history
from tudatpy.astro.time_representation import date_time_components_to_epoch
from tudatpy.dynamics.environment_setup import get_default_body_settings, create_system_of_bodies
from tudatpy.interface import spice
from tudatpy.estimation.observations.observations_processing import observation_parser
import pandas as pd
import os

# Constants
TEST_ASTEROID_MPC = (433, 779) # Asteroids used in subsequent tests
LOCAL_FILE_PATH = 'data/SsoObservation/'
GAIA_SPKID = -139479

# ── FIXTURES ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope='session')
def spice_kernels():
    """Load spice kernels, necessary for setting up a SystemOfBodies"""
    spice.load_kernel('data/ephemerides/de441.bsp')
    spice.load_kernel('data/ephemerides/pck00011.tpc')
    spice.load_kernel('data/ephemerides/gm_de440.tpc')


@pytest.fixture(scope='module')
def _gaia_query():
    """Create a GaiaQuery object for asteroids 433 and 779 from local files. Runs once per module"""
    query = GaiaQuery()
    query.retrieve_data_locally(TEST_ASTEROID_MPC, LOCAL_FILE_PATH)
    return query


@pytest.fixture
def gaia_query(_gaia_query):
    """Copy of query object to avoid mutation in subsequent tests"""
    return _gaia_query.copy()


@pytest.fixture
def observation_table(gaia_query):
    """Shortcut for gaia_query.observation_table"""
    return gaia_query.observation_table


@pytest.fixture
def observation_collection(gaia_query, spice_kernels):
    """ObservationCollection of gaia_query constructed with to_tudat()"""
    # Must create bodies object to use to_tudat(), contents do not matter
    body_settings = get_default_body_settings(['Sun'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = gaia_query.get_gaia_ephemeris(False)
    bodies =  create_system_of_bodies(body_settings)
    observation_collection = gaia_query.to_tudat(bodies)
    return observation_collection

# ── TESTS ──────────────────────────────────────────────────────────────────

@pytest.mark.slow
def test_local_astroquery_consistency(observation_table):
    """Test whether data obtained from the local catalog and astroquery are consistent"""
    query_from_astroquery = GaiaQuery()
    query_from_astroquery.retrieve_data(TEST_ASTEROID_MPC)
    observation_table_from_astroquery = query_from_astroquery.observation_table

    pdt.assert_frame_equal(observation_table, observation_table_from_astroquery, check_dtype=False)


@pytest.mark.slow
def test_local_cached_archive_consistency():
    """Test consistency between local archive and cached archive"""
    # Check if cached archive exists
    cached_archive = LOCAL_FILE_PATH + 'gaia_astrometry_archive.parquet'
    assert os.path.exists(cached_archive), 'No cached archive found'

    # Retrieve from pkl
    query_from_cached = GaiaQuery()
    query_from_cached.retrieve_data_locally(TEST_ASTEROID_MPC, LOCAL_FILE_PATH)
    table_from_cached = query_from_cached.observation_table
    os.remove(cached_archive) # Allow to rebuild observation table from .csv files

    # Retrieve from csv
    query_from_csv = GaiaQuery()
    query_from_csv.retrieve_data_locally(TEST_ASTEROID_MPC, LOCAL_FILE_PATH)
    table_from_csv = query_from_csv.observation_table

    pdt.assert_frame_equal(table_from_cached, table_from_csv, check_dtype=False)

@pytest.mark.slow
def test_data_retrieval_no_data_astroquery():
    """Test that no data returned raises an error for the astroquery retrieve_data"""
    query = GaiaQuery()
    asteroid_no_data = 250000

    with pytest.raises(LookupError):
        query.retrieve_data([asteroid_no_data])

def test_data_retrieval_no_data_locally():
    """Test that local variant of retrieve_data raises error if no data is  found"""
    query = GaiaQuery()
    asteroid_no_data = 250000

    with pytest.raises(LookupError):
        query.retrieve_data_locally([asteroid_no_data], LOCAL_FILE_PATH)


def test_unavailable_catalog_raises_error():
    """Error must be raised if catalog is not available on Gaia archives"""
    query = GaiaQuery()

    with pytest.raises(ValueError):
        query.retrieve_data([433], catalog='DR5')


def test_observation_table_angles(observation_table):
    """Test units and range of angles in observation table"""
    # Check if angles are in radians and in correct range
    ra, dec, pa = observation_table[['ra', 'dec', 'position_angle_scan']].to_numpy().T

    assert np.all((ra < np.pi) & (ra > -np.pi))
    assert np.all((dec < np.pi/2) & (dec > -np.pi/2))
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

def test_gaia_barycentric_state_with_jpl_ephemeris(observation_table):
    """Test that the barycentric states of Gaia in observation_table is consistent with those from JPL Horizons"""
    # Retrieve observation table for one asteroid only
    test_asteroid = TEST_ASTEROID_MPC[0]
    observation_table_test_asteroid = observation_table[observation_table['number_mp'] == test_asteroid]
    epoch_list = observation_table_test_asteroid['epoch'].to_list()
    states_from_horizons,_ = read_jpl_table('tests/fixtures/gaia_ephemeris_jpl_barycentric_12h.txt')
    states_from_horizons = interpolate_state_history(states_from_horizons, epoch_list) # Interpolate to same times
    states_from_horizons = states_from_horizons[:, 1:]
    states_from_table = observation_table_test_asteroid[['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']]
    states_from_table = states_from_table.to_numpy()

    np.testing.assert_allclose(states_from_horizons, states_from_table, rtol=1e-05, atol=10)
    # We use looser tolerances here because a large deviation can exist due to different adopted planetary ephemerides

def test_gaia_geocentric_state_with_jpl_ephemeris(observation_table):
    """Test that the geocentric states of Gaia in observation_table is consistent with those from JPL Horizons"""
    # Retrieve observation table for one asteroid only
    test_asteroid = TEST_ASTEROID_MPC[0]
    observation_table_test_asteroid = observation_table[observation_table['number_mp'] == test_asteroid]
    epoch_list = observation_table_test_asteroid['epoch'].to_list()
    states_from_horizons,_ = read_jpl_table('tests/fixtures/gaia_ephemeris_jpl_geocentric_12h.txt')
    states_from_horizons = interpolate_state_history(states_from_horizons, epoch_list)  # interpolate to same times
    states_from_horizons = states_from_horizons[:, 1:] # Remove times from state history
    states_from_table = observation_table_test_asteroid[['x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                                           'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']]
    states_from_table = states_from_table.to_numpy()

    np.testing.assert_allclose(states_from_horizons, states_from_table, rtol=1e-05, atol=10)
    # We use looser tolerances here because a large deviation can exist due to different adopted planetary ephemerides

def test_observation_table_epoch_filter(gaia_query):
    """Test filtering of epochs with float bounds"""
    filter_start = date_time_components_to_epoch(year=2017, month=1, day=1, hour=0, minute=0, seconds=0)
    filter_end = date_time_components_to_epoch(year=2018, month=1, day=1, hour=0, minute=0, seconds=0)

    # Check first that observations outside of this time range exist
    assert gaia_query.epoch_start < filter_start
    assert gaia_query.epoch_end > filter_end

    gaia_query.filter(epoch_start=filter_start, epoch_end=filter_end)

    assert gaia_query.epoch_start >= filter_start
    assert gaia_query.epoch_end <= filter_end



def test_correct_observations_photocenter(gaia_query, mocker):
    """Tests that correct_observations correctly applies a photocenter offset to one asteroid's observations"""
    get_obs = lambda table, mpc: table.loc[table['number_mp'] == mpc, ['ra', 'dec']].to_numpy()

    test_asteroid = TEST_ASTEROID_MPC[0] # We add offset to this asteroid, other unchanged
    observation_table = gaia_query.observation_table
    observations_test_asteroid_uncorr = get_obs(observation_table, TEST_ASTEROID_MPC[0])
    observations_other_uncorr = get_obs(observation_table, TEST_ASTEROID_MPC[1])

    # Make photocenter offset function return a fixed offset
    mock = mocker.patch('Code.data.gaia_astrometry.photocenter_offset_spherical')
    number_of_obs = len(observations_test_asteroid_uncorr)
    photocenter_offset_fake = np.full((number_of_obs, 2), 1e-9) # 1e-9 radians in RA and DEC
    mock.return_value = photocenter_offset_fake

    # Now check that test asteroid observations were corrected, other's observations unchanged
    gaia_query.correct_observations(test_asteroid, bodies=mocker.MagicMock(), light_deflection=(), correct_photocenter=True)
    observation_table = gaia_query.observation_table
    observations_test_asteroid_corr = get_obs(observation_table, TEST_ASTEROID_MPC[0])
    observations_other_corr = get_obs(observation_table, TEST_ASTEROID_MPC[1])

    np.testing.assert_allclose(observations_test_asteroid_corr-observations_test_asteroid_uncorr, photocenter_offset_fake, rtol=1e-7, atol=0)
    np.testing.assert_array_equal(observations_other_corr, observations_other_uncorr)


def test_correct_observations_light_deflection(gaia_query, mocker):
    """Test that correct_observations correctly applies a light deflection offset to one asteroid's observations"""
    get_obs = lambda table, mpc: table.loc[table['number_mp'] == mpc, ['ra', 'dec']].to_numpy()

    test_asteroid = TEST_ASTEROID_MPC[0] # We add offset to this asteroid, other unchanged
    observation_table = gaia_query.observation_table
    observations_test_asteroid_uncorr = get_obs(observation_table, TEST_ASTEROID_MPC[0])
    observations_other_uncorr = get_obs(observation_table, TEST_ASTEROID_MPC[1])

    # Make photocenter offset function return a fixed offset
    mock = mocker.patch('Code.data.gaia_astrometry.relativistic_light_deflection')
    number_of_obs = len(observations_test_asteroid_uncorr)
    offset_fake = np.full((number_of_obs, 2), 1e-9) # 1e-9 radians in RA and DEC
    mock.return_value = offset_fake

    # Now check that test asteroid observations were corrected, other's observations unchanged
    gaia_query.correct_observations(test_asteroid, bodies=mocker.MagicMock(), light_deflection=('Sun'), correct_photocenter=False)
    observation_table = gaia_query.observation_table
    observations_test_asteroid_corr = get_obs(observation_table, TEST_ASTEROID_MPC[0])
    observations_other_corr = get_obs(observation_table, TEST_ASTEROID_MPC[1])

    np.testing.assert_allclose(observations_test_asteroid_corr-observations_test_asteroid_uncorr, offset_fake, rtol=1e-7, atol=0)
    np.testing.assert_array_equal(observations_other_corr, observations_other_uncorr)

def test_get_gaia_ephemeris_geocentric(gaia_query, spice_kernels):
    """Test if states in catalog and ephemeris match (geocentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_query.get_gaia_ephemeris(geocentric=True)
    body_settings = get_default_body_settings(['Sun', 'Earth'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get('Gaia').ephemeris

    # Compare state vectors from Tudat and table
    observation_table = gaia_query.observation_table
    states_from_tudat = [gaia_ephemeris.cartesian_state(epoch) for epoch in observation_table.epoch]
    states_from_tudat = np.array(states_from_tudat)
    states_from_table = gaia_query.observation_table[['x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                                                      'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']]
    states_from_table = np.array(states_from_table)

    np.testing.assert_array_equal(states_from_table, states_from_tudat)

def test_get_gaia_ephemeris_barycentric(gaia_query, spice_kernels):
    """Test if states in catalog and ephemeris match (barycentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_query.get_gaia_ephemeris(geocentric=False)
    body_settings = get_default_body_settings(['Sun', 'Earth'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get('Gaia').ephemeris

    # Compare state vectors from Tudat and table
    observation_table = gaia_query.observation_table
    states_from_tudat = [gaia_ephemeris.cartesian_state(epoch) for epoch in observation_table.epoch]
    states_from_tudat = np.array(states_from_tudat)
    states_from_table = gaia_query.observation_table[['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']]
    states_from_table = np.array(states_from_table)

    np.testing.assert_array_equal(states_from_table, states_from_tudat)


def test_summary_does_not_raise(gaia_query, capsys):
    """summary() should run without error and print a summary header"""
    gaia_query.summary()

    assert 'Summary' in capsys.readouterr().out


def test_weight_matrix_symmetry(observation_collection):
    """Test that the weight matrix is symmetric"""
    weight_matrix = observation_collection.get_concatenated_weight_matrix().toarray()

    np.testing.assert_allclose(weight_matrix, weight_matrix.T, rtol=1e-10, atol=0)


def test_covariance_matrix_variance_consistency(observation_collection, observation_table):
    """Test that the covariance matrix diagonal is consistent with the observation table"""
    # Get variances from observation collection
    parser = observation_parser(str(TEST_ASTEROID_MPC[0])) # We check for one asteroid
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
    parser = observation_parser(str(TEST_ASTEROID_MPC[0])) # We test for one asteroid
    weight_matrix = observation_collection.get_concatenated_weight_matrix(parser).toarray()
    covariance = np.linalg.inv(weight_matrix)
    nonzeros_covariance_matrix = np.count_nonzero(covariance)

    # Calculate expected number of nonzero entries for a block-diagonal structure
    table = observation_table[observation_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    transit_ids = pd.unique(table['transit_id'])
    transit_lengths = [2 * np.count_nonzero(table['transit_id'] == id) for id in transit_ids] # times 2 because of RA, DEC
    nonzeros_expected = sum([n ** 2 for n in transit_lengths])

    assert nonzeros_covariance_matrix == nonzeros_expected


def test_observation_consistency(observation_collection, observation_table):
    """Test if observations in observation_table are consistent with those in observation_collection"""
    # Get observations and epochs from observationCollection
    parser = observation_parser(str(TEST_ASTEROID_MPC[0])) # We check for one asteroid
    observations_from_collection = observation_collection.get_concatenated_observations(parser)
    epochs_from_collection = observation_collection.get_concatenated_observation_times(parser)

    # Get observations and epochs from table
    table = observation_table[observation_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    epochs_from_table = table.epoch.to_numpy()
    observations_from_table = np.ravel(table[['ra', 'dec']])

    np.testing.assert_array_equal(epochs_from_collection, epochs_from_table)
    np.testing.assert_array_equal(observations_from_collection, observations_from_table)