"""
Tests for the GaiaAstrometry and GaiaAsteroids classes in tudatpy.data.gaia.gaia.
"""
from pathlib import Path
from unittest import mock
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from tudatpy.data.gaia.gaia import GaiaAstrometry, GaiaAsteroids
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

TEST_ASTEROID_MPC = (673, 779)  # Asteroids used in subsequent tests (673 Edda and 779 Nina)
TEST_DIR = Path(__file__).parent
ASTROMETRY_ARCHIVE_PATH = TEST_DIR / 'gaia_astro_archive_for_tests.parquet'  # Astrometry for 673 and 779 only
SOURCE_ARCHIVE_PATH = TEST_DIR / 'gaia_source_archive_for_tests.parquet'  # Orbit data for 673 and 779 only


##############################
# GaiaAstrometry tests
##############################

@pytest.fixture(scope='module')
def spice_kernels():
    """Load spice kernels, necessary for setting up a SystemOfBodies"""
    spice.load_standard_kernels()


@pytest.fixture(scope='module')
def _gaia_astrometry():
    """Create a GaiaAstrometry object for asteroids 673 and 779 from the local test archive.
    Runs once per module"""
    return GaiaAstrometry.load_from_local_archive(ASTROMETRY_ARCHIVE_PATH, TEST_ASTEROID_MPC)


@pytest.fixture
def gaia_astrometry(_gaia_astrometry):
    """Copy of the GaiaAstrometry object to avoid mutation in subsequent tests"""
    return _gaia_astrometry.copy()


@pytest.fixture
def astrometry_table(gaia_astrometry):
    """Shortcut for gaia_astrometry.table"""
    return gaia_astrometry.table


@pytest.fixture
def observation_collection(gaia_astrometry, spice_kernels):
    """ObservationCollection of gaia_astrometry with observations for 673 and 779"""
    # Gaia must be loaded in bodies with its ephemeris to use to_observation_collection()
    body_settings = get_default_body_settings(['Sun'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = gaia_astrometry.get_gaia_ephemeris(geocentric=False)
    bodies = create_system_of_bodies(body_settings)
    return gaia_astrometry.to_observation_collection(bodies)


@pytest.mark.remote_data
def test_local_astroquery_consistency(astrometry_table):
    """Test whether data obtained from the local archive and astroquery are consistent"""
    astrometry_from_astroquery = GaiaAstrometry.load_from_astroquery(TEST_ASTEROID_MPC)
    astrometry_table_from_astroquery = astrometry_from_astroquery.table

    pdt.assert_frame_equal(astrometry_table, astrometry_table_from_astroquery, check_dtype=False)


@pytest.mark.remote_data
def test_astrometry_retrieval_no_data():
    """Test that no data available raises error"""
    asteroid_no_data = 250000

    # Online variant
    with pytest.raises(LookupError):
        GaiaAstrometry.load_from_astroquery([asteroid_no_data])
    # Local variant
    with pytest.raises(LookupError):
        GaiaAstrometry.load_from_local_archive(ASTROMETRY_ARCHIVE_PATH, [asteroid_no_data])


def test_observation_table_angles(astrometry_table):
    """Test units and range of angles in observation table"""
    # Check if angles are in radians and in correct range
    ra, dec, pa = astrometry_table[['ra', 'dec', 'position_angle_scan']].to_numpy().T

    assert np.all((ra < np.pi) & (ra > -np.pi))
    assert np.all((dec < np.pi / 2) & (dec > -np.pi / 2))
    assert np.all((pa < 2 * np.pi) & (pa > 0))

def test_table_ordering(astrometry_table):
    """Test ordering of table entries (first by MPC, then by epoch)"""
    # MPC ordering
    assert astrometry_table['number_mp'].is_monotonic_increasing

    # Epoch ordering
    for mpc in TEST_ASTEROID_MPC:
        assert astrometry_table.loc[astrometry_table['number_mp'] == mpc, 'epoch'].is_monotonic_increasing


def test_observation_table_uncertainty(astrometry_table):
    """Check that uncertainty and correlation are in the correct ranges"""
    unc = astrometry_table[['ra_error_random', 'dec_error_random', 'ra_error_systematic',
                             'dec_error_systematic']].to_numpy()
    corr = astrometry_table[['ra_dec_correlation_random', 'ra_dec_correlation_systematic']].to_numpy()

    assert np.all(unc >= 0)
    assert np.all((corr >= -1) & (corr <= 1))


def test_observation_table_epochs(astrometry_table):
    """Test if the epochs are in the range of the epochs for Gaia FPR"""
    earliest_epoch = date_time_components_to_epoch(year=2014, month=7, day=26, hour=0, minute=0, seconds=0)
    final_epoch = date_time_components_to_epoch(year=2020, month=1, day=20, hour=0, minute=0, seconds=0)
    epochs = astrometry_table['epoch'].to_numpy()

    assert np.all((epochs >= earliest_epoch) & (epochs <= final_epoch))


# NOTE: Gaia ephemeris tests are loose because possibly different planetary ephemerides models were used for
# the JPL/Gaia orbit estimations
def test_gaia_barycentric_state_with_jpl_ephemeris(astrometry_table):
    """Test that the barycentric states of Gaia in the table are consistent with those from JPL Horizons"""
    columns = ['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']
    states_from_archive = astrometry_table.sort_values(by='epoch')[columns].to_numpy()
    states_from_jpl = np.loadtxt(TEST_DIR / 'gaia_ephemeris_barycentric.txt')[:, 1:]

    np.testing.assert_allclose(states_from_archive, states_from_jpl,  rtol=1e-4, atol=1e-4)


def test_gaia_geocentric_state_with_jpl_ephemeris(astrometry_table):
    """Test that the geocentric states of Gaia in the table are consistent with those from JPL Horizons"""
    columns = ['x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
               'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']
    states_from_archive = astrometry_table.sort_values(by='epoch')[columns].to_numpy()
    states_from_jpl = np.loadtxt(TEST_DIR / 'gaia_ephemeris_geocentric.txt')[:, 1:]

    np.testing.assert_allclose(states_from_archive, states_from_jpl, rtol=1e-4, atol=1e-4)


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
    """Tests that correct_observations correctly modifies observation by calculated photocenter offset"""
    get_obs = lambda table, mpc: table.loc[table['number_mp'] == mpc, ['ra', 'dec']].to_numpy()

    astrometry_table = gaia_astrometry.table
    observations_uncorr = {mpc: get_obs(astrometry_table, mpc) for mpc in TEST_ASTEROID_MPC}

    # Make photocenter offset function return a fixed offset
    offset = 1e-9  # 1e-9 radians in RA and DEC
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), offset)

    # Corrections are applied to all loaded asteroids
    diameters = {mpc: 1e3 for mpc in TEST_ASTEROID_MPC}
    with mock.patch('tudatpy.data.gaia.gaia.photocenter_corrections_from_observations',
                    side_effect=fake_correction):
        gaia_astrometry.correct_observations(bodies=None, diameters=diameters,
                                             light_deflection_bodies=[], correct_photocenter=True)
    astrometry_table = gaia_astrometry.table

    for mpc in TEST_ASTEROID_MPC:
        observations_corr = get_obs(astrometry_table, mpc)
        np.testing.assert_allclose(observations_corr - observations_uncorr[mpc],
                                   np.full_like(observations_corr, offset), rtol=1e-7, atol=0)


def test_correct_observations_light_deflection(gaia_astrometry):
    """Test that correct_observations correctly applies a light deflection offset to the observations"""
    get_obs = lambda table, mpc: table.loc[table['number_mp'] == mpc, ['ra', 'dec']].to_numpy()

    astrometry_table = gaia_astrometry.table
    observations_uncorr = {mpc: get_obs(astrometry_table, mpc) for mpc in TEST_ASTEROID_MPC}

    # Make light deflection function return a fixed offset
    offset = 1e-9  # 1e-9 radians in RA and DEC
    fake_correction = lambda observations, **kwargs: np.full((len(observations), 2), offset)

    # Corrections are applied to all loaded asteroids
    with mock.patch('tudatpy.data.gaia.gaia.relativistic_light_deflection_from_observations',
                    side_effect=fake_correction):
        gaia_astrometry.correct_observations(bodies=None, light_deflection_bodies=['Sun'],
                                             correct_photocenter=False)
    astrometry_table = gaia_astrometry.table

    for mpc in TEST_ASTEROID_MPC:
        observations_corr = get_obs(astrometry_table, mpc)
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
    """Test if states in catalog and those retrieved from ephemeris match (geocentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_astrometry.get_gaia_ephemeris(geocentric=True)
    body_settings = get_default_body_settings(['Sun', 'Earth'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get('Gaia').ephemeris

    # Compare state vectors from Tudat and table
    astrometry_table = gaia_astrometry.table
    states_from_tudat = np.array([gaia_ephemeris.cartesian_state(epoch) for epoch in astrometry_table.epoch])
    states_from_table = astrometry_table[['x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                                'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']].to_numpy()

    np.testing.assert_array_equal(states_from_table, states_from_tudat)


def test_get_gaia_ephemeris_barycentric(gaia_astrometry, spice_kernels):
    """Test if states in catalog and those retrieved from ephemeris match (barycentric case)"""
    # Construct Tudat ephemeris
    ephemeris_settings = gaia_astrometry.get_gaia_ephemeris(geocentric=False)
    body_settings = get_default_body_settings(['Sun', 'Earth'], 'SSB', 'J2000')
    body_settings.add_empty_settings('Gaia')
    body_settings.get('Gaia').ephemeris_settings = ephemeris_settings
    bodies = create_system_of_bodies(body_settings)
    gaia_ephemeris = bodies.get('Gaia').ephemeris

    # Compare state vectors from Tudat and table
    astrometry_table = gaia_astrometry.table
    states_from_tudat = np.array([gaia_ephemeris.cartesian_state(epoch) for epoch in astrometry_table.epoch])
    states_from_table = astrometry_table[['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']].to_numpy()

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


def test_covariance_matrix_variance_consistency(gaia_astrometry):
    """Test that the covariance matrix diagonal matches the per-observation variances in the
    astrometry table, ordered by epoch (RA then DEC per observation)."""
    covariance = gaia_astrometry.get_observation_covariance_matrix(TEST_ASTEROID_MPC[0])
    variance_from_matrix = np.diag(covariance)

    table = gaia_astrometry.table_for_mpc(TEST_ASTEROID_MPC[0])

    # By design, the systematic error is constant per transit: use the first entry of each transit_id
    systematic = table.groupby('transit_id', sort=False)[
        ['ra_error_systematic', 'dec_error_systematic']].transform('first')

    var_ra = table['ra_error_random'] ** 2 + systematic['ra_error_systematic'] ** 2
    var_dec = table['dec_error_random'] ** 2 + systematic['dec_error_systematic'] ** 2
    variance_from_table = np.column_stack([var_ra, var_dec]).ravel()

    np.testing.assert_array_equal(variance_from_matrix, variance_from_table)


def test_covariance_weight_matrix_consistency(observation_collection, gaia_astrometry):
    """Test that the covariance matrix computed from the astrometry table is consistent with weight matrix
    passed to observation_collection"""
    # Get variances from observation collection
    parser = observation_parser(str(TEST_ASTEROID_MPC[0]))  # We check for one asteroid
    weight_matrix = observation_collection.get_concatenated_weight_matrix(parser).toarray()
    covariance_from_obscol = np.linalg.inv(weight_matrix)
    covariance_from_table = gaia_astrometry.get_observation_covariance_matrix(TEST_ASTEROID_MPC[0])

    # NOTE: loose tolerance because of high condition number so much precision is lost during 2 inversions
    np.testing.assert_allclose(covariance_from_obscol, covariance_from_table, rtol=1e-5, atol=0)

def test_covariance_matrix_nonzero_entries(observation_collection, astrometry_table):
    """Test block diagonal structure by checking number of nonzero entries in covariance matrix"""
    # Count number of nonzero entries in covariance matrix
    parser = observation_parser(str(TEST_ASTEROID_MPC[0]))  # We test for one asteroid
    weight_matrix = observation_collection.get_concatenated_weight_matrix(parser).toarray()
    covariance = np.linalg.inv(weight_matrix)
    nonzeros_covariance_matrix = np.count_nonzero(covariance)

    # Calculate expected number of nonzero entries for a block-diagonal structure
    table = astrometry_table[astrometry_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    transit_ids = pd.unique(table['transit_id'])
    transit_lengths = [2 * np.count_nonzero(table['transit_id'] == id) for id in transit_ids]  # times 2 because of RA, DEC
    nonzeros_expected = sum([n ** 2 for n in transit_lengths])

    assert nonzeros_covariance_matrix == nonzeros_expected


def test_observation_consistency(observation_collection, astrometry_table):
    """Test if observations in the table are consistent with those in observation_collection"""
    # Get observations and epochs from ObservationCollection
    parser = observation_parser(str(TEST_ASTEROID_MPC[0]))  # We check for one asteroid
    observations_from_collection = observation_collection.get_concatenated_observations(parser)
    epochs_from_collection = observation_collection.get_concatenated_observation_times(parser)

    # Get observations and epochs from table
    table = astrometry_table[astrometry_table['number_mp'] == TEST_ASTEROID_MPC[0]]
    epochs_from_table = table.epoch.to_numpy()
    observations_from_table = np.ravel(table[['ra', 'dec']])

    np.testing.assert_array_equal(epochs_from_collection, epochs_from_table)
    np.testing.assert_array_equal(observations_from_collection, observations_from_table)


##############################
# GaiaAsteroids tests
##############################

# Reference J2000 Heliocentric state vectors retrieved from JPL Horizons at the same epoch as Gaia catalog
STATE_EDDA_JPL = [ 4.128085144236671E+08 , 7.098887570553194E+07 , 4.413226270395131E+07,
 -3.659135280911876E+00 , 1.621253557540478E+01 , 6.231297081940196E+00] # KM, KM/s
STATE_NINA_JPL = [2.143720589223085E+08 , 2.310721906203461E+08 , 1.782520574512776E+08,
 -1.254695764130771E+01 , 1.502661150744586E+01 , 4.097474113261074E+00]


@pytest.fixture(scope='module')
def _gaia_asteroids():
    """Create an instance of the GaiaAsteroids class with the test asteroids loaded"""
    return GaiaAsteroids.load_from_local_archive(SOURCE_ARCHIVE_PATH, TEST_ASTEROID_MPC)

@pytest.fixture
def gaia_asteroids(_gaia_asteroids):
    """Get a copy of the GaiaAsteroids instance"""
    return _gaia_asteroids.copy()


@pytest.mark.remote_data
def test_asteroid_local_astroquery_consistency(gaia_asteroids):
    """Test whether asteroid data obtained from the local archive and astroquery are consistent"""
    table_local = gaia_asteroids.table.sort_values('number_mp').reset_index(drop=True)
    asteroids_from_astroquery = GaiaAsteroids.load_from_astroquery(TEST_ASTEROID_MPC)
    table_astroquery = asteroids_from_astroquery.table.sort_values('number_mp').reset_index(drop=True)

    pdt.assert_frame_equal(table_local, table_astroquery, check_dtype=False)


@pytest.mark.remote_data
def test_asteroid_data_retrieval_no_data():
    """Test that no data returned raises an error for GaiaAsteroids.load_from_astroquery"""
    asteroid_no_data = 250000

    # Astroquery variant
    with pytest.raises(LookupError):
        GaiaAsteroids.load_from_astroquery([asteroid_no_data])

    # Local variant
    with pytest.raises(LookupError):
        GaiaAsteroids.load_from_local_archive(SOURCE_ARCHIVE_PATH, [asteroid_no_data])


def test_asteroid_epochs(gaia_asteroids):
    """Test that the state vector epochs fall within the Gaia FPR observation period"""
    earliest_epoch = date_time_components_to_epoch(year=2014, month=7, day=26, hour=0, minute=0, seconds=0)
    final_epoch = date_time_components_to_epoch(year=2020, month=1, day=20, hour=0, minute=0, seconds=0)
    epochs = gaia_asteroids.table['epoch_state_vector'].to_numpy()

    assert np.all((epochs >= earliest_epoch) & (epochs <= final_epoch))


def test_asteroid_state_vector_with_jpl_horizons(gaia_asteroids):
    """Test match between Gaia state vectors and that of JPL Horizons"""
    ref_states = [np.array(STATE_EDDA_JPL) * 1e3,
                  np.array(STATE_NINA_JPL) * 1e3] # JPL Horizons vectors
    for mpc_number, ref_state in zip(TEST_ASTEROID_MPC, ref_states):
        # Heliocentric state vectors:
        epoch, state = gaia_asteroids.get_state_for_mpc(mpc_number)

        np.testing.assert_allclose(state, ref_state, rtol=1e-7, atol=1e-7)


def test_asteroid_orbit_class(gaia_asteroids):
    """Both test asteroids are main-belt asteroids"""
    assert (gaia_asteroids.table['orbit_class'] == 'MB').all()


def test_asteroid_covariance_shape_and_symmetry(gaia_asteroids):
    """Covariance matrices must be reconstructed from the raw upper triangle as symmetric 6x6 matrices"""
    for _, asteroid_data in gaia_asteroids.table.iterrows():
        for column in ['orbital_elements_var_covar_matrix', 'h_state_vector_var_covar_matrix']:
            covariance = asteroid_data[column]

            assert covariance.shape == (6, 6)
            np.testing.assert_allclose(covariance, covariance.T, rtol=1e-12, atol=0)


def test_asteroid_covariance_values(gaia_asteroids):
    """Sanity checks on the covariance matrices: variances must be positive, and the semi-major axis
    uncertainty must be in a range plausible for Gaia orbit solutions (order of meters, not AU)"""
    for _, asteroid_data in gaia_asteroids.table.iterrows():
        for column in ['orbital_elements_var_covar_matrix', 'h_state_vector_var_covar_matrix']:
            assert np.all(np.diag(asteroid_data[column]) > 0)

        sma_uncertainty = np.sqrt(asteroid_data['orbital_elements_var_covar_matrix'][0, 0])
        assert 1e-2 < sma_uncertainty < 1e6  # In meters; catches missing or double unit scaling
