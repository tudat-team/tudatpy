"""
Retrieve Gaia FPR astrometry from the archives
"""
import numpy as np
from astroquery.gaia import Gaia
import pandas as pd
from tudatpy.astro.time_representation import julian_day_to_seconds_since_epoch, TCB_to_TDB, DateTime
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import model_settings, links
from tudatpy.constants import ASTRONOMICAL_UNIT,JULIAN_DAY
from tudatpy.dynamics.environment_setup import ephemeris
from datetime import datetime
from scipy.linalg import block_diag
from scipy.constants import arcsec
from tudatpy.estimation.observable_models_setup.biases.photocenter_correction import photocenter_offset_spherical
from tudatpy.estimation.observable_models_setup.biases.light_deflection_correction import relativistic_light_deflection
import copy
from pathlib import Path
from tudatpy.astro.element_conversion import j2000_to_eclipj2000, cartesian_to_keplerian
import ast

# Constants
_J2010 = 2455197.5  # Reference time J2010.0 in Julian days
_TIME_SCALE_CORRECTION = 1 - 1.550519768e-8 # See e.g. Klioner (2003)
_ASTROMETRY_TABLE_COLUMNS = ['epoch', 'ra', 'dec', 'number_mp', 'transit_id', 'position_angle_scan',
                   'ra_error_systematic', 'dec_error_systematic', 'ra_dec_correlation_systematic',
                   'ra_error_random', 'dec_error_random', 'ra_dec_correlation_random',
                   'x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia',
                   'x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                   'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']
_ASTROMETRY_CATALOG_NAMES = {'DR2': 'gaiadr2.sso_observation', 'DR3': 'gaiadr3.sso_observation', 'FPR': 'gaiafpr.sso_observation'}
_ASTEROID_TABLE_COLUMNS = [
    'number_mp', # MPC number
    'epoch_state_vector', # Epoch in TCB
    'orbital_elements_var_covar_matrix', # Cov. matrix of heliocentric orbital elements as (a, e, i, node, peri, M)
    'h_state_vector',
    'h_state_vector_var_covar_matrix'
]

_STATE_SCALING_FACTOR = 1.0000000051686297 # account for Gaia FPR state vector inconsistency
_ASTEROID_CATALOG_NAMES = {'DR2': 'gaiadr2.sso_source',
                 'DR3': 'gaiadr3.sso_source',
                 'FPR': 'gaiafpr.sso_source'}

class GaiaAstrometry:

    def __init__(self,):
        """
        Create a Gaia query object
        """
        self._table = pd.DataFrame() # Holds astrometry and metadata after calling retrieve_data
        self._measurement_covariance = {}
        self._corrected = False # Indicates if corrections have been applied

    @property
    def table(self):
        """
        Return a copy of the observation table

        :return:
        """
        return self._table.copy()

    @property
    def epoch_start(self):
        """
        Return first epoch in the observation table (any object)

        :return:
        """
        return self.table['epoch'].iloc[0] # Table is ordered by epoch

    @property
    def epoch_end(self):
        """
        Return last epoch in the observation table (any object)

        :return:
        """
        return self.table['epoch'].iloc[-1] # Table is ordered by epoch

    @property
    def mpc_numbers(self):
        """
        Which asteroid MPC numbers appear in the observations table
        """
        return pd.unique(self._table['number_mp'])

    def copy(self):
        """
        Return a copy of the query object

        Returns:
            GaiaAstrometry: Deep copy of the object
        """
        return copy.deepcopy(self)

    def copy_for_mpc(self,
                     mpc_numbers: int|list|tuple):
        """
        Return a copy of the query object for a selection of asteroids. Can be useful to create multiple query
        objects of individual asteroids without calling retrieve_data every time.

        Args:
            mpc_numbers (list): MPC numbers to keep in copy

        Returns:
            GaiaAstrometry: Copy containing only data for mpc_numbers given
        """
        if isinstance(mpc_numbers, int):
            mpc_numbers = [mpc_numbers]

        if any(m not in self.mpc_numbers for m in mpc_numbers):
            raise ValueError('Requested MPC numbers are not in loaded observation table')

        new = self.copy()
        new._table = new._table.query('number_mp in @mpc_numbers')
        new._measurement_covariance = {}

        return new

    def to_tudat(self,
                 bodies: SystemOfBodies):
        """
        Convert observations into Tudat ObservationCollection

        Args:
            bodies: SystemOfBodies object

        Returns:
            ObservationCollection: Tudat ObservationCollection containing observations of all asteroids
        """
        if self._table.empty:
            raise RuntimeError('No observations loaded')

        # Check if Gaia is in bodies
        if not bodies.does_body_exist('Gaia') or bodies.get('Gaia').ephemeris is None:
            raise ValueError('Gaia satellite and associated ephemeris must be loaded in SystemOfBodies')

        observation_set_list = []
        for mpc_number in self.mpc_numbers:

            # Add asteroids to bodies
            if not bodies.does_body_exist(str(mpc_number)):
                bodies.create_empty_body(str(mpc_number))

            # Create link ends for each asteroid
            link_ends = {}
            link_ends[links.transmitter] = links.body_origin_link_end_id(str(mpc_number))
            link_ends[links.receiver] = links.body_origin_link_end_id('Gaia')
       #     link_definition = links.link_definition(link_ends)

            # Get the data for current asteroid
            temp = self.table.query('number_mp == @mpc_number')
            observation_angles = temp.loc[:, ['ra', 'dec']].to_numpy()
            observation_times = temp['epoch'].to_numpy()

            observation_set = observations.create_single_observation_set(
                model_settings.angular_position_type,
                link_ends,
                observation_angles,
                observation_times,
                links.receiver
            )

            # Transit IDs for current asteroid
            assert temp['transit_id'].is_monotonic_increasing # Ensure no faulty transit_ids
            transit_ids_unique = pd.unique(temp['transit_id'])

            # Build covariance matrix
            measurement_covariance_matrix = [] # Matrix over full observation set
            for id in transit_ids_unique:

                transit = temp[temp['transit_id'] == id] # Observations in transit
         #       assert 2 <= len(transit) <= 9, 'Transit length off'

                # Error components
                uncertainty_random = transit[
                    ['ra_error_random',
                     'dec_error_random',
                     'ra_dec_correlation_random']
                ]
                sigma_ra_r, sigma_dec_r, corr_r = uncertainty_random.to_numpy().T

                sigma_ra_s = np.mean(transit['ra_error_systematic']) # Take means since there are small diffs due to numerical instability (?)
                sigma_dec_s = np.mean(transit['dec_error_systematic'])
                corr_s = np.mean(transit['ra_dec_correlation_systematic'])

                # Random component for observation AF1-9
                covariance_random = [
                    np.array([[sigma_ra_r[ii] ** 2, corr_r[ii]*sigma_ra_r[ii]*sigma_dec_r[ii]],
                              [corr_r[ii] * sigma_ra_r[ii] * sigma_dec_r[ii], sigma_dec_r[ii]**2]])
                    for ii in range(len(transit))
                ]
                covariance_random = block_diag(*covariance_random)
          #      assert len(covariance_random) == 2 * len(transit)

                # Systematic component over transit
                covariance_systematic_sub = np.array(
                    [[sigma_ra_s**2, corr_s * sigma_ra_s * sigma_dec_s],
                     [corr_s * sigma_ra_s * sigma_dec_s, sigma_dec_s ** 2]]
                )
                covariance_systematic = np.tile(covariance_systematic_sub, (len(transit), len(transit)))

             #   assert covariance_random.size == covariance_systematic.size, 'Size mismatch'
                measurement_covariance_matrix.append(covariance_random + covariance_systematic)

            # Form full weight matrix for all observations of asteroid
            measurement_covariance_matrix = block_diag(*measurement_covariance_matrix)
            weight_matrix = np.linalg.inv(measurement_covariance_matrix)
            self._measurement_covariance[mpc_number] = measurement_covariance_matrix

            observation_set.set_full_weight_matrix(weight_matrix)

            # Add SingleObservationSet to list
            observation_set_list.append(observation_set)

        observation_collection = observations.ObservationCollection(
            observation_set_list
        )

        return observation_collection

    def correct_observations(self,
                             mpc_numbers: list[int],
                             bodies: SystemOfBodies,
                             light_deflection: tuple | list = ('Sun',),
                             correct_photocenter: bool = True):
        """
        Apply corrections to the observations (in-place).

        Args:
            mpc_numbers: List of asteroids to apply corrections to
            bodies (SystemOfBodies): Bodies object which must have appropriate ephemerides loaded
            light_deflection (list or tuple): Body objects which exert relativistic light bending
            correct_photocenter (bool): Apply a photocenter correction

        Returns:
            None
        """
        if self._corrected: # Check if this function has already been used for current object
            raise RuntimeError('Observations: corrections already applied.')

        if not isinstance(mpc_numbers, (list,tuple)):
            mpc_numbers = [mpc_numbers]

        # Apply correction to specified asteroid bodies
        for mpc_number in mpc_numbers:
            # Check if correction is possible for this asteroid
            if mpc_number not in self.mpc_numbers:
                raise ValueError(f'Correction not possible for {mpc_number}, no observations loaded')
            if not bodies.does_body_exist(str(mpc_number)) or bodies.get(str(mpc_number)).ephemeris is None:
                raise ValueError(f'Correction not possible for {mpc_number}. Body needs to be in SystemOfBodies'
                                     f'with loaded ephemeris')

            mask = self.table['number_mp'] == mpc_number
            # Apply photocenter corrections to the table
            if correct_photocenter:
                photocenter_corrections = photocenter_offset_spherical(mpc_number, self.table, bodies)
                self._table.loc[mask, ['ra', 'dec']] += photocenter_corrections

            # Apply light-bending corrections
            if light_deflection:
                if not all(bodies.does_body_exist(body) for body in light_deflection):
                    raise ValueError('Light deflection bodies missing from bodies object')
                light_bending_offsets = relativistic_light_deflection(mpc_number, self.table,
                                                                      bodies, light_deflection)
                self._table.loc[mask, ['ra', 'dec']] += light_bending_offsets

        self._corrected = True

    def retrieve_data(self,
                      mpc_numbers: tuple[int] | list[int],
                      catalog: str = 'FPR',
                      username:str = None,
                      password:str = None,):

        """
        Retrieve the astrometric observations through astroquery. Observations are stored in the observation table
        attribute.

        Args:
            mpc_numbers (tuple, list): List of asteroid MPC numbers to retrieve
            catalog (str): Which catalog to use. Options: DR2, DR3, FPR
            username (str): Username for the Gaia archives (optional)
            password (str): Password for the Gaia archives (optional)
        """
        if not isinstance(mpc_numbers, (list,tuple)):
            mpc_numbers = [mpc_numbers]

        if catalog not in _ASTROMETRY_CATALOG_NAMES.keys():
            raise ValueError(f'Catalog not available. Catalog options are: {", ".join(_ASTROMETRY_CATALOG_NAMES.keys())}')

        # Define query to database
        query_mpc_numbers = ", ".join(str(mpc_number) for mpc_number in mpc_numbers)
        query_columns = ', '.join(_ASTROMETRY_TABLE_COLUMNS)
        query_catalog = _ASTROMETRY_CATALOG_NAMES[catalog]

        login_provided = username is not None and password is not None

        if login_provided:
            Gaia.login(user=username, password=password)

        try:
            query = f"""
            SELECT {query_columns}
            FROM {query_catalog}
            WHERE number_mp IN ({query_mpc_numbers}) 
            AND astrometric_outcome_ccd = 1
            AND astrometric_outcome_transit = 1
            ORDER BY epoch ASC
            """
            job = Gaia.launch_job_async(query)
            table = job.get_results()

        except Exception as err:
            raise RuntimeError(f'Error while retrieving astrometric observations: \n{err}') from err

        table = table.to_pandas() # Convert astropy table to dataframe
        if table.empty:
            raise LookupError(f'No observations found for mpc numbers {mpc_numbers}')

        # Pre-process data
        table = self._convert_units(table)
        table = table.reset_index(drop=True)
        assert table['epoch'].is_monotonic_increasing # Sanity check for ordering by epoch

        # Store
        self._table = table

    def retrieve_data_locally(self,
                              mpc_numbers: tuple[int] | list[int],
                              archive_file_path: str | Path):
        """
        Retrieve astrometry locally from .parquet file (generated by generate_parquet). Mirrors retrieve_data.

        Args:
            mpc_numbers (tuple, list): List of asteroid MPC numbers to retrieve
            archive_file_path (str, Path): Path to archive .parquet file
        """
        if not isinstance(mpc_numbers, (list,tuple)):
            mpc_numbers = [mpc_numbers]

        # Read from parquet
        filters = [
            ('number_mp', 'in', mpc_numbers),
            ('astrometric_outcome_ccd', '==', 1),
            ('astrometric_outcome_transit', '==', 1),
        ]
        table = pd.read_parquet(archive_file_path, filters=filters)

        if table.empty:
            raise LookupError(f'No observations found for mpc numbers {mpc_numbers}')

        table = table[_ASTROMETRY_TABLE_COLUMNS]
        table = table.sort_values(by='epoch')

        # Convert units and save
        table = self._convert_units(table)
        table = table.reset_index(drop=True)

        self._table = table

    @staticmethod
    def generate_parquet(archive_dir: Path | str,
                         dir_to_save: Path | str):
        """
        Generate a parquet file from CSV files

        Args:
            archive_dir (Path | str): Path to the archive directory
            dir_to_save (Path | str): Path to the directory to save the parquet file

        Returns:
            None
        """
        print('Loading archive from CSV files\n')

        table = pd.DataFrame()

        for i in range(20):
            print(f'Loading file number {i}/19...')
            # Get current file path
            file_number = str(i) if i >= 10 else '0' + str(i)
            file_name = 'SsoObservation_' + file_number + '.csv'
            csv_file = Path(archive_dir) / file_name
            # Concatenate into one dataframe
            table_chunk = pd.read_csv(csv_file,
                                      comment='#')
            table = pd.concat([table, table_chunk], ignore_index=True)

        # Save to parquet
        table.to_parquet(Path(dir_to_save) / 'gaia_astrometry_archive.parquet')


    def _convert_units(self,
                       table: pd.DataFrame) -> pd.DataFrame:
        """
        Convert the table columns into correct format for Tudat

        :return:
        """
        # Convert epoch to seconds since J2000
        func = lambda jd: julian_day_to_seconds_since_epoch(jd + _J2010)
        table['epoch'] = table['epoch'].apply(func)

        # Convert TCB to TDB epoch
        table['epoch'] = table['epoch'].apply(TCB_to_TDB)

        # Convert angles to rad
        table['ra'] = (np.deg2rad(table['ra']) + np.pi) % (2 * np.pi) - np.pi
        table['dec'] = np.deg2rad(table['dec'])
        table['position_angle_scan'] = np.deg2rad(table['position_angle_scan'])

        # Convert mas to radians
        table[['ra_error_random','dec_error_random','ra_error_systematic','dec_error_systematic']] *= arcsec/1e3

        # Remove the cos delta factor from the right ascension uncertainty values
        table['ra_error_random'] /= np.cos(table['dec'])
        table['ra_error_systematic'] /= np.cos(table['dec'])

        # Convert Gaia state vectors to SI, apply correction to position vectors due to time scale change
        pos_names = ['x_gaia', 'y_gaia', 'z_gaia', 'x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric']
        table.loc[:, pos_names] *= ASTRONOMICAL_UNIT * _TIME_SCALE_CORRECTION
        vel_names = ['vx_gaia', 'vy_gaia', 'vz_gaia', 'vx_gaia_geocentric', 'vy_gaia_geocentric',
                     'vz_gaia_geocentric']
        table.loc[:, vel_names] *= ASTRONOMICAL_UNIT / DAY_IN_S

        return table


    def filter(self,
               epoch_start: float | datetime | DateTime,
               epoch_end: float | datetime | DateTime) -> None:
        """
        Filter the observations after they have been loaded (in-place)
        :return:
        """
        # Check if observations have been loaded
        if len(self._table) == 0:
            raise Exception('No observations loaded')

        # Convert parameters to seconds since J2000
        if isinstance(epoch_start, datetime) and isinstance(epoch_end, datetime):
            epoch_start = DateTime.from_python_datetime(epoch_start)
            epoch_end = DateTime.from_python_datetime(epoch_end)

        if isinstance(epoch_start, DateTime) and isinstance(epoch_end, DateTime):
            epoch_start = epoch_start.to_epoch()
            epoch_end = epoch_end.to_epoch()

        # Find observations in time span
        obs_in_timespan = self._table.query('@epoch_start <= epoch <= @epoch_end')

        if not obs_in_timespan.empty:
            self._table = obs_in_timespan

        else:
            raise Exception('No observations left after filtering')


    def summary(self) -> None:
        """
        Print a convenient summary of the astrometric observations
        :return:
        """
        if len(self._table) == 0:
            print('Observations not loaded')

        # Check mpc numbers in table (after filtering)
        print('Summary:')
        print(f'Observations for {len(self.mpc_numbers)} objects:')

        first_epoch = DateTime.from_epoch(self.epoch_start)
        final_epoch = DateTime.from_epoch(self.epoch_end)
        print(f'First observation yy/mm/dd: {first_epoch.year}, {first_epoch.month}, {first_epoch.day}')
        print(f'Final observation yy/mm/dd: {final_epoch.year}, {final_epoch.month}, {final_epoch.day}')

        for mpc_number in self.mpc_numbers:

            print(f'\nMinor planet {mpc_number}:')
            table_single_obj = self.table.query('number_mp == @mpc_number')

            nr_of_observations = len(table_single_obj)
            print(f'Number of observations: {nr_of_observations}')

            epochs_as_list = table_single_obj['epoch'].to_list()
            first_epoch = DateTime.from_epoch(epochs_as_list[0])
            final_epoch = DateTime.from_epoch(epochs_as_list[-1])

            print(f'First observation yy/mm/dd: {first_epoch.year}, {first_epoch.month}, {first_epoch.day}')
            print(f'Final observation yy/mm/dd: {final_epoch.year}, {final_epoch.month}, {final_epoch.day}')


    def get_gaia_ephemeris(self,
                           geocentric = True):
        """
        Get tabulated ephemeris settings generated from the achived Gaia state vectors

        Args:
            geocentric (bool): If true, use geocentric Gaia states (recommended). If false, use barycentric states

        Returns:
            EphemerisSettings: Tabulated ephemeris settings
        """
        # Variable names for the state vector
        state_vector_labels = ['x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia']
        if geocentric:
            state_vector_labels = [label + '_geocentric' for label in state_vector_labels]

        # Create dict of state vectors
        table = self.table
        epochs = table['epoch'].to_numpy()
        states = table[state_vector_labels].to_numpy()
        gaia_state_history = dict(zip(epochs, states))

        settings = ephemeris.tabulated(
            gaia_state_history,
            frame_origin='Earth' if geocentric else 'SSB',
            frame_orientation = 'J2000'
        )
        return settings


class GaiaAsteroids:

    def __init__(self):
        """
        GaiaOrbits collects states and covariance as calculated by Gaia DPAC
        """
        self._table = pd.DataFrame()

    @property
    def table(self) -> pd.DataFrame:
        """
        Return a copy of orbits table
        """
        return self._table.copy()

    def get_state_for_mpc(self,
                          mpc_number: int):
        """
        Retrieve state vector from the table for a single MPC number

        Args:
            mpc_number (int): Asteroid to retrieve data for

        Returns:
            float: epoch,
            np.ndarray: State vector (heliocentric J2000)
        """
        filtered = self.table.query('number_mp == @mpc_number').iloc[0]
        return filtered.epoch_state_vector, filtered.h_state_vector


    def retrieve_data(self,
                      mpc_numbers: int | list | tuple,
                      catalog: str = 'FPR'):
        """
        Retrieve orbital solutions from astroquery (SLOW)

        Args:
            mpc_numbers (int, list, tuple): MPC numbers to retrieve orbits for
            catalog (str): Catalog to retrieve data from

        Returns:
            None
        """
        if not isinstance(mpc_numbers, (list, tuple)):
            mpc_numbers = [mpc_numbers]
        if catalog not in _ASTEROID_CATALOG_NAMES.keys():
            raise ValueError(f'Catalog name invalid. Options are: {list(_ASTROMETRY_CATALOG_NAMES.keys())}')

        # Set up query to DB
        query_columns = ', '.join(_ASTEROID_TABLE_COLUMNS)
        query_mpc_numbers = ', '.join([str(mpc) for mpc in mpc_numbers])
        query_catalog = _ASTEROID_CATALOG_NAMES[catalog]
        query = f"SELECT {query_columns} FROM {query_catalog} WHERE number_mp IN ({query_mpc_numbers})"

        # Retrieve from astroquery
        try:
            job = Gaia.launch_job_async(query)
            table = job.get_results()

        except Exception as e:
            raise RuntimeError('Error while querying Gaia archives') from e

        table = table.to_pandas()
        table = table[table['epoch_state_vector'] != 0] # Filter entries with missing data
        table = self._to_tudat_format(table)
        table = self._add_orbital_elements(table)
        table = self._add_convenience_columns(table)

        self._table = table

    def retrieve_data_locally(self,
                              archive_file_path: Path | str,
                              mpc_numbers: int | list | tuple = None,):
        """
        Retrieve orbit data locally from .parquet file. Mirrors retrieve_data, but mpc_numbers can be empty, in which
        case all data is retrieved.

        Args:
            archive_file_path (Path | str): Path to the archive .parquet file.
            mpc_numbers (int, list, tuple): MPC numbers to retrieve orbits for. If none, select all objects in catalog

        Returns:
            None
        """
        if isinstance(mpc_numbers, int): # Convert single mpc to list
            mpc_numbers = [mpc_numbers]

        table = pd.read_parquet(archive_file_path)

        table = table[_ASTEROID_TABLE_COLUMNS]
        if mpc_numbers is not None: # Filter out MPC numbers
            table = table[table['number_mp'].isin(mpc_numbers)]
        table = table[table['epoch_state_vector'] != 0] # Filter empty datasets

        # Read literal lists and convert to arrays
        str_to_array = lambda s: np.array(ast.literal_eval(s))
        table['orbital_elements_var_covar_matrix'] = table['orbital_elements_var_covar_matrix'].apply(
            str_to_array
        )
        table['h_state_vector'] = table['h_state_vector'].apply(
            str_to_array
        )
        table['h_state_vector_var_covar_matrix'] = table['h_state_vector_var_covar_matrix'].apply(
            str_to_array
        )

        table = self._to_tudat_format(table)
        table = self._add_orbital_elements(table)
        table = self._add_convenience_columns(table)
        self._table = table

    @staticmethod
    def generate_parquet(archive_dir: Path | str,
                         dir_to_save: Path | str,):
        """
        Generate a parquet file from CSV files

        Args:
            archive_dir (Path | str): Path to the archive directory
            dir_to_save (Path | str): Path to the directory to save the parquet file

        Returns:
            None
        """
        print('Loading archive from CSV files...\n')

        table = pd.DataFrame()

        for i in range(20):
            # Get file path
            print(f'Loading file number {i}/19...')
            file_number = str(i) if i >= 10 else '0' + str(i)
            file_name = 'SsoSource_' + file_number + '.csv'
            csv_file = Path(archive_dir) / file_name
            # Concatenate into one table
            table_chunk = pd.read_csv(csv_file,
                                      comment='#')
            table = pd.concat([table, table_chunk], ignore_index=True)

        # Save into parquet
        table.to_parquet(Path(dir_to_save) / 'gaia_source_archive.parquet')

    def _to_tudat_format(self,
                         table):
        """
        Convert units and format to a tudat compatible format and apply scaling corrections to state vectors.

        Args:
            table (pd.DataFrame): Dataframe to format

        Returns:
            pd.DataFrame: Formatted dataframe
        """
        # Convert epoch to seconds since J2000
        jd_to_s = lambda jd: julian_day_to_seconds_since_epoch(jd + _J2010)
        table['epoch_state_vector'] = table['epoch_state_vector'].apply(jd_to_s)

        # Convert TCB to TDB epoch
        table['epoch_state_vector'] = table['epoch_state_vector'].apply(TCB_to_TDB)

        # Scaling factors
        length_conversion = ASTRONOMICAL_UNIT * _TIME_SCALE_CORRECTION * _STATE_SCALING_FACTOR # AU to SI
        velocity_conversion = ASTRONOMICAL_UNIT * _STATE_SCALING_FACTOR / 86400 # AU/day to SI

        # Upper triangle components to matrix
        to_matrix = lambda a: np.array(
            [
                [a[0], a[1], a[2], a[3], a[4], a[5]],
                [a[1], a[6], a[7], a[8], a[9], a[10]],
                [a[2], a[7], a[11], a[12], a[13], a[14]],
                [a[3], a[8], a[12], a[15], a[16], a[17]],
                [a[4], a[9], a[13], a[16], a[18], a[19]],
                [a[5], a[10], a[14], a[17], a[19], a[20]]
            ]
        )

        table['orbital_elements_var_covar_matrix'] = table['orbital_elements_var_covar_matrix'].apply(to_matrix)
        table['h_state_vector_var_covar_matrix'] = table['h_state_vector_var_covar_matrix'].apply(to_matrix)

        # Convert orbital element covariance
        a_scaling = np.identity(6)
        a_scaling[0,0] = length_conversion # SMA from AU -> m, rest are angles
        scale_cov_oe = lambda x: a_scaling @ x @ a_scaling
        table['orbital_elements_var_covar_matrix'] = table['orbital_elements_var_covar_matrix'].apply(scale_cov_oe)

        # Convert state vector elements
        scale_state = lambda x: np.concatenate((length_conversion * x[:3], velocity_conversion * x[3:]))
        table['h_state_vector'] = table['h_state_vector'].apply(scale_state)

        # Convert cartesian covariance
        state_vector_scaling = np.diag((length_conversion, length_conversion, length_conversion,
                                        velocity_conversion, velocity_conversion, velocity_conversion))
        scale_cov_cart = lambda x: state_vector_scaling @ x @ state_vector_scaling
        table['h_state_vector_var_covar_matrix'] = table['h_state_vector_var_covar_matrix'].apply(scale_cov_cart)

        return table

    def _add_orbital_elements(self, table):
        """
        Calculate and add orbital elements (heliocentric ecliptic)

        Args:
            table (pd.DataFrame): Table

        Returns:
            pd.DataFrame: Table with added orbital elements as [a,e,i,w,node,nu]
        """
        to_ecliptic = lambda x: np.concatenate((j2000_to_eclipj2000() @ x[:3], j2000_to_eclipj2000() @ x[3:]))

        table = table.assign(
            orbital_elements=lambda x: x['h_state_vector']
            .apply(to_ecliptic)
            .apply(cartesian_to_keplerian, gravitational_parameter=1.32712440042e20),
        )

        return table


    def _add_convenience_columns(self, table) -> pd.DataFrame:
        """
        Add some columns for to make accessing elements more convenient
        """
        # Important orbital elements
        orbital_elements = np.vstack(table.orbital_elements)
        a, e, i, _, _, _ = orbital_elements.T
        table['a'] = a
        table['e'] = e
        table['i'] = i
        table['q'] = a * (1 -e) # Perihelion distance
        table['Q'] = a * (1 + e) # Aphelion distance

        # Uncertainty in semi major axis
        cov_stacked = np.vstack(table['orbital_elements_var_covar_matrix'])
        table['sigma_a'] = np.sqrt(cov_stacked[0:-1:6, 0])

        return table





