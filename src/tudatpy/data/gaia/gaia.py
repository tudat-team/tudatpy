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
from scipy.linalg import block_diag
from scipy.constants import arcsec
from tudatpy.estimation.observable_models_setup.biases.photocenter_correction import photocenter_corrections_from_observations
from tudatpy.estimation.observable_models_setup.biases.light_deflection_correction import relativistic_light_deflection_from_observations
import copy
from pathlib import Path
from tudatpy.astro.element_conversion import j2000_to_eclipj2000, cartesian_to_keplerian
import ast
from collections.abc import Iterable

# Constants
_J2010 = 2455197.5  # Reference time J2010.0 in Julian days
_TIME_SCALE_CORRECTION = 1 - 1.550519768e-8 # See e.g. Klioner (2003)
_ASTROMETRY_TABLE_COLUMNS = ['epoch', 'ra', 'dec', 'number_mp', 'transit_id', 'position_angle_scan', 'denomination',
                   'ra_error_systematic', 'dec_error_systematic', 'ra_dec_correlation_systematic',
                   'ra_error_random', 'dec_error_random', 'ra_dec_correlation_random',
                   'x_gaia', 'y_gaia', 'z_gaia', 'vx_gaia', 'vy_gaia', 'vz_gaia',
                   'x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric',
                   'vx_gaia_geocentric', 'vy_gaia_geocentric', 'vz_gaia_geocentric']
_ASTEROID_TABLE_COLUMNS = [
    'number_mp', # MPC number
    'epoch_state_vector', # Epoch in TCB
    'orbital_elements_var_covar_matrix', # Cov. matrix of heliocentric orbital elements as (a, e, i, node, peri, M)
    'h_state_vector',
    'h_state_vector_var_covar_matrix'
]

_STATE_SCALING_FACTOR = 1.0000000051686297 # account for Gaia FPR state vector inconsistency (see Gaia Collaboration 2023)
_ASTROMETRY_CATALOG = 'gaiafpr.sso_observation' # Latest Gaia data release
_ASTEROID_CATALOG = 'gaiafpr.sso_source'


class GaiaAstrometry:
    """Handles retrieval and processing of Gaia solar-system astrometry."""
    def __init__(self,
                 observations_and_metadata: pd.DataFrame):
        """Create an empty GaiaAstrometry object."""
        self._table = observations_and_metadata
        self._corrected = False # Flag indicating if observation corrections have been applied


    @property
    def table(self) -> pd.DataFrame:
        """Copy of the observation table."""
        return self._table.copy()


    def table_for_mpc(self, mpc_number: int)-> pd.DataFrame:
        """Return a copy of the observations and metadata filtered for one asteroid by MPC number"""
        if mpc_number not in self.mpc_numbers:
            raise ValueError(f'Observations requested for {mpc_number}, but no observations were found')

        return self.table.query('number_mp == @mpc_number')


    @property
    def mpc_numbers(self) -> np.ndarray:
        """Array of asteroid MPC numbers that appear in the observation table."""
        return pd.unique(self._table['number_mp'])


    def copy(self) -> "GaiaAstrometry":
        """Get a copy of the GaiaAstrometry object.

        Returns
        -------
        GaiaAstrometry
            Deep copy of the current GaiaAstrometry object.
        """
        return copy.deepcopy(self)


    def copy_for_mpc(self,
                     mpc_numbers: int | Iterable) -> "GaiaAstrometry":
        """Return a copy of the query object for a selection of asteroids.

        Useful to create query objects for individual asteroids without calling
        retrieve_data every time.

        Parameters
        ----------
        mpc_numbers : int | list | tuple
            MPC numbers to keep in the copy.

        Returns
        -------
        GaiaAstrometry
            Copy containing only data for the given mpc_numbers.
        """
        if isinstance(mpc_numbers, Iterable):
            mpc_numbers = [mpc_numbers]

        if any(m not in self.mpc_numbers for m in mpc_numbers):
            raise ValueError('Requested MPC numbers are not in loaded observation table')

        new = self.copy()
        new._table = new._table.query('number_mp in @mpc_numbers')

        return new

    def get_observation_covariance_matrix(self,
                                          mpc_number: int) -> np.ndarray:
        """
        From the observation set, build the observation covariance matrix for a certain asteroid in the observations.

        Parameters
        ----------
        mpc_number : int
            MPC Number of asteroid to obtain covariance for

        Returns
        -------
        np.ndarray
            Observation covariance matrix
        """
        table = self.table_for_mpc(mpc_number)

        to_cov = lambda ra, dec, corr: np.array([[ra**2, ra*dec*corr], [ra*dec*corr, dec**2]])
        observation_covariance_matrix = []

        # Transit ID must be increasing with epoch to build a matrix consistent with the observations:
        if not table['transit_id'].is_monotonic_increasing:
            raise RuntimeError('Error while building observation covariance matrix: possible broken observation entries')

        for transit_id, transit_obs in table.groupby('transit_id', sort=False):

            transit_length = len(transit_obs)

            # Sigma's and correlation for current transit:
            uncertainty_random = transit_obs[['ra_error_random', 'dec_error_random', 'ra_dec_correlation_random']]
            uncertainty_sys = transit_obs[['ra_error_systematic', 'dec_error_systematic', 'ra_dec_correlation_systematic']]

            # Random and systematic components of covariance for current transit
            covariance_random = block_diag(*[to_cov(ra, dec, corr) for ra, dec, corr in uncertainty_random.to_numpy()])
            covariance_sys_sub = to_cov(*uncertainty_sys.iloc[0]) # Systematic component is constant over transit
            covariance_sys = np.tile(covariance_sys_sub, (transit_length, transit_length))

            observation_covariance_matrix.append(covariance_random + covariance_sys)

        # Combine individual transits into one block diagonal matrix, with blocks being the transits.
        observation_covariance_matrix = block_diag(*observation_covariance_matrix)

        return observation_covariance_matrix


    def to_observation_collection(self,
                 bodies: SystemOfBodies) -> observations.ObservationCollection:
        """Collect all Gaia observations into an ObservationCollection and apply the
        observation weights according to the Gaia weighting scheme. Any filtering or corrections must be done before
        constructing the observation collection.

        Parameters
        ----------
        bodies : SystemOfBodies
            The SystemOfBodies object. Must have the object 'Gaia' loaded along with its
            ephemeris, retrieved from get_gaia_ephemeris.

        Returns
        -------
        ObservationCollection
            Tudat ObservationCollection containing observations of all asteroids organized in SingleObservationSets
                by link-ends.
        """
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
            temp = self.table_for_mpc(mpc_number)
            observation_angles = temp.loc[:, ['ra', 'dec']].to_numpy()
            observation_times = temp['epoch'].to_numpy()

            observation_set = observations.create_single_observation_set(
                model_settings.angular_position_type,
                link_ends,
                observation_angles,
                observation_times,
                links.receiver
            )

            measurement_covariance_matrix = self.get_observation_covariance_matrix(mpc_number)
            weight_matrix = np.linalg.inv(measurement_covariance_matrix)

            observation_set.set_full_weight_matrix(weight_matrix)

            # Add SingleObservationSet to list
            observation_set_list.append(observation_set)

        observation_collection = observations.ObservationCollection(
            observation_set_list
        )

        return observation_collection

    def correct_observations(self,
                             bodies: SystemOfBodies,
                             diameters: dict = None,
                             light_deflection_bodies: list = ['Sun'],
                             correct_photocenter: bool = True) -> None:
        """Apply photocenter and/or light-deflection corrections to the observations.

        Corrections are made in-place, modifying the observations themselves. Applies the corrections to all asteroid
        data currently loaded.

        Parameters
        ----------
        bodies : SystemOfBodies
            Bodies object. Must have the Gaia body and its ephemeris loaded, and a
            reference ephemeris loaded for the asteroids over the complete timespan of
            observations.
        diameters : dict, optional
            Diameters of the asteroids as dict with keys the MPC numbers and values the diameters in meter.
        light_deflection_bodies : list, optional
            Names of perturber bodies involved in light bending, by default 'Sun'.
        correct_photocenter : bool, optional
            Whether photocenter-barycenter corrections should be applied, by default True.
        """
        # Input validation
        if self._corrected: # Check if this function has already been used for current object
            raise RuntimeError('correct_observations was called, but corrections were already applied on this instance')

        if correct_photocenter and diameters is None:
            raise ValueError('Diameters must be specified when photocenter correction is applied.')

        # Apply corrections to all asteroid data currently loaded
        for mpc_number in self.mpc_numbers:

            mask = self.table['number_mp'] == mpc_number
            observations_array = self.table.loc[mask, ['epoch', 'ra', 'dec']].to_numpy()

            # Photocenter
            if correct_photocenter:
                photocenter_corrections = photocenter_corrections_from_observations(
                    observations=observations_array,
                    diameter=diameters[mpc_number],
                    bodies=bodies,
                    body_name=str(mpc_number),
                    observer_body_name='Gaia'
                )
                self._table.loc[mask, ['ra', 'dec']] += photocenter_corrections

            # Apply light-bending corrections
            if light_deflection_bodies:

                light_bending_corrections = relativistic_light_deflection_from_observations(
                    observations=observations_array,
                    bodies=bodies,
                    body_name=str(mpc_number),
                    observer_body_name='Gaia',
                    perturbing_bodies_list=light_deflection_bodies,
                )
                self._table.loc[mask, ['ra', 'dec']] += light_bending_corrections

        # Wrap RA
        self._table['ra'] = (self._table['ra'] + np.pi) % (2 * np.pi) - np.pi

        # Set correction flag
        self._corrected = True

    @classmethod
    def load_from_astroquery(cls,
                             mpc_numbers: int | Iterable[int],
                             username: str | None = None,
                             password: str | None = None) -> "GaiaAstrometry":

        """Retrieve the astrometric observations through astroquery.

        Requires an internet connection. Login to the Gaia archive website is optional.
        Observations and metadata are stored on the table attribute.

        Parameters
        ----------
        mpc_numbers : int | Iterable[int]
            MPC numbers of the asteroids to retrieve data for.
        username : str, optional
            Username for the Gaia archives, by default None.
        password : str, optional
            Password for the Gaia archives, by default None.

        Returns
        -------
        GaiaAstrometry
            A GaiaAstrometry object with observations loaded
        """
        if not isinstance(mpc_numbers, Iterable):
            mpc_numbers = [mpc_numbers]

        # Define query to database
        query_mpc_numbers = ", ".join(str(mpc_number) for mpc_number in mpc_numbers)
        query_columns = ', '.join(_ASTROMETRY_TABLE_COLUMNS)

        login_provided = username and password

        if login_provided:
            Gaia.login(user=username, password=password)

        try:
            query = f"""
            SELECT {query_columns}
            FROM {_ASTROMETRY_CATALOG}
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

        # Convert units and reset index
        table = cls._convert_units(table)
        table = table.reset_index(drop=True)
        assert table['epoch'].is_monotonic_increasing # Sanity check for ordering by epoch

        return cls(table)


    @classmethod
    def load_from_local_archive(cls,
                                archive_file_path: str | Path,
                                mpc_numbers: int | Iterable[int],) -> "GaiaAstrometry":
        """Retrieve astrometry locally from a .parquet file (generated by generate_parquet).

        Mirrors load_from_astroquery. This method of loading data can be faster when loading a lot of data or when the Gaia
        archive connection is slow.

        Parameters
        ----------
        mpc_numbers : tuple[int] | list[int]
            List of asteroid MPC numbers to retrieve data for.
        archive_file_path : str | Path
            Path to the .parquet file.

        Returns
        -------
        GaiaAstrometry
            A GaiaAstrometry object with observations loaded
        """
        if not isinstance(mpc_numbers, Iterable):
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
        table = cls._convert_units(table)
        table = table.reset_index(drop=True)

        return cls(table)


    @staticmethod
    def generate_parquet(archive_dir: Path | str,
                         dir_to_save: Path | str) -> None:
        """Generate a .parquet file of the Gaia archive, to be used to retrieve data locally.

        Requires all .csv files to be stored in the same directory. CSV files can be
        downloaded from: https://cdn.gea.esac.esa.int/?prefix=Gaia/

        Parameters
        ----------
        archive_dir : Path | str
            Path to the archive directory where the CSV files are stored.
        dir_to_save : Path | str
            Directory where the .parquet file should be saved.
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

    @staticmethod
    def _convert_units(table: pd.DataFrame) -> pd.DataFrame:
        """Convert raw table values into a tudat-compatible format."""
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
        table.loc[:, vel_names] *= ASTRONOMICAL_UNIT / JULIAN_DAY

        return table


    def filter(self,
               epoch_start: float | DateTime,
               epoch_end: float | DateTime) -> None:
        """Filter observations (in-place) by epoch.

        Parameters
        ----------
        epoch_start : float | DateTime
            Lower bound for epoch.
        epoch_end : float | DateTime
            Upper bound for epoch.
        """
        if isinstance(epoch_start, DateTime):
            epoch_start = epoch_start.epoch()
        if isinstance(epoch_end, DateTime):
            epoch_end = epoch_end.epoch()

        # Find observations in time span
        epoch_filters = (self._table.epoch >= epoch_start) & (self._table.epoch <= epoch_end)
        filtered_table = self._table[epoch_filters]

        if filtered_table.empty:
            raise RuntimeError('No observations left after filtering by epochs')

        self._table = filtered_table


    def print_summary(self) -> None:
        """Print a summary of the loaded observations."""
        print('GAIA OBSERVATIONS SUMMARY: \n')
        print(f'Observations loaded for {len(self.mpc_numbers)} object(s):')
        print('MPC | DENOMINATION | NUMBER OF OBSERVATIONS | FIRST OBSERVATION TIME | LAST OBSERVATION TIME')

        for mpc_number in self.mpc_numbers:

            table = self.table_for_mpc(mpc_number)
            denom = table['denomination'].iloc[0]
            number_of_obs = len(table)
            dt_first = DateTime.from_epoch(table['epoch'].min())
            dt_final = DateTime.from_epoch(table['epoch'].max())
            first_epoch = f'{dt_first.year}/{dt_first.month}/{dt_first.day}'
            last_epoch = f'{dt_final.year}/{dt_final.month}/{dt_final.day}'

            asteroid_data = f'{mpc_number}   {denom}   {number_of_obs}   {first_epoch}   {last_epoch}\n'
            print(asteroid_data)


    def get_gaia_ephemeris(self,
                           geocentric: bool = True) -> ephemeris.EphemerisSettings:
        """Get tabulated ephemeris settings generated from the archived Gaia state vectors.

        Uses all reported state vectors in the loaded observation table.

        Parameters
        ----------
        geocentric : bool, optional
            If True, use the geocentric variant of the Gaia state vectors. Recommended
            because it reduces the dependency on the adopted planetary ephemeris model, by default True.

        Returns
        -------
        TabulatedEphemerisSettings
            Tabulated ephemeris settings of Gaia.
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
    """
    Class that retrieves Gaia-derived asteroid orbits and covariance data, found in the sso_source tables on the Gaia archives
    """
    def __init__(self,
                 asteroid_orbits_and_covariance: pd.DataFrame):
        """Construct a GaiaAsteroids object"""
        self._table = asteroid_orbits_and_covariance

    @property
    def table(self) -> pd.DataFrame:
        """Return a copy of the asteroid data table which contains orbit and covariance data."""
        return self._table.copy()

    def get_state_for_mpc(self,
                          mpc_number: int) -> tuple[float, np.ndarray]:
        """Retrieve the state vector from the table for a single MPC number. State is heliocentric in the J2000 frame

        Parameters
        ----------
        mpc_number : int
            MPC number of asteroid to retrieve state for

        Returns
        -------
        float
            Epoch of the state vector.
        np.ndarray
            State vector (heliocentric J2000).
        """
        asteroid_data = self.table.query('number_mp == @mpc_number').iloc[0]
        return asteroid_data.epoch_state_vector, asteroid_data.h_state_vector

    @classmethod
    def load_from_astroquery(cls,
                             mpc_numbers: int | Iterable,) -> "GaiaAsteroids":
        """Retrieve Gaia-derived asteroid data from astroquery. Requires an internet connection.

        Parameters
        ----------
        mpc_numbers : int | Iterable
            MPC numbers of asteroids to retrieve data for

        Returns
        -------
        GaiaAsteroids
            Instance with observations loaded from astroquery.
        """
        if not isinstance(mpc_numbers, Iterable):
            mpc_numbers = [mpc_numbers]

        # Set up query to DB
        query_columns = ', '.join(_ASTEROID_TABLE_COLUMNS)
        query_mpc_numbers = ', '.join([str(mpc) for mpc in mpc_numbers])
        query = f"SELECT {query_columns} FROM {_ASTEROID_CATALOG} WHERE number_mp IN ({query_mpc_numbers})"

        # Retrieve from astroquery
        try:
            job = Gaia.launch_job_async(query)
            table = job.get_results()

        except Exception as e:
            raise RuntimeError('Error while querying Gaia archives') from e

        table = table.to_pandas()
        if table.empty:
            raise LookupError(f'No asteroid data could be found for {mpc_numbers}')

        table = table[table['epoch_state_vector'] != 0] # Filter entries with missing data
        table = cls._convert_units(table)
        table = cls._add_orbital_elements(table)
        table = cls._add_orbit_class(table)

        return cls(table)

    @classmethod
    def load_from_local_archive(cls,
                                archive_file_path: Path | str,
                                mpc_numbers: int | Iterable | None = None) -> "GaiaAsteroids":
        """Retrieve orbit data locally from a .parquet file of the Gaia archive.

        Mirrors load_from_astroquery, but mpc_numbers can be empty, in which case the data is retrieved for all asteroids
        in the catalog.

        Parameters
        ----------
        archive_file_path : Path | str
            Path to the archive .parquet file.
        mpc_numbers : int | list | tuple, optional
            MPC numbers to retrieve orbits for. If None, select all objects in the
            catalog, by default None.
        """
        if not isinstance(mpc_numbers, Iterable):
            mpc_numbers = [mpc_numbers]

        filter = [('number_mp', 'in', mpc_numbers)] if mpc_numbers is not None else None
        table = pd.read_parquet(archive_file_path,filters=filter)
        if table.empty:
            raise LookupError(f'No asteroid data could be found for {mpc_numbers}')

        table = table[_ASTEROID_TABLE_COLUMNS]
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

        table = cls._convert_units(table)
        table = cls._add_orbital_elements(table)
        table = cls._add_orbit_class(table)

        return cls(table)

    @staticmethod
    def generate_parquet(archive_dir: Path | str,
                         dir_to_save: Path | str) -> None:
        """Generate a .parquet file from CSV files. This file can be passed to retrieve_data_locally. All CSV files
        labeled SsoSource_ must be in the archive directory. Files can be downloaded from
        https://cdn.gea.esac.esa.int/?prefix=Gaia/

        Parameters
        ----------
        archive_dir : Path | str
            Path to the archive directory.
        dir_to_save : Path | str
            Path to the directory to save the parquet file.
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

    @staticmethod
    def _convert_units(table: pd.DataFrame) -> pd.DataFrame:
        """Convert units and format to a tudat-compatible format and apply scaling
        corrections to state vectors."""
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

    @staticmethod
    def _add_orbital_elements(table: pd.DataFrame) -> pd.DataFrame:
        """Calculate and add orbital elements (heliocentric ecliptic)"""
        to_ecliptic = lambda x: np.concatenate((j2000_to_eclipj2000() @ x[:3], j2000_to_eclipj2000() @ x[3:]))

        table = table.assign(
            orbital_elements=lambda x: x['h_state_vector']
            .apply(to_ecliptic)
            .apply(cartesian_to_keplerian, gravitational_parameter=1.32712440042e20),
        )

        return table

    @staticmethod
    def _add_orbit_class(table: pd.DataFrame) -> pd.DataFrame:
        """Add some columns to make accessing elements more convenient."""
        # Important orbital elements
        orbital_elements = np.vstack(table.orbital_elements)
        a, e, i, _, _, _ = orbital_elements.T
        q = a * (1 - e) # Perihelion
        Q = a * (1 + e) # Aphelion

        au = ASTRONOMICAL_UNIT
        conditions = [
            Q < 0.983 * au,  # Atira
            (a < au) & (Q > 0.983 * au),  # Aten
            (a > au) & (q < 1.017 * au),  # Apollo
            (q > 1.017 * au) & (q < 1.3 * au),  # Amor
            (q > 1.3 * au) & (q < 1.666 * au) & (a < 3.2 * au),  # MCA
            (a < 2 * au) & (q > 1.666 * au),  # IMB
            (a > 2 * au) & (a < 3.2 * au) & (q > 1.666 * au),  # MB
            (a > 3.2 * au) & (a < 4.6 * au),  # OMB
            (a > 4.6 * au) & (a < 5.5 * au) & (e < 0.3),  # Trojan
            (a > 5.5 * au) & (a < 30.1 * au),  # Centaur
            a > 30.1 * au,  # TNO
        ]
        labels = ['Atira', 'Aten', 'Apollo', 'Amor', 'MCA',
                  'IMB', 'MB', 'OMB', 'Trojan', 'Centaur', 'TNO']
        table['orbit_class'] = np.select(conditions, labels, default='unknown')

        return table





