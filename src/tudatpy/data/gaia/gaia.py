"""
Retrieve Gaia FPR astrometry from the archives
"""
import numpy as np
import pandas as pd
from tudatpy.astro.time_representation import julian_day_to_seconds_since_epoch, TCB_to_TDB, DateTime
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import model_settings, links
from tudatpy.constants import ASTRONOMICAL_UNIT,JULIAN_DAY
from tudatpy.dynamics.environment_setup import ephemeris
from scipy.linalg import block_diag
from scipy.constants import arcsec
from tudatpy.estimation.observations.observation_corrections import (light_deflection_correction_angular_observations,
                                                                     photocenter_correction_angular_observations)
from warnings import warn
from tudatpy.interface import spice
import copy
from pathlib import Path
from tudatpy.astro.element_conversion import j2000_to_eclipj2000, cartesian_to_keplerian
import ast
from collections.abc import Iterable

# Constants
_J2010 = 2455197.5  # Reference time J2010.0 in Julian days
_TIME_SCALE_CORRECTION = 1 - 1.550519768e-8 # See e.g. Klioner (2003)
_STATE_SCALING_FACTOR = 1.0000000051686297 # account for Gaia FPR state vector inconsistency (see Gaia Collaboration 2023)
_ASTROMETRY_CATALOG = 'gaiafpr.sso_observation' # Latest Gaia data release as of sep. 2026
_ASTEROID_CATALOG = 'gaiafpr.sso_source'


def _warn_missing_entries(
        entries: np.ndarray | None,
        expected: np.ndarray | None,):
    """Output warning if there are missing entries in retrieved dataset"""
    if entries is not None and expected is not None:
        missing_entries = np.setxor1d(entries, expected)
        if len(missing_entries) > 0:
            warn(f'No data found for {missing_entries}')

def _as_iterable(data) -> Iterable:
    """Make sure scalars are iterable"""
    return [data] if np.isscalar(data) else data


def _generate_parquet(archive_dir: Path | str,
                      dir_to_save: Path | str,
                      csv_prefix: str,
                      parquet_name: str,
                      literal_array_columns: list | None = None) -> None:
    """Concatenate the 20 numbered Gaia archive CSV files into a single .parquet file."""
    print('Loading archive from CSV files\n')

    chunks = []

    for i in range(20):
        print(f'Loading file number {i}/19...')
        csv_file = Path(archive_dir) / f'{csv_prefix}_{i:02d}.csv'
        chunks.append(pd.read_csv(csv_file,
                                  comment='#'))

    table = pd.concat(chunks, ignore_index=True)
    del chunks

    # Literal lists need to be converted to Python lists
    str_to_list = lambda entry: ast.literal_eval(entry) if isinstance(entry, str) else None
    for column in literal_array_columns or []:
        table[column] = table[column].apply(str_to_list)

    # Sort for faster loading
    table = table.sort_values('number_mp', ignore_index=True)

    # Save to .parquet
    table.to_parquet(Path(dir_to_save) / parquet_name, row_group_size=250000)


def generate_astrometry_parquet(archive_dir: Path | str,
                                dir_to_save: Path | str) -> None:
    """
    Generate a .parquet file of the Gaia astrometry archive, to be used to retrieve data locally.

    The resulting file is saved as ``gaia_astrometry_archive.parquet`` and can be passed to
    :meth:`~tudatpy.data.gaia.GaiaAstrometry.load_from_local_archive`.

    Requires all SsoObservation_*.csv files to be stored in the same directory. CSV files can be
    downloaded from: https://cdn.gea.esac.esa.int/?prefix=Gaia/

    Parameters
    ----------
    archive_dir : Path | str
        Path to the archive directory where the CSV files are stored.
    dir_to_save : Path | str
        Directory where the .parquet file should be saved.
    """
    _generate_parquet(archive_dir,
                      dir_to_save,
                      csv_prefix='SsoObservation',
                      parquet_name='gaia_astrometry_archive.parquet')


def generate_asteroid_parquet(archive_dir: Path | str,
                              dir_to_save: Path | str) -> None:
    """
    Generate a .parquet file of the Gaia asteroid archive, to be used to retrieve data locally.

    The resulting file is saved as ``gaia_source_archive.parquet`` and can be passed to
    :meth:`~tudatpy.data.gaia.GaiaAsteroids.load_from_local_archive`.

    Requires all SsoSource_*.csv files to be stored in the same directory. CSV files can be
    downloaded from: https://cdn.gea.esac.esa.int/?prefix=Gaia/

    Parameters
    ----------
    archive_dir : Path | str
        Path to the archive directory where the CSV files are stored.
    dir_to_save : Path | str
        Directory where the .parquet file should be saved.
    """
    _generate_parquet(archive_dir,
                      dir_to_save,
                      csv_prefix='SsoSource',
                      parquet_name='gaia_source_archive.parquet',
                      literal_array_columns=['orbital_elements_var_covar_matrix',
                            'h_state_vector', 'h_state_vector_var_covar_matrix'])


class GaiaAstrometry:
    """The class acts as a container for all Gaia astrometric observations and its metadata. It takes care of
    retrieval, filtering/correcting of observations, applying weights and conversion to a tudat-compatible format."""
    def __init__(self,
                 observations_and_metadata: pd.DataFrame) -> None:
        """Create an empty GaiaAstrometry object.

        Usually the GaiaAstrometry class is instantiated via the classmethods:
        :meth:`~tudatpy.data.gaia.GaiaAstrometry.load_from_astroquery` or
        :meth:`~tudatpy.data.gaia.GaiaAstrometry.load_from_local_archive`.
        """
        self._table = observations_and_metadata
        self._corrected = False # Flag to prevent apply_corrections from being called twice


    @property
    def table(self) -> pd.DataFrame:
        """Read-only copy of the astrometry table."""
        return self._table.copy()


    @property
    def mpc_numbers(self) -> np.ndarray:
        """Array of asteroid MPC numbers that have data in the astrometry table."""
        return pd.unique(self._table['number_mp'])


    def copy(self) -> "GaiaAstrometry":
        """Get a copy of the ``GaiaAstrometry`` object.

        Returns
        -------
        GaiaAstrometry
            Deep copy of the current ``GaiaAstrometry`` object.
        """
        return copy.deepcopy(self)

    def _table_for_single_object(self,
                                 mpc_number: int)-> pd.DataFrame:
        """
        Retrieve the astrometry table for a single object, queried by MPC number.
        """
        return self._table[self._table['number_mp'] == mpc_number].reset_index(drop=True)


    def get_observation_covariance_matrix(self,
                                          mpc_number: int) -> np.ndarray:
        """
        From the retrieved astrometry and metadata, build the observation covariance matrix for one asteroid.

        This matrix is also passed to the :class:`~tudatpy.estimation.observations.ObservationCollection` upon calling
        :meth:`~tudatpy.data.gaia.GaiaAstrometry.to_observation_collection`.

        Parameters
        ----------
        mpc_number : int
            MPC Number of asteroid to obtain covariance for

        Returns
        -------
        np.ndarray
            Observation covariance matrix in the same order as the observations
        """
        if mpc_number not in self.mpc_numbers:
            raise RuntimeError(f'Requested observation covariance for {mpc_number}, but no observations were found.')
        table = self._table_for_single_object(mpc_number)

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
        """Collect all Gaia observations into an :class:`~tudatpy.estimation.observations.ObservationCollection` and apply the
        observation weights according to the Gaia weighting scheme. Any filtering or corrections must be done before
        constructing the observation collection.

        Parameters
        ----------
        bodies : SystemOfBodies
            The SystemOfBodies object. Must have the object 'Gaia' loaded along with its
            ephemeris, retrieved from :meth:`~tudatpy.data.gaia.GaiaAstrometry.get_gaia_ephemeris_settings`.

        Returns
        -------
        ObservationCollection
            Tudat ObservationCollection containing observations of all asteroids organized in SingleObservationSets
            by link-ends. Asteroids are named by their MPC number.
        """
        # Check if Gaia is in bodies
        if not bodies.does_body_exist('Gaia') or bodies.get('Gaia').ephemeris is None:
            raise ValueError('Gaia satellite and associated ephemeris must be loaded in SystemOfBodies')

        # Force the weight matrix to be completely symmetric (due to possible numerical error introduced in inversion)
        force_symmetric = lambda mat: (mat + mat.T) / 2

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
            temp = self._table_for_single_object(mpc_number)
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

            observation_set.set_full_weight_matrix(force_symmetric(weight_matrix))

            # Add SingleObservationSet to list
            observation_set_list.append(observation_set)

        observation_collection = observations.ObservationCollection(
            observation_set_list
        )

        return observation_collection

    def apply_corrections(self,
                          bodies: SystemOfBodies,
                          light_deflection_bodies: Iterable | None = ('Sun',),
                          photocenter_body_dimensions: dict | None = None) -> None:
        """
        Apply photocenter and/or light-deflection corrections to the observations.

        Apply photocenter and/or light-deflection corrections to the observations. Calls the functions
        :func:`~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations`
        and :func:`~tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations`.

        To apply corrections, the following must be loaded into the environment (SystemOfBodies):

        * Reference ephemeris of each of the objects for which astrometry exists on this instance, covering at least
          the entire span of observations (typically july 2014 - january 2020, for Gaia FPR and DR4);
        * Ephemeris of Gaia, retrieved from :meth:`~tudatpy.data.gaia.GaiaAstrometry.get_gaia_ephemeris_settings`;
        * Ephemeris of the Sun (for the photocenter correction);
        * Ephemeris of each body in ``light_deflection_bodies``;
        * Rotational model of each object for which the ellipsoid photocenter correction is used. The body-fixed frame
          must be aligned with the principal axis frame as specified in
          :func:`~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations`;

        After corrections are applied, they are stored on the table property for inspection and post-processing. Photocenter
        corrections are stored on ``photocenter_corr_ra`` and ``photocenter_corr_dec``, light-deflection corrections
        are stored on ``light_deflection_corr_ra`` and ``light_deflection_corr_dec``.

        Parameters
        ----------
        bodies : SystemOfBodies
            The SystemOfBodies object
        light_deflection_bodies : list, optional
            Names of perturber bodies involved in light bending. If None, no correction is applied. By default, 'Sun'.
        photocenter_body_dimensions : dict, optional
            Dimensions of bodies (in m) for which photocenter correction must be applied, keyed by MPC number. If the
            dimension is scalar, the body is assumed spherical with a radius equal to the input. If the dimension is an
            array of 3 floats as ``[a, b, c]``, it is assumed to be an ellipsoid with semi-axes a, b, and c. If None,
            no correction is applied.
        """
        if self._corrected:
            raise RuntimeError('correct_observations cannot be called more than once on the same instance')

        for mpc_number in self.mpc_numbers:
            observation_mask = self.table['number_mp'] == mpc_number
            observations_array = self.table.loc[observation_mask, ['epoch', 'ra', 'dec']].to_numpy()

            # Photocenter corrections
            if photocenter_body_dimensions is not None:
                photocenter_corrections = photocenter_correction_angular_observations(
                    observations=observations_array,
                    body_dimensions=photocenter_body_dimensions[mpc_number],
                    bodies=bodies,
                    body_name=str(mpc_number),
                    observer_body_name='Gaia',
                )
            else:
                photocenter_corrections = np.zeros((len(observations_array),2))

            # Light-bending
            if light_deflection_bodies is not None:
                light_bending_corrections = light_deflection_correction_angular_observations(
                    observations=observations_array,
                    bodies=bodies,
                    body_name=str(mpc_number),
                    observer_body_name='Gaia',
                    perturbing_bodies_list=light_deflection_bodies,
                )
            else:
                light_bending_corrections = np.zeros((len(observations_array),2))

            self._table.loc[observation_mask, ['ra', 'dec']] += (photocenter_corrections + light_bending_corrections)

            # Store corrections on table
            self._table.loc[
                observation_mask, ['photocenter_corr_ra', 'photocenter_corr_dec']] = photocenter_corrections
            self._table.loc[observation_mask, ['light_deflection_corr_ra', 'light_deflection_corr_dec']] = light_bending_corrections


        # Wrap RA
        self._table['ra'] = (self._table['ra'] + np.pi) % (2 * np.pi) - np.pi

        self._corrected = True

    @classmethod
    def load_from_astroquery(cls,
                             mpc_numbers: int | Iterable[int],
                             username: str | None = None,
                             password: str | None = None) -> "GaiaAstrometry":

        """
        Retrieve the astrometric observations through astroquery.

        Requires an internet connection. Login to the Gaia archive website is optional. Note this method of loading the
        data may be slow or unreliable. For loading large batches of data, it is recommended to use
        :meth:`~tudatpy.data.gaia.GaiaAstrometry.load_from_local_archive` instead.

        Observations and metadata are stored on the :meth:`~tudatpy.data.gaia.GaiaAstrometry.table` attribute.

        Parameters
        ----------
        mpc_numbers : int | list[int]
            List of asteroid MPC numbers to retrieve data for.
        username : str, optional
            Username for the Gaia archives, by default None.
        password : str, optional
            Password for the Gaia archives, by default None.

        Returns
        -------
        GaiaAstrometry
            A ``GaiaAstrometry`` object with observations loaded
        """
        from astroquery.gaia import Gaia # late import because it tends to be slow

        mpc_numbers = _as_iterable(mpc_numbers)

        # Define query to database
        query_mpc_numbers = ", ".join(str(mpc_number) for mpc_number in mpc_numbers)

        login_provided = username and password

        if login_provided:
            Gaia.login(user=username, password=password)

        try:
            query = f"""
            SELECT *
            FROM {_ASTROMETRY_CATALOG}
            WHERE number_mp IN ({query_mpc_numbers}) 
            """
            job = Gaia.launch_job_async(query)
            table = job.get_results()

        except Exception as err:
            raise RuntimeError(f'Error while retrieving astrometric observations: \n{err}') from err

        table = table.to_pandas() # Convert astropy table to dataframe
        if table.empty:
            raise LookupError(f'No observations found for mpc numbers {mpc_numbers}')
        _warn_missing_entries(np.unique(table['number_mp']), mpc_numbers)

        # Convert units and reset index
        table = cls._convert_units(table)
        table = table.sort_values(by=['number_mp', 'epoch']).reset_index(drop=True)

        return cls(table)


    @classmethod
    def load_from_local_archive(cls,
                                archive_file_path: str | Path,
                                mpc_numbers: int | Iterable[int],) -> "GaiaAstrometry":
        """
        Retrieve astrometry locally from a .parquet file (generated by :func:`~tudatpy.data.gaia.generate_astrometry_parquet`).

        Mirrors :meth:`~tudatpy.data.gaia.GaiaAstrometry.load_from_astroquery`. This method of loading data is typically
        much faster, at a small one-time cost of generating parquet files.

        Parameters
        ----------
        archive_file_path : str | Path
            Path to the .parquet file.
        mpc_numbers : int | list[int]
            List of asteroid MPC numbers to retrieve data for.

        Returns
        -------
        GaiaAstrometry
            A ``GaiaAstrometry`` object with observations loaded
        """
        mpc_numbers = _as_iterable(mpc_numbers)

        # Read from parquet
        filters = [
            ('number_mp', 'in', mpc_numbers),
        ]
        table = pd.read_parquet(archive_file_path, filters=filters)

        if table.empty:
            raise LookupError(f'No observations found for mpc numbers {mpc_numbers}')
        _warn_missing_entries(np.unique(table['number_mp']), mpc_numbers)

        # Convert units and save
        table = cls._convert_units(table)
        table = table.sort_values(by=['number_mp', 'epoch']).reset_index(drop=True)

        return cls(table)


    @staticmethod
    def _convert_units(table: pd.DataFrame) -> pd.DataFrame:
        """Convert raw table values into a tudat-compatible format."""
        # Convert epoch to seconds since J2000
        jd2010_to_epoch = lambda jd: julian_day_to_seconds_since_epoch(jd + _J2010)
        table['epoch'] = table['epoch'].apply(jd2010_to_epoch)

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
        position_labels = ['x_gaia', 'y_gaia', 'z_gaia', 'x_gaia_geocentric', 'y_gaia_geocentric', 'z_gaia_geocentric']
        table.loc[:, position_labels] *= ASTRONOMICAL_UNIT * _TIME_SCALE_CORRECTION
        velocity_labels = ['vx_gaia', 'vy_gaia', 'vz_gaia', 'vx_gaia_geocentric', 'vy_gaia_geocentric',
                     'vz_gaia_geocentric']
        table.loc[:, velocity_labels] *= ASTRONOMICAL_UNIT / JULIAN_DAY

        return table


    def apply_filters(self,
                      exclude_poor_observations: bool = True,
                      mpc_numbers: int | Iterable[int] | None = None,
                      epoch_start: float | DateTime | None = None,
                      epoch_end: float | DateTime | None = None, ) -> None:
        """
        Apply filters to the set of observations in-place. This modifies the dataset.

        Parameters
        ----------
        exclude_poor_observations : bool
            Exclude observations which have ``astrometric_outcome_ccd != 1`` or ``astrometric_outcome_transit != 1``.
            These observations may be affected by one of several issues and can be unreliable. See Gaia documentation
            for more information.
        mpc_numbers : int | Iterable[int], optional
            Retain only objects with these MPC numbers
        epoch_start : float | DateTime, optional
            Remove observations before this epoch
        epoch_end : float | DateTime, optional
            Remove observations later than this epoch.
        """
        # Filter observations
        if exclude_poor_observations:
            poor_observation_filter = (
                    (self._table['astrometric_outcome_ccd'] != 1) |
                    (self._table['astrometric_outcome_transit'] != 1))
            self._table = self._table[~poor_observation_filter]

        # Filter by epoch
        if isinstance(epoch_start, DateTime):
            epoch_start = epoch_start.epoch()
        if isinstance(epoch_end, DateTime):
            epoch_end = epoch_end.epoch()

        if epoch_start is not None:
            epoch_start_filter = (self._table['epoch'] >= epoch_start)
            self._table = self._table[epoch_start_filter]
        if epoch_end is not None:
            epoch_end_filter = (self._table['epoch'] <= epoch_end)
            self._table = self._table[epoch_end_filter]

        # Filter MPC numbers
        if mpc_numbers is not None:
            mpc_numbers = _as_iterable(mpc_numbers)
            mpc_number_filter = self._table['number_mp'].isin(mpc_numbers)
            self._table = self._table[mpc_number_filter]

        if self._table.empty:
            raise RuntimeError('No observations left after applying filters')


    def print_summary(self) -> None:
        """Print a summary of the loaded observations."""
        print('GAIA OBSERVATIONS SUMMARY: \n')
        print(f'Observations loaded for {len(self.mpc_numbers)} object(s):')
        print('MPC | DENOMINATION | NUMBER OF OBSERVATIONS | FIRST OBSERVATION TIME | LAST OBSERVATION TIME')

        for mpc_number in self.mpc_numbers:

            table = self._table_for_single_object(mpc_number)
            denom = table['denomination'].iloc[0]
            number_of_obs = len(table)
            dt_first = DateTime.from_epoch(table['epoch'].min())
            dt_final = DateTime.from_epoch(table['epoch'].max())
            first_epoch = f'{dt_first.year}/{dt_first.month}/{dt_first.day}'
            last_epoch = f'{dt_final.year}/{dt_final.month}/{dt_final.day}'

            asteroid_data = f'{mpc_number}   {denom}   {number_of_obs}   {first_epoch}   {last_epoch}\n'
            print(asteroid_data)


    def get_gaia_ephemeris_settings(self,
                                    geocentric: bool = True) -> ephemeris.EphemerisSettings:
        """
        Get tabulated ephemeris settings generated from the archived Gaia state vectors.

        Uses all reported state vectors in the loaded astrometry table.

        Parameters
        ----------
        geocentric : bool, optional
            If True, use the geocentric variant of the Gaia state vectors, if false, use the barycentric variant.
            It is recommended to use the geocentric variant, because it reduces the dependency on the adopted planetary
            ephemeris model, by default True.

        Returns
        -------
        EphemerisSettings
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
    Class that acts as a container for Gaia-derived asteroid data, such as state vectors and covariance. It takes care of
    the retrieval and conversion to a tudat-compatible format.
    """
    def __init__(self,
                 asteroid_orbits_and_covariance: pd.DataFrame) -> None:
        """Create an empty GaiaAsteroids object.

        Usually, this object is created by the classmethods
        :meth:`~tudatpy.data.gaia.GaiaAsteroids.load_from_astroquery` or
        :meth:`~tudatpy.data.gaia.GaiaAsteroids.load_from_local_archive`.
        """
        self._table = asteroid_orbits_and_covariance

    @property
    def table(self) -> pd.DataFrame:
        """Read-only copy of the asteroid data table which contains orbit and covariance data."""
        return self._table.copy()


    def copy(self) -> 'GaiaAsteroids':
        """Return a copy of the object itself"""
        return copy.deepcopy(self)


    @property
    def mpc_numbers(self)->np.ndarray:
        """List of all asteroid MPC numbers for which data is retrieved."""
        return pd.unique(self._table['number_mp'])


    def get_state_for_object(self,
                             mpc_number: int,
                             frame_origin: str = 'Sun',
                             frame_orientation: str = 'J2000') -> tuple[float, np.ndarray]:
        """
        Retrieve state vector for an object queried by MPC number. An appropriate spice kernel must be loaded to
        translate the ``frame_origin``.

        Parameters
        ----------
        mpc_number : int
            MPC number of asteroid to retrieve state for
        frame_origin : str
            Origin of the state vector, by default 'Sun'
        frame_orientation : str
            Orientation of the state vector, by default 'J2000'

        Returns
        -------
        float
            Epoch of the state vector in seconds since J2000.
        np.ndarray
            State vector in SI units.
        """
        if mpc_number not in self.mpc_numbers:
            raise LookupError(f'State requested for {mpc_number}, but no asteroid data found in fetched archive')

        asteroid_row = self.table.query('number_mp == @mpc_number').iloc[0]
        epoch, state = asteroid_row['epoch_state_vector'], asteroid_row['h_state_vector'] # Heliocentric J2000

        # Translate origin
        if frame_origin != 'Sun':
            origin_state = spice.get_body_cartesian_state_at_epoch(
                target_body_name = 'Sun',
                observer_body_name = frame_origin,
                reference_frame_name = 'J2000',
                aberration_corrections = 'NONE',
                ephemeris_time = epoch
            )
            state = state + origin_state # Not in-place, to avoid mutating the stored state vector

        # Rotate frame
        if frame_orientation == 'J2000':
            pass
        elif frame_orientation == 'ECLIPJ2000':
            to_ecliptic = lambda x: np.concatenate((j2000_to_eclipj2000() @ x[:3], j2000_to_eclipj2000() @ x[3:]))
            state = to_ecliptic(state)
        else:
            raise ValueError('frame_orientation must be J2000 or ECLIPJ2000')

        return epoch, state

    @classmethod
    def load_from_astroquery(cls,
                             mpc_numbers: int | Iterable | None = None) -> "GaiaAsteroids":
        """
        Retrieve Gaia-derived asteroid data from astroquery. Requires an internet connection. Note this method of
        loading the data may be slow or unreliable. For loading large batches of data, it is recommended to use
        :meth:`~tudatpy.data.gaia.GaiaAsteroids.load_from_local_archive` instead.

        Parameters
        ----------
        mpc_numbers : int | list[int], optional
            MPC numbers of asteroids to retrieve data for. If None, the complete archive of Gaia-derived orbits
            is retrieved (~150,000 objects for Gaia FPR). This may take several minutes.

        Returns
        -------
        GaiaAsteroids
            Instance with data loaded from astroquery.
        """
        from astroquery.gaia import Gaia  # late import because it tends to be slow

        mpc_numbers = _as_iterable(mpc_numbers)

        # Set up query to DB
        query = f"SELECT * FROM {_ASTEROID_CATALOG}"
        if mpc_numbers is not None:
            query_mpc_numbers = ', '.join([str(mpc) for mpc in mpc_numbers])
            query += f" WHERE number_mp IN ({query_mpc_numbers})"

        # Retrieve from astroquery
        try:
            job = Gaia.launch_job_async(query)
            table = job.get_results()

        except Exception as e:
            raise RuntimeError('Error while querying Gaia archives') from e

        table = table.to_pandas()
        if table.empty:
            raise LookupError(f'No asteroid data could be found for {mpc_numbers}')
        _warn_missing_entries(np.unique(table['number_mp']), mpc_numbers)

        table = table[table['epoch_state_vector'] != 0] # Filter entries with missing data
        table = cls._convert_units(table)
        table = cls._add_orbital_elements(table)
        table = cls._add_orbit_class(table)
        table = table.sort_values(by='number_mp').reset_index(drop=True)

        return cls(table)

    @classmethod
    def load_from_local_archive(cls,
                                archive_file_path: Path | str,
                                mpc_numbers: int | Iterable | None = None) -> "GaiaAsteroids":
        """
        Retrieve orbit data locally from a .parquet file of the Gaia archive.

        Mirrors :meth:`~tudatpy.data.gaia.GaiaAsteroids.load_from_astroquery`. This method of loading data is
        typically much faster, at a small one-time cost of generating the .parquet file.

        Orbit and covariance data are stored on the :meth:`~tudatpy.data.gaia.GaiaAsteroids.table` attribute.

        Parameters
        ----------
        archive_file_path : Path | str
            Path to the archive .parquet file.
        mpc_numbers : int | list[int], optional
            MPC numbers of asteroids to retrieve data for. If None, the complete archive of Gaia-derived orbits
            is retrieved (~150,000 objects for Gaia FPR).

        Returns
        -------
        GaiaAsteroids
            Instance with data loaded from a .parquet file of the Gaia archive.
        """
        mpc_numbers = _as_iterable(mpc_numbers)

        filter = [('number_mp', 'in', mpc_numbers)] if mpc_numbers is not None else None
        table = pd.read_parquet(archive_file_path, filters=filter)
        if table.empty:
            raise LookupError(f'No asteroid data could be found for {mpc_numbers}')
        _warn_missing_entries(np.unique(table['number_mp']), mpc_numbers)

        table = table[table['epoch_state_vector'] != 0] # Filter empty datasets

        table = cls._convert_units(table)
        table = cls._add_orbital_elements(table)
        table = cls._add_orbit_class(table)
        table = table.sort_values(by='number_mp').reset_index(drop=True)

        return cls(table)

    @staticmethod
    def _convert_units(table: pd.DataFrame) -> pd.DataFrame:
        """Convert units and format to a tudat-compatible format and apply scaling
        corrections to state vectors."""
        # Convert epoch to seconds since J2000
        jd2010_to_epoch = lambda jd: julian_day_to_seconds_since_epoch(jd + _J2010)
        table['epoch_state_vector'] = table['epoch_state_vector'].apply(jd2010_to_epoch)

        # Convert TCB to TDB epoch
        table['epoch_state_vector'] = table['epoch_state_vector'].apply(TCB_to_TDB)

        # Scaling factors
        length_conversion = ASTRONOMICAL_UNIT * _TIME_SCALE_CORRECTION * _STATE_SCALING_FACTOR # AU to SI
        velocity_conversion = ASTRONOMICAL_UNIT * _STATE_SCALING_FACTOR / JULIAN_DAY # AU/day to SI

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
        scale_covariance_orbit_elements = lambda x: a_scaling @ x @ a_scaling
        table['orbital_elements_var_covar_matrix'] = table['orbital_elements_var_covar_matrix'].apply(
            scale_covariance_orbit_elements)

        # Convert state vector elements
        scale_state = lambda x: np.concatenate((length_conversion * x[:3], velocity_conversion * x[3:]))
        table['h_state_vector'] = table['h_state_vector'].apply(scale_state)

        # Convert cartesian covariance
        state_vector_scaling = np.diag((length_conversion, length_conversion, length_conversion,
                                        velocity_conversion, velocity_conversion, velocity_conversion))
        scale_covariance_cartesian = lambda x: state_vector_scaling @ x @ state_vector_scaling
        table['h_state_vector_var_covar_matrix'] = table['h_state_vector_var_covar_matrix'].apply(scale_covariance_cartesian)

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
        """Add orbit class of each asteroid according to JPL SBDB convention."""
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





