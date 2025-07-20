import typing
import astropy
import datetime
import numpy as np
import pandas as pd
from _typeshed import Incomplete
from tudatpy.dynamics import environment as environment, environment_setup as environment_setup
from tudatpy.dynamics.environment_setup import add_gravity_field_model as add_gravity_field_model
from tudatpy.dynamics.environment_setup.gravity_field import central_sbdb as central_sbdb
from tudatpy.estimation import observations as observations
from tudatpy.estimation.observable_models_setup import links as links, model_settings as model_settings
BIAS_LOWRES_FILE: Incomplete
DEFAULT_CATALOG_FLAGS: Incomplete

def load_bias_file(filepath: str, Nside: None | int=None, catalog_flags: list=...) -> tuple[pd.DataFrame, int]:
    """Loads a healpix star catalog debias file and processes it into a dataframe. Automatically retrieves NSIDE parameter.
    
    Parameters
    ----------
    filepath : str
        Filepath of debias file.
    Nside : Union[None, int], optional
        NSIDE value, to be left None in most cases as this is retrieved automatically by the function, by default None
    catalog_flags : Union[None, list], optional
        list of catalog flags, should be left default in most cases, by default None
    
    Returns
    -------
    Tuple[pd.DataFrame, int]
        Dataframe with biases in multiindex format ((Npix x Ncat) x Nvals), the numpix value
    
    Raises
    ------
    ValueError
        If NSIDE cannot be retrieved automatically."""

def get_biases_EFCC18(RA: float | np.ndarray | list, DEC: float | np.ndarray | list, epochJ2000secondsTDB: float | np.ndarray | list, catalog: str | np.ndarray | list, bias_file: str | None=..., Nside: int | None=None, catalog_flags: list[str]=...) -> tuple[np.ndarray, np.ndarray]:
    """Calculate and return star catalog bias values as described in:
    "Star catalog position and proper motion corrections in asteroid astrometry II: The Gaia era" by Eggl et al. (2018).
    Uses the regular bias set by default. A high res version of the bias map can be retrieved from the paper.
    This can then be selected using the bias_file paramater.
    
    Parameters
    ----------
    RA : Union[float, np.ndarray, list]
        Right Ascension value in radians
    DEC : Union[float, np.ndarray, list]
        Declination value in radians
    epochJ2000secondsTDB : Union[float, np.ndarray, list]
        Time in seconds since J2000 TDB.
    catalog : Union[str, np.ndarray, list]
        Star Catalog code as described by MPC: https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
    bias_file : Tuple[str, None], optional
        Optional bias file location, used to load in alternative debias coefficients. By default coefficients are retrieved from Tudat resources, by default None
    Nside : Tuple[int, None], optional
        Optional Nside value, should be left None in most cases, by default None
    catalog_flags : Tuple[List[str], None], optional
        List of catalog values to use, should be left None in most cases, by default None
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Right Ascencion Corrections, Declination corrections
    
    Raises
    ------
    ValueError
        If all mandatory inputs are not matching in size."""

def get_weights_VFCC17(MPC_codes: pd.Series | list | np.ndarray | None=None, epochUTC: pd.Series | list | np.ndarray | None=None, observation_type: pd.Series | list | np.ndarray | None=None, observatory: pd.Series | list | np.ndarray | None=None, star_catalog: pd.Series | list | np.ndarray | None=None, mpc_table: pd.DataFrame | None=None, return_full_table: bool=False) -> np.ndarray | pd.DataFrame:
    """Retrieves observation weights using the weighting scheme presented in
    "Statistical analysis of astrometric errors for the most productive
    asteroid surveys" by Veres et al. (2017). Input may be provided using
    either a full MPC table (e.g. from BatchMPC) or using the individual
    variables.
    
    Observation types: "x", "X", "V", "v", "W", "w", "R", "r", "Q", "q", "O",
    are not described by the paper and receive a placeholder weight of 1/100
    if provided.
    
    Parameters
    ----------
    MPC_codes : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the MPC target codes, e.g. 433 for Eros. Size must match
        other iterables, by default None
    epochUTC : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with UTC times. Size must match other iterables, by default None
    observation_type : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the observation types in MPC format.
        See the NOTE2 section of the MPC format description for the exact encoding:
        https://minorplanetcenter.net/iau/info/OpticalObs.html.
        Size must match other iterables, by default None
    observatory : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the MPC target codes, e.g. 433 for Eros.
        Size must match other iterables, by default None
    star_catalog : Union[pd.Series, list, np.ndarray, None], optional
        Iterable with the star catalog codes.
        See the MPC catalog codes page for the exact encoding:
        https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html.
        Size must match other iterables, by default None
    mpc_table : Union[pd.DataFrame, None], optional
        Table retrieved by calling the mpc.BatchMPC.table property.
        Set None when using iterable input.
        Set others None when using table, by default None
    return_full_table : bool, optional
        Return the table with all intermediate calculations if True,
        return a numpy array if False, by default False
    
    Returns
    -------
    np.ndarray
        If `return_full_table` is False, numpy array with weights with same size as input.
    pd.DataFrame
        If `return_full_table` is True, pandas table with all intermediate calculations.
    
    Raises
    ------
    ValueError
        MPC_codes, epochUTC, observation_type, observatory and star_catalog must all
        be not None and the same size. mpc_table must be None.
        If table input is used, the remaining input parameters must be done."""

class BatchMPC:
    """This class provides an interface between observations
    in the Minor Planet Centre database and Tudat.
    
    Notes
    ----------
    Currently, transformations between reference frames are not implemented.
    As such observations are only returned in the J2000 Frame.
    
    Examples
    ----------
    Basic sequence of usage:
    
    Initialise and retrieve data:
    
    >>> MPCcodes = [1, 4] # Ceres and Vesta
    >>> batch = BatchMPC()
    >>> batch.get_observations(MPCcodes)
    >>> batch.filter(epoch_start=datetime.datetime(2016, 1, 1))
    
    Transform to Tudat format:
    >>> ...
    >>> bodies = environment_setup.create_system_of_bodies(body_settings)
    >>> observation_collection = batch.to_tudat(bodies=bodies, included_satellites=None)"""

    def __init__(self) -> None:
        """Create an empty MPC batch."""

    def __copy__(self):
        ...

    def copy(self) -> BatchMPC:
        """Create a copy of the batch, equivalent to copy.copy(batchMPC())

        Returns
        -------
        BatchMPC
            Copy of batch.
        """

    @property
    def table(self) -> pd.DataFrame:
        """Pandas dataframe with observation data"""

    @property
    def observatories(self) -> list[str]:
        """List of observatories in batch"""

    @property
    def space_telescopes(self) -> list[str]:
        """List of satellite_observatories in batch"""

    @property
    def MPC_objects(self) -> list[str]:
        """List of MPC objects"""

    @property
    def size(self) -> int:
        """Number of observations in batch"""

    @property
    def bands(self) -> list[str]:
        """List of bands in batch"""

    @property
    def epoch_start(self) -> float:
        """Epoch of oldest observation in batch in seconds since J2000 TDB"""

    @property
    def epoch_end(self) -> float:
        """Epoch of latest observation in batch in seconds since J2000 TDB"""

    @property
    def bodies_created(self) -> dict:
        """Dictionary with the bodies created by to_tudat and details."""

    def __len__(self) -> int:
        ...

    def __add__(self, other):
        ...

    def get_observations(self, MPCcodes: list[str | int], drop_misc_observations: bool=True) -> None:
        """Retrieve all observations for a set of MPC listed objeccts.
        This method uses astroquery to retrieve the observations from the MPC.
        An internet connection is required, observations are cached for faster subsequent retrieval.
        Removes duplicate and irrelevant observation data by default (see `drop_misc_observations`).

        Parameters
        ----------
        MPCcodes : List[int]
            List of integer MPC object codes for minor planets or and comets.
        drop_misc_observations : bool, default True
            Drops observations made by method: radar and offset (natural satellites).
            Drops observations made by roaming observers.
            Drops duplicate listings to denote first observation.
        """

    def from_astropy(self, table: astropy.table.QTable, in_degrees: bool=True, frame: str='J2000'):
        """Manually input an astropy table with observations into the batch. 
        Usefull when manual filtering from astroquery is required
        Table must contain the following columns:
        number - MPC code of the minor planet\\\n        epoch - in julian days\\\n        RA - right ascension in either degrees or radians (degrees is default)\\\n        DEC - declination in either degrees or radians (degrees is default)\\\n        band - band of the observation (currently unused)\\\n        observatory - MPC code of the observatory

        Parameters
        ----------
        table : astropy.table.QTable
            Astropy table with the observations
        in_degrees : bool, optional
            if True RA and DEC are assumed in degrees, else radians, by default True
        frame : str, optional
            Reference frame. Please note that only J2000 is currently supported
            , by default "J2000"
        """

    def from_pandas(self, table: pd.DataFrame, in_degrees: bool=True, frame: str='J2000'):
        """Manually input an pandas dataframe with observations into the batch. 
        Usefull when manual filtering from astroquery is required
        Table must contain the following columns:
        number - MPC code of the minor planet\\\n        epoch - in julian days\\\n        RA - right ascension in either degrees or radians (degrees is default)\\\n        DEC - declination in either degrees or radians (degrees is default)\\\n        band - band of the observation (currently unused)\\\n        observatory - MPC code of the observatory

        Parameters
        ----------
        table : astropy.table.QTable
            Astropy table with the observations
        in_degrees : bool, optional
            if True RA and DEC are assumed in degrees, else radians, by default True
        frame : str, optional
            Reference frame. Please not that only J2000 is currently supported
            , by default "J2000"
        """

    def set_weights(self, weights: list | np.ndarray | pd.Series):
        """Manually set weights per observation. Weights are passed to
        observation collection when `.to_tudat()` is called. Set the
        `apply_weights_VFCC17` parameter in `.to_tudat()` to `False` to avoid
        overwriting. The order of the weights should match the order found in
        the `.table` parameter.

        Parameters
        ----------
        weights : Union[list, np.ndarray, pd.Series]
            Iterable with weights per observation.

        Raises
        ------
        ValueError
            If the size of weights does not match the number of observations in
            the batch table.
        """

    def filter(self, bands: list[str] | None=None, catalogs: list[str] | None=None, observation_types: list[str] | None=None, observatories: list[str] | None=None, observatories_exclude: list[str] | None=None, epoch_start: float | datetime.datetime | None=None, epoch_end: float | datetime.datetime | None=None, in_place: bool=True) -> None | BatchMPC:
        """Filter out observations from the batch.

        Parameters
        ----------
        bands : Union[List[str], None], optional
            List of observation bands to keep in the batch, by default None
        catalogs : Union[List[str], None], optional
            List of star catalogs to keep in the batch. See 
            https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
            for the exact encodings used by MPC, by default None
        observation_types : Union[List[str], None], optional
            List of observation types to keep in the batch, e.g. CCD,
            Photographic etc. See the Note2 section of the MPC format description:
            https://www.minorplanetcenter.net/iau/info/OpticalObs.html
            for the exact encoding, by default None
        observatories : Union[List[str], None], optional
            List of observatories to keep in the batch, by default None
        observatories_exclude : Union[List[str], None], optional
            List of observatories to remove from the batch, by default None
        epoch_start : Union[float, datetime.datetime, None], optional
            Start date to include observations from, can be in python datetime in utc                 or the more conventional tudat seconds since j2000 in TDB if float,                     by default None
        epoch_end : Union[float, datetime.datetime, None], optional
            Final date to include observations from, can be in python datetime in utc                 or the more conventional tudat seconds since j2000 in TDB if float,                     by default None
        in_place : bool, optional
            If true, modify the current batch object.                  If false returns a new object that is                      filtered, currect batch remains unmodified, by default True

        Raises
        ------
        ValueError
            Is raised if bands, observatories, or observatories_exclude are not list or 
            None.
        ValueError
            Is raised if both observations_exclude and observatories are not None.

        Returns
        -------
        None or BatchMPC
            Returns a new instance of BatchMPC that is filtered.
        """

    def to_tudat(self, bodies: environment.SystemOfBodies, included_satellites: dict[str, str] | None, station_body: str='Earth', add_sbdb_gravity_model: bool=False, apply_weights_VFCC17: bool=True, apply_star_catalog_debias: bool=True, debias_kwargs: dict=...) -> observations.ObservationCollection:
        """Converts the observations in the batch into a Tudat compatible format and
          sets up the relevant Tudat infrastructure to support estimation.
        This method does the following:\\\n            1. (By Default) Applies star catalog debiasing and estimation
            weight calculation.\\\n            2. Creates an empty body for each minor planet with their MPC code as a 
            name.\\\n            3. Adds this body to the system of bodies inputted to the method.\\\n            4. Retrieves the global position of the terrestrial observatories in 
            the batch and adds these stations to the Tudat environment.\\\n            5. Creates link definitions between each unique terrestrial 
            observatory/ minor planet combination in the batch.\\\n            6. (Optionally) creates a link definition between each 
            space telescope / minor planet combination in the batch. 
            This requires an addional input.\\\n            7. Creates a `SingleObservationSet` object for each unique link that 
            includes all observations for that link.\\\n            8. (By Default) Add the relevant weights to the `SingleObservationSet`
            per observation.\\\n            8. Returns the observations


        Parameters
        ----------
        bodies : environment.SystemOfBodies
            SystemOfBodies containing at least the earth to allow for the placement of 
            terrestrial telescopes
        included_satellites : Union[Dict[str, str], None], optional
            A dictionary that links the name of a space telescope used by the user with 
            the observatory code in MPC. Used when utilising observations from space 
            telescopes like TESS and WISE. The keys should be the MPC observatory 
            codes. The values should be the bodies\' in the user\'s code. The relevant 
            observatory code can be retrieved using 
            the .observatories_table() method, by default None
        station_body : str, optional
            Body to attach ground stations to. Does not need to be changed unless the 
            `Earth` body has been renamed, by default "Earth"
        station_body : bool, optional
            Adds a central_sbdb gravity model to the object, generated using JPL\'s small body database. 
            This option is only available for a limited number of bodies and raises an error if unavailable. 
            See tudatpy.dynamics.environment_setup.gravity_field.central_sbdb for more info.
            Enabled if True, by default False
        apply_weights_VFCC17 : bool, optional
            Applies the weighting scheme as described in: "Statistical analysis of astrometric errors
            for the most productive asteroid surveys" by Veres et al. (2017). Overwrites custom weights
            set through the `.set_weights()` method if True, by default True
        apply_star_catalog_debias : bool, optional
            Applies star catalog debiasing as described in: "Star catalog position and proper motion corrections 
            in asteroid astrometry II: The Gaia era" by Eggl et al. (2018), by default True
        debias_kwargs : dict, optional
            Additional options when applying star catalog debiasing. A different debias file 
            can be set here. Options are set as kwargs using a dictionary, see data._biases.get_biases_EFCC18() 
            for more info, by default dict()

        Returns
        -------
        observations.ObservationCollection
            ObservationCollection with the observations found in the batch

        Examples
        ----------
        Using Space Telescopes

        Create dictionary to link name in tudat with mpc code in batch:

        >> sats_dict = {"C57":"TESS"}
        >> obs = batch1.to_tudat(bodies, included_satellites=sats_dict)


        """

    def plot_observations_temporal(self, objects: list[str] | None=None, figsize: tuple[float]=(9.0, 6.0)):
        """Generates a matplotlib figure with the declination and right ascension
        over time.


        Parameters
        ----------
        objects : Union[List[str], None], optional
            List of specific MPC objects in batch to plot, None to plot all
            , by default None
        projection : str, optional
            projection of the figure options are: \'aitoff\', \'hammer\',
            \'lambert\' and \'mollweide\', by default "aitoff"
        figsize : Tuple[float], optional
            size of the matplotlib figure, by default (15, 7)

        Returns
        -------
        Matplotlib figure
            Matplotlib figure with the observations
        """

    def plot_observations_sky(self, objects: list[str] | None=None, projection: str | None=None, figsize: tuple[float]=(14.0, 7.0)):
        """Generates a matplotlib figure with the observations\'
        right ascension and declination over time.

        Parameters
        ----------
        objects : Union[List[str], None], optional
            List of specific MPC objects in batch to plot, None to plot all
            , by default None
        projection : str, optional
            projection of the figure options are: \'aitoff\', \'hammer\',
            \'lambert\' and \'mollweide\', by default "aitoff"
        figsize : Tuple[float], optional
            size of the matplotlib figure, by default (15, 7)

        Returns
        -------
        Matplotlib figure
            Matplotlib figure with the observations
        """

    def summary(self) -> None:
        """Produce a quick summary of the content of the batch"""

    def observatories_table(self, only_in_batch: bool=True, only_space_telescopes: bool=False, exclude_space_telescopes: bool=False, include_positions: bool=False) -> pd.DataFrame:
        """Returns a pandas DataFrame with information about all MPC observatories,
        Carthesian positions are only available after running the `to_tudat()` method.

        Parameters
        ----------
        only_in_batch : bool, optional
            Filter out observatories that are not in the batch, by default True
        only_space_telescopes : bool, optional
            Filter out all observatories except space telescopes, by default False
        only_space_telescopes : bool, optional
            Filter out all space telescopes, by default False
        include_positions : bool, optional
            Include cartesian positions of the terrestrial telescopes only available
            after running to_tudat(), by default False

        Returns
        -------
        pd.DataFrame
            Dataframe with information about the observatories.
        """