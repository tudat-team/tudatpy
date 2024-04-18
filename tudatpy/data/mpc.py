from tudatpy.numerical_simulation import environment_setup  # type:ignore
from tudatpy.numerical_simulation import estimation, environment  # type:ignore
from tudatpy.numerical_simulation.estimation_setup import observation  # type:ignore
from tudatpy.numerical_simulation.environment_setup import add_gravity_field_model
from tudatpy.numerical_simulation.environment_setup.gravity_field import central_sbdb
from tudatpy.data._biases import get_biases_EFCC18, get_weights_VFCC17

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
import datetime

from astroquery.mpc import MPC
from astropy.time import Time
import astropy

from typing import Union, Tuple, List, Dict

import copy


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
    >>> observation_collection = batch.to_tudat(bodies=bodies, included_satellites=None)

    """

    def __init__(self) -> None:
        """Create an empty MPC batch."""
        self._table: pd.DataFrame = pd.DataFrame()
        self._observatories: List[str] = []
        self._space_telescopes: List[str] = []
        self._bands: List[str] = []
        self._MPC_codes: List[str] = []
        self._size: int = 0

        self._observatory_info: Union[pd.DataFrame, None] = None
        self._MPC_space_telescopes: List[str] = []

        self._get_station_info()

        self._epoch_start: float = 0.0
        self._epoch_end: float = 0.0

        self._bodies_created = {}

        self._EFCC18_applied = False
        self._custom_weights_set = False

        # for manual additions of table (from_pandas, from_astropy)
        self._req_cols = ["number", "epoch", "RA", "DEC", "band", "observatory"]

    def __copy__(self):
        new = BatchMPC()

        new._table = copy.deepcopy(self._table)
        new._refresh_metadata()

        return new

    def copy(self) -> "BatchMPC":
        """Create a copy of the batch, equivalent to copy.copy(batchMPC())

        Returns
        -------
        BatchMPC
            Copy of batch.
        """
        return copy.copy(self)

    # getters to make everything read-only
    @property
    def table(self) -> pd.DataFrame:
        """Pandas dataframe with observation data"""
        return self._table

    @property
    def observatories(self) -> List[str]:
        """List of observatories in batch"""
        return self._observatories

    @property
    def space_telescopes(self) -> List[str]:
        """List of satellite_observatories in batch"""
        return self._space_telescopes

    @property
    def MPC_objects(self) -> List[str]:
        """List of MPC objects"""
        return self._MPC_codes

    @property
    def size(self) -> int:
        """Number of observations in batch"""
        return self._size

    @property
    def bands(self) -> List[str]:
        """List of bands in batch"""
        return self._bands

    @property
    def epoch_start(self) -> float:
        """Epoch of oldest observation in batch in seconds since J2000 TDB"""
        return self._epoch_start

    @property
    def epoch_end(self) -> float:
        """Epoch of latest observation in batch in seconds since J2000 TDB"""
        return self._epoch_end

    @property
    def bodies_created(self) -> dict:
        """Dictionary with the bodies created by to_tudat and details."""
        return self._bodies_created

    def __len__(self):
        return self._size

    def __add__(self, other):
        temp = BatchMPC()

        temp._table = (
            pd.concat([self._table, other._table])
            .sort_values("epoch")
            .drop_duplicates()
        )

        temp._refresh_metadata()

        return temp

    # helper functions
    def _refresh_metadata(self) -> None:
        """Internal. Update batch metadata"""
        self._table.drop_duplicates()

        self._observatories = list(self._table.observatory.unique())
        self._space_telescopes = [
            x for x in self._observatories if x in self._MPC_space_telescopes
        ]
        if "bands" in self._table.columns:
            self._bands = list(self._table.band.unique())
        self._MPC_codes = list(self._table.number.unique())
        self._size = len(self._table)

        self._epoch_start = self._table.epochJ2000secondsTDB.min()
        self._epoch_end = self._table.epochJ2000secondsTDB.max()

    def _get_station_info(self) -> None:
        """Internal. Retrieve data on MPC listed observatories"""
        try:
            temp = MPC.get_observatory_codes().to_pandas()  # type: ignore
            # This query checks if Longitude is Nan: non-terretrial telescopes
            sats = list(temp.query("Longitude != Longitude").Code.values)
            self._observatory_info = temp
            self._MPC_space_telescopes = sats
        except Exception as e:
            print("An error occured while retrieving observatory data")
            print(e)

    def _add_observatory_positions(
        self, bodies: environment.SystemOfBodies, earth_name
    ) -> None:
        """Internal. Add observatory cartesian postions to station data"""
        temp = self._observatory_info

        if temp is None:
            # This error will probably never occur
            txt = """Observatory positions can not be assigned.
                  This is likely due to a failure in retrieving the observatories."""
            raise ValueError(txt)

        r_earth = bodies.get(earth_name).shape_model.average_radius

        # Add geocentric cartesian positions
        temp = (
            temp.assign(X=lambda x: x.cos * r_earth * np.cos(np.radians(x.Longitude)))
            .assign(Y=lambda x: x.cos * r_earth * np.sin(np.radians(x.Longitude)))
            .assign(Z=lambda x: x.sin * r_earth)
        )
        self._observatory_info = temp

    def _apply_EFCC18(
        self,
        bias_file=None,
        Nside=None,
        catalog_flags=None,
    ):
        if self._EFCC18_applied:
            pass
        else:
            RA_corr, DEC_corr = get_biases_EFCC18(
                RA=self.table.RA.to_numpy(),
                DEC=self.table.DEC.to_numpy(),
                epochJ2000secondsTDB=self.table.epochJ2000secondsTDB.to_numpy(),
                catalog=self.table.catalog.to_numpy(),
                bias_file=bias_file,
                Nside=Nside,
                catalog_flags=catalog_flags,
            )

            self._table = self._table.assign(RA_EFCC18=lambda x: x.RA - RA_corr)
            self._table = self._table.assign(DEC_EFCC18=lambda x: x.DEC - DEC_corr)
            self._table = self._table.assign(corr_RA_EFCC18=lambda x: RA_corr)
            self._table = self._table.assign(corr_DEC_EFCC18=lambda x: DEC_corr)

            self._EFCC18_applied = True

    # methods for data retrievels
    def get_observations(
        self, MPCcodes: List[Union[str, int]], drop_misc_observations: bool = True
    ) -> None:
        """Retrieve all observations for a set of MPC listed objeccts.
        This method uses astroquery to retrieve the observations from the MPC.
        An internet connection is required, observations are cached for faster subsequent retrieval.
        Removes duplicate and irrelevant observation data by default (see `drop_misc_observations`).

        Parameters
        ----------
        MPCcodes : List[int]
            List of integer MPC object codes for minor planets or and comets.
        drop_misc_observations : List[int]
            Drops observations made by method: radar and offset (natural satellites).
            Drops observations made by roaming observers.
            Drops duplicate listings to denote first observation.
        """

        if not isinstance(MPCcodes, list):
            raise ValueError("MPCcodes parameter must be list of integers/strings")
        for code in MPCcodes:
            if not (isinstance(code, int) or isinstance(code, str)):
                raise ValueError(
                    "All codes in the MPCcodes parameter must be integers or string"
                )

            try:
                obs = MPC.get_observations(code).to_pandas()  # type: ignore

                # convert JD to J2000 and UTC, convert deg to rad
                obs = (
                    obs.assign(
                        epochJ2000secondsTDB=lambda x: (
                            Time(x.epoch, format="jd", scale="utc").tdb.value
                            - 2451545.0
                        )
                        * 86400
                    )
                    .assign(RA=lambda x: np.radians(x.RA))
                    .assign(DEC=lambda x: np.radians(x.DEC))
                    .assign(epochUTC=lambda x: Time(x.epoch, format="jd").to_datetime())
                )

                # drop miscellenous observations:
                if drop_misc_observations:
                    first_discoveries = ["x", "X"]
                    roaming = ["V", "v", "W", "w"]
                    radar = ["R", "r", "Q", "q"]
                    offset = ["O"]
                    observation_types_to_drop = (
                        first_discoveries + roaming + radar + offset
                    )
                    obs = obs.query("note2 != @observation_types_to_drop")

                # convert object mpc code to string
                if "comettype" in obs.columns:
                    # for the case where we have a comet
                    obs.loc[:, "number"] = obs.desig
                else:
                    obs.loc[:, "number"] = obs.number.astype(str)

                self._table = pd.concat([self._table, obs])

            except Exception as e:
                print(f"An error occured while retrieving observations of MPC: {code}")
                print(e)

        self._refresh_metadata()

    def _add_table(self, table: pd.DataFrame, in_degrees: bool = True):
        """Internal. Formats a manually entered table of observations, use from_astropy or from_pandas."""
        obs = table
        obs = obs.assign(
            epochJ2000secondsTDB=lambda x: (
                Time(x.epoch, format="jd", scale="utc").tdb.value - 2451545.0
            )
            * 86400
        ).assign(epochUTC=lambda x: Time(x.epoch, format="jd").to_datetime())
        if in_degrees:
            obs = obs.assign(RA=lambda x: np.radians(x.RA)).assign(
                DEC=lambda x: np.radians(x.DEC)
            )

        # convert object mpc code to string
        obs["number"] = obs.number.astype(str)
        self._table = pd.concat([self._table, obs])
        self._refresh_metadata()

    def from_astropy(
        self, table: astropy.table.QTable, in_degrees: bool = True, frame: str = "J2000"
    ):
        """Manually input an astropy table with observations into the batch. 
        Usefull when manual filtering from astroquery is required
        Table must contain the following columns:
        number - MPC code of the minor planet\\
        epoch - in julian days\\
        RA - right ascension in either degrees or radians (degrees is default)\\
        DEC - declination in either degrees or radians (degrees is default)\\
        band - band of the observation (currently unused)\\
        observatory - MPC code of the observatory

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
        if not (
            isinstance(table, astropy.table.QTable)
            or isinstance(table, astropy.table.Table)
        ):
            raise ValueError(
                "Table must be of type astropy.table.QTable or astropy.table.Table"
            )
        if frame != "J2000":
            txt = "Only observations in J2000 are supported currently"
            raise NotImplementedError(txt)

        # check if all mandatory names are present
        if not set(self._req_cols).issubset(set(table.colnames)):
            txt = f"Table must include a set of mandatory columns: {self._req_cols}"
            raise ValueError(txt)

        self._add_table(table=table.to_pandas(), in_degrees=in_degrees)

    def from_pandas(
        self, table: pd.DataFrame, in_degrees: bool = True, frame: str = "J2000"
    ):
        """Manually input an pandas dataframe with observations into the batch. 
        Usefull when manual filtering from astroquery is required
        Table must contain the following columns:
        number - MPC code of the minor planet\\
        epoch - in julian days\\
        RA - right ascension in either degrees or radians (degrees is default)\\
        DEC - declination in either degrees or radians (degrees is default)\\
        band - band of the observation (currently unused)\\
        observatory - MPC code of the observatory

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
        if not isinstance(table, pd.DataFrame):
            raise ValueError("Table must be of type pandas.DataFrame")
        if frame != "J2000":
            txt = "Only observations in J2000 are supported currently"
            raise NotImplementedError(txt)

        # check if all mandatory names are present
        if not set(self._req_cols).issubset(set(table.columns)):
            txt = f"Table must include a set of mandatory columns: {self._req_cols}"
            raise ValueError(txt)

        self._add_table(table=table, in_degrees=in_degrees)

    def set_weights(
        self,
        weights: Union[list, np.ndarray, pd.Series],
    ):
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
        if len(weights) != len(self.table):
            raise ValueError(
                "Weights vector must have same length as number of observations (BatchMPC.table)."
            )
        if isinstance(weights, pd.Series):
            # this is to avoid errors with indexing in pandas.
            weights = weights.to_numpy()

        self._table.loc["weight"] = weights
        self._custom_weights_set = True

    def filter(
        self,
        bands: Union[List[str], None] = None,
        catalogs: Union[List[str], None] = None,
        observation_types: Union[List[str], None] = None,
        observatories: Union[List[str], None] = None,
        observatories_exclude: Union[List[str], None] = None,
        epoch_start: Union[float, datetime.datetime, None] = None,
        epoch_end: Union[float, datetime.datetime, None] = None,
        in_place: bool = True,
    ) -> Union[None, "BatchMPC"]:
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
            Start date to include observations from, can be in python datetime in utc\
                 or the more conventional tudat seconds since j2000 in TDB if float,\
                     by default None
        epoch_end : Union[float, datetime.datetime, None], optional
            Final date to include observations from, can be in python datetime in utc\
                 or the more conventional tudat seconds since j2000 in TDB if float,\
                     by default None
        in_place : bool, optional
            If true, modify the current batch object.\
                  If false returns a new object that is\
                      filtered, currect batch remains unmodified, by default True

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

        # basic user input handling
        assert not isinstance(
            observatories, int
        ), "stations parameter must be of type 'str' or 'List[str]'"

        if not (isinstance(observatories, list) or (observatories is None)):
            raise ValueError("observatories parameter must be list of strings or None")

        if not (isinstance(catalogs, list) or (catalogs is None)):
            raise ValueError("catalogs parameter must be list of strings or None")

        if not (isinstance(observation_types, list) or (observation_types is None)):
            raise ValueError(
                "observation_types parameter must be list of strings or None"
            )

        if not (
            isinstance(observatories_exclude, list) or (observatories_exclude is None)
        ):
            raise ValueError(
                "observatories_exclude parameter must be list of strings or None"
            )

        if not (isinstance(bands, list) or (bands is None)):
            raise ValueError("bands parameter must be list of strings or None")

        if (observatories is not None) and (observatories_exclude is not None):
            txt = "Include or exclude observatories, not both at the same time."
            raise ValueError(txt)

        if in_place:
            if bands is not None:
                self._table = self._table.query("band == @bands")
            if catalogs is not None:
                self._table = self._table.query("catalog == @catalogs")
            if observation_types is not None:
                self._table = self._table.query("note2 == @observation_types")
            if observatories is not None:
                self._table = self._table.query("observatory == @observatories")
            if observatories_exclude is not None:
                self._table = self._table.query("observatory != @observatories_exclude")
            if epoch_start is not None:
                if isinstance(epoch_start, float) or isinstance(epoch_start, int):
                    self._table = self._table.query(
                        "epochJ2000secondsTDB >= @epoch_start"
                    )
                elif isinstance(epoch_start, datetime.datetime):
                    self._table = self._table.query("epochUTC >= @epoch_start")
            if epoch_end is not None:
                if isinstance(epoch_end, float) or isinstance(epoch_end, int):
                    self._table = self._table.query(
                        "epochJ2000secondsTDB <= @epoch_end"
                    )
                elif isinstance(epoch_end, datetime.datetime):
                    self._table = self._table.query("epochUTC <= @epoch_end")

            self._refresh_metadata()
            return None
        else:
            new = self.copy()
            new.filter(
                bands=bands,
                catalogs=catalogs,
                observation_types=observation_types,
                observatories=observatories,
                observatories_exclude=observatories_exclude,
                epoch_start=epoch_start,
                epoch_end=epoch_end,
                in_place=True,
            )
            return new

    def to_tudat(
        self,
        bodies: environment.SystemOfBodies,
        included_satellites: Union[Dict[str, str], None],
        station_body: str = "Earth",
        add_sbdb_gravity_model: bool = False,
        apply_weights_VFCC17: bool = True,
        apply_star_catalog_debias: bool = True,
        debias_kwargs: dict = dict(),
    ) -> estimation.ObservationCollection:
        """Converts the observations in the batch into a Tudat compatible format and
          sets up the relevant Tudat infrastructure to support estimation.
        This method does the following:\\
            1. (By Default) Applies star catalog debiasing and estimation
            weight calculation.\\
            2. Creates an empty body for each minor planet with their MPC code as a 
            name.\\
            3. Adds this body to the system of bodies inputted to the method.\\
            4. Retrieves the global position of the terrestrial observatories in 
            the batch and adds these stations to the Tudat environment.\\
            5. Creates link definitions between each unique terrestrial 
            observatory/ minor planet combination in the batch.\\
            6. (Optionally) creates a link definition between each 
            space telescope / minor planet combination in the batch. 
            This requires an addional input.\\
            7. Creates a `SingleObservationSet` object for each unique link that 
            includes all observations for that link.\\
            8. (By Default) Add the relevant weights to the `SingleObservationSet`
            per observation.\\
            8. Returns the observations


        Parameters
        ----------
        bodies : environment.SystemOfBodies
            SystemOfBodies containing at least the earth to allow for the placement of 
            terrestrial telescopes
        included_satellites : Union[Dict[str, str], None], optional
            A dictionary that links the name of a space telescope used by the user with 
            the observatory code in MPC. Used when utilising observations from space 
            telescopes like TESS and WISE. The keys should be the MPC observatory 
            codes. The values should be the bodies' in the user's code. The relevant 
            observatory code can be retrieved using 
            the .observatories_table() method, by default None
        station_body : str, optional
            Body to attach ground stations to. Does not need to be changed unless the 
            `Earth` body has been renamed, by default "Earth"
        station_body : bool, optional
            Adds a central_sbdb gravity model to the object, generated using JPL's small body database. 
            This option is only available for a limited number of bodies and raises an error if unavailable. 
            See tudatpy.numerical_simulation.environment_setup.gravity_field.central_sbdb for more info.
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
        estimation.ObservationCollection
            ObservationCollection with the observations found in the batch

        Examples
        ----------
        Using Space Telescopes

        Create dictionary to link name in tudat with mpc code in batch:

        >> sats_dict = {"C57":"TESS"}
        >> obs = batch1.to_tudat(bodies, included_satellites=sats_dict)


        """
        # apply EFCC18 star catalog corrections:
        if apply_star_catalog_debias:
            self._apply_EFCC18(**debias_kwargs)
            RA_col = "RA_EFCC18"
            DEC_col = "DEC_EFCC18"
        else:
            RA_col = "RA"
            DEC_col = "DEC"

        # Calculate observation weights and update table:
        if apply_weights_VFCC17:
            temp_table: pd.DataFrame = get_weights_VFCC17(  # type:ignore
                mpc_table=self.table,
                return_full_table=True,
            )

            self._table.loc[:, "weight"] = temp_table.loc[:, "weight"].values
        else:
            temp_table = self.table.copy()

        # start user input validation
        # Ensure that Earth is in the SystemOfBodies object
        try:
            bodies.get(station_body)
        except Exception as e:
            print(
                f"Body {station_body} is not in bodies, if you have renamed Earth, "
                + "set the station_body paramater to the new name."
            )
            raise e

        # Add positions of the observatories
        self._add_observatory_positions(bodies, station_body)

        # get satellites to include and exclude
        if included_satellites is not None:
            sat_obs_codes_included = list(included_satellites.keys())
            # this appears unused but is used in a pandas query:
            sat_obs_codes_excluded = list(
                set(self._space_telescopes) - set(sat_obs_codes_included)
            )

            # Ensure that the satellite is in the SystemOfBodies object
            for sat in list(included_satellites.values()):
                try:
                    bodies.get(sat)
                except Exception as e:
                    print(f"Body {sat} is not in bodies")
                    raise e
        else:
            sat_obs_codes_included = []
            # this appears unused but is used in a pandas query:
            sat_obs_codes_excluded = self._space_telescopes

        # end user input validation

        # get relevant stations positions
        tempStations = self._observatory_info.query("Code == @self.observatories").loc[
            :, ["Code", "X", "Y", "Z"]
        ]

        # add station positions to the observations
        observations_table = pd.merge(
            left=temp_table,  # type: ignore
            right=tempStations,
            left_on="observatory",
            right_on="Code",
            how="left",
        )
        # remove the observations from unused satellite observatories
        observations_table = observations_table.query("Code != @sat_obs_codes_excluded")

        # add asteroid bodies to SystemOfBodies object
        for body in self._MPC_codes:
            if not bodies.does_body_exist(body):
                bodies.create_empty_body(str(body))
                # save the name of the body added
                self._bodies_created[str(body)] = "empty body"

                if add_sbdb_gravity_model:
                    add_gravity_field_model(bodies, str(body), central_sbdb(str(body)))
                    self._bodies_created[str(body)] = "empty body + sbdb gravity"

        # add ground stations to the earth body
        for idx in range(len(tempStations)):
            station_name = tempStations.iloc[idx].Code

            # skip if it is a satellite observatory
            if station_name in self._space_telescopes:
                continue

            ground_station_settings = environment_setup.ground_station.basic_station(
                station_name=station_name,
                station_nominal_position=[
                    tempStations.iloc[idx].X,
                    tempStations.iloc[idx].Y,
                    tempStations.iloc[idx].Z,
                ],
            )

            # Check if station already exists
            if station_name not in bodies.get(station_body).ground_station_list:
                # Add the ground station to the environment
                environment_setup.add_ground_station(
                    bodies.get_body(station_body), ground_station_settings
                )

        # get unique combinations of mpc bodies and observatories
        unique_link_combos = (
            observations_table.loc[:, ["number", "observatory"]].drop_duplicates()
        ).values

        observation_set_list = []
        for combo in unique_link_combos:
            MPC_number = combo[0]
            station_name = combo[1]

            # CREATE LINKS
            link_ends = dict()

            # observed body link
            link_ends[observation.transmitter] = observation.body_origin_link_end_id(
                MPC_number
            )

            if station_name in sat_obs_codes_included:
                # link for a satellite
                sat_name = included_satellites[station_name]
                link_ends[observation.receiver] = observation.body_origin_link_end_id(
                    sat_name
                )
                link_definition = observation.link_definition(link_ends)
            else:
                # link for a ground station
                link_ends[observation.receiver] = (
                    observation.body_reference_point_link_end_id(
                        station_body, station_name
                    )
                )
                link_definition = observation.link_definition(link_ends)

            # get observations, angles and times for this specific link
            observations_for_this_link = observations_table.query(
                "number == @MPC_number"
            ).query("observatory == @station_name")

            observation_angles = observations_for_this_link.loc[
                :, [RA_col, DEC_col]
            ].to_numpy()

            observation_times = observations_for_this_link.loc[
                :, ["epochJ2000secondsTDB"]
            ].to_numpy()[:, 0]

            # create a set of obs for this link
            observation_set = estimation.single_observation_set(
                observation.angular_position_type,
                link_definition,
                observation_angles,
                observation_times,
                observation.receiver,
            )

            # apply weights if apply_weights is True or set_weights() has been used.
            if apply_weights_VFCC17 or self._custom_weights_set:
                observation_weights = observations_for_this_link.loc[
                    :, ["weight"]
                ].to_numpy()[:, 0]
                observation_weights = np.concatenate(
                    [observation_weights, observation_weights]
                ).flatten()
                observation_set.weights_vector = observation_weights

            observation_set_list.append(observation_set)

        observation_collection = estimation.ObservationCollection(observation_set_list)
        return observation_collection

    def plot_observations_temporal(
        self,
        objects: Union[List[str], None] = None,
        figsize: Tuple[float] = (9.0, 6.0),
    ):
        """Generates a matplotlib figure with the declination and right ascension
        over time.


        Parameters
        ----------
        objects : Union[List[str], None], optional
            List of specific MPC objects in batch to plot, None to plot all
            , by default None
        projection : str, optional
            projection of the figure options are: 'aitoff', 'hammer',
            'lambert' and 'mollweide', by default "aitoff"
        figsize : Tuple[float], optional
            size of the matplotlib figure, by default (15, 7)

        Returns
        -------
        Matplotlib figure
            Matplotlib figure with the observations
        """
        fig, (axRA, axDEC) = plt.subplots(2, 1, figsize=figsize, sharex=True)

        if objects is None:
            objs = self.MPC_objects
        else:
            objs = [str(o) for o in objects]

        for obj in objs:
            tab = self._table.query("number == @obj")

            axRA.scatter(
                tab.epochUTC,
                np.degrees(tab.RA),
                label="MPC: " + obj,
                marker=".",
                # linestyle=None,
            )
            axDEC.scatter(
                tab.epochUTC,
                np.degrees(tab.DEC),
                label="MPC: " + obj,
                marker=".",
                # linestyle=None,
            )

        axRA.set_ylabel(r"Right Ascension $[\deg]$")
        axDEC.set_ylabel(r"Declination $[\deg]$")
        axDEC.set_xlabel(r"Time [year]")

        axRA.grid()
        axDEC.grid()
        axRA.set_title("Right Ascension")
        axDEC.set_title("Declination")
        axRA.legend(bbox_to_anchor=(1.01, 1), loc="upper left")

        buffer = 10
        axRA.set_ylim(0 - buffer, 360 + buffer)
        axDEC.set_ylim(-90 - buffer, 90 + buffer)

        fig.set_layout_engine("tight")

        return fig

    def plot_observations_sky(
        self,
        objects: Union[List[str], None] = None,
        projection: Union[str, None] = None,
        figsize: Tuple[float] = (14.0, 7.0),
    ):
        """Generates a matplotlib figure with the observations'
        right ascension and declination over time.

        Parameters
        ----------
        objects : Union[List[str], None], optional
            List of specific MPC objects in batch to plot, None to plot all
            , by default None
        projection : str, optional
            projection of the figure options are: 'aitoff', 'hammer',
            'lambert' and 'mollweide', by default "aitoff"
        figsize : Tuple[float], optional
            size of the matplotlib figure, by default (15, 7)

        Returns
        -------
        Matplotlib figure
            Matplotlib figure with the observations
        """
        fig, ax = plt.subplots(
            1, 1, subplot_kw={"projection": projection}, figsize=figsize
        )

        if objects is None:
            objs = self.MPC_objects
        else:
            objs = [str(o) for o in objects]

        markers = ["o", "+", "^"]

        vmin = mdates.date2num(self._table.query("number == @objs").epochUTC.min())
        vmax = mdates.date2num(self._table.query("number == @objs").epochUTC.max())

        for idx, obj in enumerate(objs):
            tab = self._table.query("number == @obj")

            a = plt.scatter(
                (tab.RA) - np.pi if projection is not None else np.degrees(tab.RA),
                (tab.DEC) if projection is not None else np.degrees(tab.DEC),
                s=30,
                marker=markers[int(idx % len(markers))],
                cmap=cm.plasma,
                label="MPC: " + obj,
                c=mdates.date2num(tab.epochUTC),
                vmin=vmin,
                vmax=vmax,
            )

        ax.legend(markerscale=2)
        ax.set_xlabel(r"Right Ascension $[\deg]$")
        ax.set_ylabel(r"Declination $[\deg]$")

        yticks = [
            f"{x}°"
            for x in (np.degrees(np.array(ax.get_yticks().tolist()))).astype(int)
        ]
        xticks = [
            f"{x}°"
            for x in (np.degrees(np.array(ax.get_xticks().tolist())) + 180).astype(int)
        ]
        if projection in ["aitoff", "hammer", "mollweide"]:
            ax.set_yticklabels(yticks)
            ax.set_xticklabels(xticks)
        elif projection == "lambert":
            ax.set_xticklabels(xticks)
        else:
            pass

        ticks = np.linspace(vmin, vmax, 7)
        labels = [mdates.num2date(t).strftime("%d %b %Y") for t in ticks]
        cbar = plt.colorbar(mappable=a, ax=ax, label="Time", ticks=ticks)

        cbar.ax.set_yticklabels(labels)

        startUTC = self._table.query("number == @objs").epochUTC.min()
        endUTC = self._table.query("number == @objs").epochUTC.max()
        ax.grid()
        fig.suptitle(f"{self.size} observations between {startUTC} and {endUTC}")

        fig.set_layout_engine("tight")

        return fig

    def summary(self):
        """Produce a quick summary of the content of the batch"""
        print()
        print("   Batch Summary:")
        print(f"1. Batch includes {len(self._MPC_codes)} minor planets:")
        print("  ", self.MPC_objects)
        satObs = len(self._table.query("observatory == @self.space_telescopes"))
        print(
            f"2. Batch includes {self.size} observations, including {satObs} "
            + "observations from space telescopes"
        )
        print(
            f"3. The observations range from {self._table.epochUTC.min()} "
            + f"to {self._table.epochUTC.max()}"
        )
        print(f"   In seconds TDB since J2000: {self.epoch_start} to {self.epoch_end}")
        print(
            f"   In Julian Days: "
            + f"{self._table.epoch.min()}"
            + f" to {self._table.epoch.max()}"
        )
        print(
            f"4. The batch contains observations from {len(self.observatories)} "
            + f"observatories, including {len(self.space_telescopes)} space telescopes"
        )
        print()

    def observatories_table(
        self,
        only_in_batch: bool = True,
        only_space_telescopes: bool = False,
        exclude_space_telescopes: bool = False,
        include_positions: bool = False,
    ) -> pd.DataFrame:
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
        temp = self._observatory_info
        temp2 = self._table

        count_observations = (
            temp2.groupby("observatory")
            .count()
            .rename(columns={"number": "count"})
            .reset_index(drop=False)
            .loc[:, ["observatory", "count"]]
        )

        temp = pd.merge(
            left=temp,
            right=count_observations,
            left_on="Code",
            right_on="observatory",
            how="left",
        )
        if only_in_batch:
            temp = temp.query("Code == @self.observatories")
        if only_space_telescopes:
            temp = temp.query("Code == @self._MPC_space_telescopes")
        if exclude_space_telescopes:
            temp = temp.query("Code != @self._MPC_space_telescopes")
        if not include_positions:
            temp = temp.loc[:, ["Code", "Name", "count"]]
        return temp
