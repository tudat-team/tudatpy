import astropy
from astroquery.jplhorizons import Horizons
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.table import vstack

import math
import numpy as np
import pandas as pd

from typing import Union, List, Tuple
import datetime

from tudatpy import constants
from tudatpy.numerical_simulation import environment_setup

import re

# NOTE when updating, also update tudatpy.numerical_simulation.environment_setup.ephemeris.horizons

class HorizonsQuery:
    """This class provides an interface to JPL's Horizon System. 
    JPL Horizons provides access to highly accurate ephemerides for many solar system objects,
    including asteriuds, comets, planets, moons and select spacecraft.
    The class extends astroquery's to cater to the needs of Tudat users,
    while maintaining compatibility with all of astroquery's features.

    There are some notable differences:

    Time input has been simplified to reduce ambiguities:
    - List of times are given in seconds since J2000 TDB.
    - Start can be given in datetime format or seconds since J2000 TDB.
    - Timesteps like months and years are not permitted.

    And some additional features:
    - Extended query allows data retrieval limits
    to be broken by automatically splitting up a query into multiple subqueries
    and combining the data.
    - Ephemeris settings can automatically be generated using Vectors API.


    """

    query_limit = 90024

    # 50 is chosen to be conservative,
    # 50 is chosen to prevent scenarios where higher precision is used,
    # to prevent undefined behaviour when very futuristic dates are used,
    # and because it is a nice human friendly number.
    # Testing found a limit of 77 but this will not be conservative.
    query_limit_list = 50

    def __init__(
        self,
        query_id: str,
        location: str,
        query_type: str = "default",
        epoch_start: Union[datetime.datetime, float, None] = None,
        epoch_end: Union[datetime.datetime, float, None] = None,
        epoch_step: Union[str, None] = None,
        epoch_list: Union[list, None] = None,
        extended_query: bool = False,
    ) -> None:
        """Query object to retrieve Horizons data for a singular body.

        Parameters
        ----------
        query_id : str
            Query term for object to retrieve data for. The query may be ambiguous, 
            in which case suggestions for exact objects are presented. 
            The behaviour for this parameter changes based on the `query_type` 
            parameter. 
            The default behaviour is the same as the official Horizons API.

            Here are some examples for the default behaviour: 

            - `Earth` - ambiguous, will 
            suggest 3 (Earth-Moon Barycentre) and 399 (Earth).

            - `3` - will retrieve Earth-Moon Barycentre.
            
            - `3;` - semi-colon searches for minor planets.\
            In this case it will search for the minor planet with MPC code 3: Juno.

            - `-3` - A minus sign searches for spacecraft. 
            In this case the Mars Orbiter Mission.

            See the Horizons System manual 
            for an extensive explanation on this parameter and the timespans 
            page for some examples:
            https://ssd.jpl.nasa.gov/horizons/manual.html#select 

            https://ssd.jpl.nasa.gov/horizons/time_spans.html 

        location : str
            Coordinate centre for the data with syntax `site@body`. 

            In general, the syntax is: `site@body` but there are several shorthands.
            The sight may be a specific location on the body, such as an observatory.
            500 specifies the geocentre, and can be left out as a shorthand 
            (500@399 = @399).
            Here are some examples:

            - `@10` or `@Sun` - if no site is given, the geocenter will be taken.

            - `500@399` or `@399` or `500` - Geocentric. The site defaults to `500`.
            The body defaults to `@399`.

            - `@0` or `@SSB` - The solar system barycentre.

            - `0` - without the `@` symbol, this location is equivalent to `0@399` 
            which is the observatory with MPC code 0 on Earth and not the SSB.
            In this case: Greenwich Observatory.

            - `Greenwich` - Equivalent to `Greenwich@399`, Greenwich Observatory.

            - `Earth` - Heaven on Earth Observatory, not equivalent to `@399`.

            - `@Earth` - ambiguous, will suggest `@399` or `@3`.

            - `Geocentric` - Equivalent to `500@399`.

            As mistakes can easily made, it is highly recommended to consult 
            the manual and use specific codes for this parameter: 

            https://ssd.jpl.nasa.gov/horizons/manual.html#center
        query_type : str, optional
            The query type constrains the search for the `query_id` parameter. 
            Can have values: 'default', 'None', 'smallbody', 'designation', 'name', 
            'asteroid_name', 'comet_name', by default "default"

            While all objects can be found with the default parameter, 
            some queries have ambiguities. Specifying the query type eliminates the 
            specific syntax. For example searching Juno with `query_type` `spacecraft` 
            retrieves data for the spacecraft while setting `smallbody` 
            retrieves data for the asteroid Juno. 
            See the Horizons system Manual and astroquery documentation for more info: 

            https://ssd.jpl.nasa.gov/horizons/manual.html#select

            https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html

        epoch_start : Union[datetime.datetime, float, None], optional
            Starting date to retrieve data for. 
            Must be either a python datetime object or a float seconds since J2000 TDB.
            Combined with `epoch_end` and 
            `epoch_step` can be used to retrieve a range of data. 
            If `epoch_list` is used, value must be None, by default None
        epoch_end : Union[datetime.datetime, float, None], optional
            Final date to retrieve data for. 
            Must be either a python datetime object or a float seconds since J2000 TDB.
            If the start stop, and step parameters dont result in an 
            integer multiple of values,the final date is not used. 
            Combined with `epoch_end` and 
            `epoch_step` can be used to retrieve a range of data. 
            If `epoch_list` is used, value must be None, by default None
        epoch_step : Union[str, None], optional
            Step for the range of epochs to retrieve data for. 
            Can be be either a specific timestep such as `15m`,
            or a number of partitions `10` for 10 equidistant steps within the range.
            For the timestep, quantifiers `m` for minute, `h` for hour and `d` for day 
            can be used. Seconds are not available, consider using the `epoch_list` 
            parameter or a number of partitions instead. 
            Month and Year are normally available in JPL Horizons but are restricted 
            here because they may produce ambiguous results (leap years etc.).
            If `epoch_list` is used, value must be None, by default None
        epoch_list : Union[list, None], optional
            List of times in seconds since J2000 TDB. Can be used 
            to retrieve specific times instead of a range. 
            Must be None if start, end and step are set, by default None
        extended_query : bool, optional
            Enables the retrieval of larger collections of data, by default False. 
            Horizons System has a limit on how much data can be output (90024) lines. 
            Additionally there is a limit on input length which may be exceeded if a 
            long list of times is given with the `epoch_list` option (50 epochs). 
            When enabled, the times are split in to multiple subqueries. 
            The data of these subqueries is then combined and presented as a single 
            query. Using a very large number of queries may make data retrieval slow, 
            especially when using the list option. Consider whether your use case 
            requires a large number of specific times, 
            or if a range can be used instead.

            The epoch step partitions (`10` not `10h`) is currently unsupported with 
            the extended_query.

        """

        self._query_id = query_id
        self._query_type = None if query_type == "default" else query_type

        # this gets set after vectors is run
        self._target_full_name = None

        # these get set if the above is not None and the name property gets called
        self._MPC_number = None
        self._name = None
        self._designation = None
        self.epoch_type = None

        # the 'subqueries' when using an extended query and their lengths.
        self.queries = []
        self.query_lengths = []

        # ##############
        # error handling
        # ##############

        # epoch range format:
        # epoch_list IS none rest is NOT none
        if epoch_list is None:
            if (epoch_start is None) or (epoch_end is None) or (epoch_step is None):
                raise ValueError(
                    "Must specify either a list of times in sec since J2000 "
                    + "or a combined start, end and step parameters"
                )
            else:
                ts_seconds, num_lines = self._interpret_timestep(
                    epoch_step, epoch_start, epoch_end
                )

                self._validate_time_range(epoch_start, epoch_end)

                if re.findall(r"\d+", epoch_step)[0] == epoch_step:
                    self.epoch_type = "partition"
                else:
                    self.epoch_type = "range"

        # epoch list format:
        # start step end IS none, list NOT none
        elif (epoch_start is None) and (epoch_end is None) and (epoch_step is None):
            if epoch_list is None:
                raise ValueError(
                    "Must specify either a list of times in sec since J2000 "
                    + "or a combined start, end and step parameters not both"
                )
            elif (
                isinstance(epoch_list, list)
                or isinstance(epoch_list, tuple)
                or isinstance(epoch_list, np.ndarray)
            ):
                self.epoch_type = "list"
                num_lines = len(epoch_list)

            else:
                raise ValueError(
                    "epoch_list must be a list of times in sec since J2000"
                )
        # Neither format indicated
        else:
            raise ValueError(
                "Must specify either a list of times in sec since J2000 "
                + "or a combined start, end and step parameters"
            )

        # ########################################
        # Check for query limit and create queries
        # ########################################
        if (
            (not extended_query)
            and (num_lines > HorizonsQuery.query_limit)
            and (self.epoch_type == "range")
        ):
            txt = (
                "The query is larger than Horizon's limit of "
                + f"{HorizonsQuery.query_limit} lines: ({num_lines}). "
                + f"Consider making the timestep higher than {epoch_step} "
                + "or use the 'extended_query' option."
            )

            raise ValueError(txt)
        elif (
            (not extended_query)
            and (num_lines > HorizonsQuery.query_limit)
            and (self.epoch_type == "partition")
        ):
            txt = (
                "The query is larger than Horizon's limit of "
                + f"{HorizonsQuery.query_limit} lines: ({num_lines}). "
                + f"Consider using fewer partitions than {epoch_step} "
                + "or use the 'extended_query' option."
            )

            raise ValueError(txt)
        elif (
            (not extended_query)
            and (num_lines > HorizonsQuery.query_limit_list)
            and (self.epoch_type == "list")
        ):
            txt = (
                "The query is larger than Horizon's input limit of "
                + f"{HorizonsQuery.query_limit_list} lines: ({num_lines}). "
                + "Consider selecting fewer epochs to poll, "
                + f"using the range option for up to {HorizonsQuery.query_limit} "
                + "epochs (or higher with the 'extended_query' option). "
                + "Using the 'extended_query' option with lists may be slow for "
                + "a high number of epochs (>200)."
            )

            raise ValueError(txt)

        # query is smaller than limit -> one batch
        # seperate check for list as the num lines is smaller
        elif (
            (self.epoch_type != "list") and (num_lines < HorizonsQuery.query_limit)
        ) or (
            (self.epoch_type == "list") and (num_lines < HorizonsQuery.query_limit_list)
        ):
            if self.epoch_type == "list":
                # convert seconds since J2000 TDB to JD TDB
                epoch_def = self._format_time_list(epoch_list)
            else:
                epoch_def = dict(
                    start=self._format_time_range(epoch_start),
                    stop=self._format_time_range(epoch_end),
                    step=epoch_step,
                )

            self.queries.append(
                Horizons(
                    id=self._query_id,
                    location=location,
                    id_type=self._query_type,
                    epochs=epoch_def,
                )
            )
            self.query_lengths.append(num_lines)

        # query is extended type -> multiple batches
        elif extended_query:
            # case where its a list -> split list
            if self.epoch_type == "list":
                num_splits = math.ceil(num_lines / HorizonsQuery.query_limit_list)

                epoch_def = self._format_time_list(epoch_list)
                splits = np.array_split(epoch_def, num_splits)

                # makes a set of queries that are at the query limit, plus 1 smaller one
                for _, split in enumerate(splits):
                    self.queries.append(
                        Horizons(
                            id=self._query_id,
                            location=location,
                            id_type=self._query_type,
                            epochs=list(split),
                        )
                    )
                    self.query_lengths.append(len(split))

            # Case where it is a range.
            elif (self.epoch_type == "range") or (self.epoch_type == "partition"):
                if self.epoch_type == "partition":
                    raise NotImplementedError(
                        "Using number of divisions for time "
                        + "range with extended queries is "
                        + "currently unsupported, please use "
                        + "timesteps like 15m instead."
                    )

                start_astro = self._convert_time_to_astropy(epoch_start)
                end_astro = self._convert_time_to_astropy(epoch_end)

                ts_seconds = ts_seconds * u.second

                max_query_step = TimeDelta(
                    ((HorizonsQuery.query_limit - 1) * ts_seconds), format="sec"
                )

                next_start = start_astro
                next_limit = start_astro + max_query_step

                formatt = r"%Y-%m-%d %H:%M:%S.%f"

                while next_limit < end_astro:
                    query_len = math.ceil((next_limit - next_start) / ts_seconds)
                    self.query_lengths.append(query_len)

                    if self.epoch_type == "partition":
                        epoch_step_batch = str(query_len)
                    else:
                        epoch_step_batch = epoch_step

                    epoch_def = dict(
                        start=next_start.strftime(formatt),
                        stop=next_limit.strftime(formatt),
                        step=epoch_step_batch,
                    )

                    self.queries.append(
                        Horizons(
                            id=self._query_id,
                            location=location,
                            id_type=self._query_type,
                            epochs=epoch_def,
                        )
                    )

                    next_start = next_limit + ts_seconds
                    next_limit = next_start + max_query_step

                # the final batch:
                query_len = math.ceil((end_astro - next_start) / ts_seconds)
                self.query_lengths.append(query_len)

                if self.epoch_type == "partition":
                    epoch_step_batch = str(query_len)
                else:
                    epoch_step_batch = epoch_step

                epoch_def = dict(
                    start=next_start.strftime(formatt),
                    stop=end_astro.strftime(formatt),
                    step=epoch_step_batch,
                )

                self.queries.append(
                    Horizons(
                        id=self._query_id,
                        location=location,
                        id_type=self._query_type,
                        epochs=epoch_def,
                    )
                )

            else:
                raise RuntimeError("Undefined behaviour, (extended)")

        else:
            raise RuntimeError("Undefined behaviour, (unknown error in init)")

        self.num_queries = len(self.queries)

        self._object_type = None

    @property
    def name(self) -> Union[str, None]:
        """Retrieve the name of the query's object.
        The name is infered from the data retrieved and will return none
        if data has not been retrieved yet.
        Unnamed minor planets will use their designation instead.
        If a name can not be infered, the raw name from Horizons will be returned.
        Please consider raising an issue on the Tudat github in such cases."""
        try:
            self._infer_name()
            return self._name
        except Exception as _:
            print(
                f"Unable to infer name, will use full designation instead: {self._target_full_name}"
            )
            return self._target_full_name

    @property
    def MPC_number(self) -> Union[str, None]:
        """Retrieve the MPC (Minor Planet Centre) number of the object.
        The MPC number is infered from data retrieved and will return none
        if data has not been retrieved yet.
        The MPC number is only relevant to minor planets such as asteroids, TNOs and
        Near-Earth Asteroids."""
        try:
            self._infer_name()
            return self._MPC_number
        except Exception as _:
            return None

    @property
    def designation(self) -> Union[str, None]:
        """Retrieve the relevant designation of the query's object.
        The designation is infered from the data retrieved and will return none
        if data has not been retrieved yet.
        Minor planets and Comets will return their provisional designation
        (1898 DQ for Eros, 1982 HG1 for Halley).
        A comets' formal designation can often be retrieved using the `name` property
        Spacecraft and Major Planets/ Moons will return their JPL number
        (-28 for JUICE, 6 for Saturn Barycentre)."""
        try:
            self._infer_name()
            return self._designation
        except Exception as _:
            return None

    def _infer_name(self) -> None:
        """Internal Method to extract the name, designation and MPC code for the
        object from the data."""
        if self._target_full_name is None:
            self._name = None
        else:
            num_between_brackets = re.findall(r"\((.*?)\)", self._target_full_name)

            # comet
            if ("/" in self._target_full_name) and (
                "S/2" not in self._target_full_name
            ):
                self._object_type = "comet"
                temp = self._target_full_name.split("/")

                self._name = temp[1]
                self._designation = temp[0]
            # planet/moon
            elif (
                not (self._target_full_name[0].isnumeric())
                and num_between_brackets[0].isnumeric()
            ):
                self._object_type = "major"
                temp = self._target_full_name.split(None, 1)
                self._name = temp[0]
                self._designation = num_between_brackets[0]
            # spacecraft
            elif "spacecraft" in self._target_full_name.lower():
                self._object_type = "spacecraft"
                self._name = re.findall(r"^[^\(]*", self._target_full_name)[0]
                try:
                    self._designation = (
                        re.findall(r"\((.*?)\)", self._target_full_name)[1]
                        .replace("(", "")
                        .replace(")", "")
                    )
                except Exception as e:
                    self._designation = None
            # asteroids
            else:
                self._object_type = "minorplanet"
                temp = self._target_full_name.split(None, 1)
                self._MPC_number = temp[0]
                self._name = re.sub(r"\(.+?\)\s*", "", temp[1]).strip()
                self._designation = (
                    re.findall(r"\((.*?)\)", temp[1])[0]
                    .replace("(", "")
                    .replace(")", "")
                )
                # this is for unnamed asteroids
                if len(self._name) == 0:
                    self._name = self._designation

    def _convert_time_to_astropy(self, time) -> astropy.time.Time:
        """Internal method to convert inputted times to astropy"""
        # time is tudat format: seconds since j2000 TDB
        if (
            isinstance(time, float)
            or isinstance(time, int)
            or isinstance(time, np.ndarray)
        ):
            # convert to julian days
            time = (time / constants.JULIAN_DAY) + constants.JULIAN_DAY_ON_J2000
            time_astro = Time(time, format="jd", scale="tdb", precision=9)
        # time is python datetime
        else:
            time_astro = Time(time, format="datetime", scale="tdb", precision=9)

        return time_astro

    def _validate_time_range(self, start, end) -> None:
        """Internal method to validate start/end user input"""
        if not (
            isinstance(start, float)
            or isinstance(start, int)
            or isinstance(start, datetime.datetime)
        ):
            raise TypeError(
                "Incorrect start time given, must be datetime "
                + "object or float seconds since J2000 TDB"
            )
        if not (
            isinstance(end, float)
            or isinstance(end, int)
            or isinstance(end, datetime.datetime)
        ):
            raise TypeError(
                "Incorrect start time given, must be datetime "
                + "object or float seconds since J2000 TDB"
            )

    def _format_time_range(self, time: Union[float, datetime.datetime]) -> str:
        """Internal method to format input time ranges to a JPL accepted format"""
        # https://ssd.jpl.nasa.gov/horizons/manual.html#time
        if isinstance(time, float) or isinstance(time, int):
            # convert sec since J2000 TDB to JD TDB
            formatt = r"%Y-%m-%d %H:%M:%S.%f"
            return self._convert_time_to_astropy(time).strftime(formatt)
        elif isinstance(time, datetime.datetime):
            # using recommended 3 letter months may have compatibility issues with locales
            formatt = r"%Y-%m-%d %H:%M:%S.%f"
            return str(time.strftime(formatt))
        else:
            raise TypeError(
                "Incorrect time value given, must be "
                + "datetime object or seconds since J2000 TDB"
            )

    def _format_time_list(self, times) -> List[astropy.time.Time]:
        """Internal method to format input time lists to a JPL accepted format"""
        # This is accurate enough considering the precision JPL takes as input
        times = self._convert_time_to_astropy((np.array(times))).jd
        return times

    def _interpret_timestep(
        self,
        timestep: str,
        start: Union[datetime.datetime, float],
        end: Union[datetime.datetime, float],
    ):
        """Internal method to determine the length of the timestep in seconds
        This is to help partition the queries in the case of an extended query"""
        numerical_part = re.findall(r"\d+", timestep)
        if len(numerical_part) == 0:
            raise ValueError("Timestep is incorrect")
        numerical_part = float(numerical_part[0])
        alpha_part = re.sub(r"\d+", "", timestep).replace(" ", "").lower()
        time_seconds = None

        ambiguoustxt = (
            "Time intervals like month and year "
            + "are ambigious, please reformulate the timestep in days instead"
        )
        if alpha_part.startswith("y") or alpha_part.startswith("mo"):
            raise ValueError(ambiguoustxt)
        elif alpha_part.startswith("d"):
            time_seconds = 86400 * numerical_part
        elif alpha_part.startswith("h"):
            time_seconds = 3600 * numerical_part
        elif alpha_part.startswith("m"):
            time_seconds = 60 * numerical_part
        # if the timestep is a partition (no hour/day suffix)
        elif len(alpha_part) == 0:
            time_seconds = None
        else:
            raise ValueError("Unrecognized time step, use '1d', '1min', '2 hours' etc.")

        start_astro = self._convert_time_to_astropy(start)
        end_astro = self._convert_time_to_astropy(end)

        duration = (end_astro - start_astro).sec

        if time_seconds is None:
            time_seconds = duration / numerical_part
            steps = numerical_part
            return time_seconds, steps
        else:
            steps = math.ceil(duration / time_seconds)
            return time_seconds, steps

    # ################
    # END USER METHODS
    # ################

    def vectors(
        self,
        frame_orientation: str = "ECLIPJ2000",
        aberations: str = "geometric",
    ) -> astropy.table.Table:
        """Retrieve Horizons Vectors api data in raw astropy format.
        For general purposes, use the `.cartesian()` method instead.

        Parameters
        ----------
        frame_orientation : str, optional
            Reference Frame Orientation, equivalent to the astroquery refplane
            parameter. Options are 'J2000'/'earth', 'ECLIPJ2000'/'ecliptic' and 'body'
            , by default "ECLIPJ2000".
            See astroquery documentation for information about the 'body' option:

            https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html
        aberations : str, optional
            Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
            'apparent', by default "geometric".

            See the Horizons System Manual for more info:

            https://ssd.jpl.nasa.gov/horizons/manual.html#output

        Returns
        -------
        astropy.table.Table
            Unprocessed vectors API data in astropy Table format.
        """
        # User input handling
        if aberations not in ["geometric", "astrometric", "apparent"]:
            raise ValueError(
                "refplane parameter must be one of: "
                + "'geometric', 'astrometric', 'apparent'"
            )
        if frame_orientation not in [
            "ECLIPJ2000",
            "J2000",
            "ecliptic",
            "earth",
            "body",
        ]:
            raise ValueError(
                "refplane parameter must be one of: "
                + '"ECLIPJ2000", "J2000", "ecliptic", "earth", "body"'
            )

        # convert refplanet to one compatible with Horizons
        if frame_orientation == "J2000":
            frame_orientation = "earth"
        elif frame_orientation == "ECLIPJ2000":
            frame_orientation = "ecliptic"
        else:
            frame_orientation = frame_orientation

        res_list = []
        for query in self.queries:
            res = query.vectors(
                refplane=frame_orientation,
                aberrations=aberations,
            )

            res_list.append(res)

        raw = vstack(res_list)

        # retrieve the object's full name from Horizons
        self._target_full_name = raw["targetname"][0]

        return raw

    def cartesian(
        self,
        frame_orientation: str = "ECLIPJ2000",
        aberations: str = "geometric",
    ) -> np.ndarray:
        """Retrieve the cartesian state using the Horizons Vector API.

        Parameters
        ----------
        frame_orientation : str, optional
            Reference Frame Orientation, equivalent to the astroquery refplane
            parameter. Options are 'J2000'/'earth', 'ECLIPJ2000'/'ecliptic' and 'body'
            , by default "ECLIPJ2000".
            See astroquery documentation for information about the 'body' option:

            https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html
        aberations : str, optional
            Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
            'apparent', by default "geometric".

            See the Horizons System Manual for more info:

            https://ssd.jpl.nasa.gov/horizons/manual.html#output

        Returns
        -------
        np.ndarray
            returns an n by 7 array with the time in seconds since J2000 TDB,
            and the cartesian position and velocities.
        """
        raw = self.vectors(frame_orientation=frame_orientation, aberations=aberations)

        # A.D. 2019-Jan-05 22:40:00.0000
        timeformatt = "A.D. %Y-%b-%d %H:%M:%S.%f"

        tab = (
            raw.to_pandas()
            # format time: first parse the time string and then into seconds since J2000
            .assign(
                epoch_dt=lambda x: pd.to_datetime(x.datetime_str, format=timeformatt)
            )
            .assign(
                epochJ2000secondsTDB=lambda x: (
                    (
                        Time(x.epoch_dt, format="datetime64").jd1
                        - constants.JULIAN_DAY_ON_J2000
                    )
                    * constants.JULIAN_DAY
                )
                + ((Time(x.epoch_dt, format="datetime64").jd2) * constants.JULIAN_DAY)
            )
            .assign(x=lambda i: i.x * constants.ASTRONOMICAL_UNIT)
            .assign(y=lambda i: i.y * constants.ASTRONOMICAL_UNIT)
            .assign(z=lambda i: i.z * constants.ASTRONOMICAL_UNIT)
            .assign(
                vx=lambda i: i.vx * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY
            )
            .assign(
                vy=lambda i: i.vy * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY
            )
            .assign(
                vz=lambda i: i.vz * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY
            )
            .loc[:, ["epochJ2000secondsTDB", "x", "y", "z", "vx", "vy", "vz"]]
        )

        return tab.to_numpy()

    def create_ephemeris_tabulated(
        self,
        frame_origin: str,
        frame_orientation: str = "ECLIPJ2000",
        aberations: str = "geometric",
    ) -> environment_setup.ephemeris.EphemerisSettings:
        """Create ephemeris settings for a body using Horizons Vector API.

        Parameters
        ----------
        frame_origin : str
            Global frame origin, should match the queries' location parameter.
        frame_orientation : str, optional
            Reference Frame Orientation, equivalent to the astroquery refplane
            parameter. Options are 'J2000' and 'ECLIPJ2000', by default "ECLIPJ2000".
        aberations : str, optional
            Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
            'apparent', by default "geometric".

            See the Horizons System Manual for more info:

            https://ssd.jpl.nasa.gov/horizons/manual.html#output

        Returns
        -------
        environment_setup.ephemeris.EphemerisSettings
            Ephemeris settings for the query's body.

        Examples
        ----------
        Add Ephemerides of JUICE to the body_settings

        >>> body_settings.add_empty_settings("JUICE")
        >>> body_settings.get("JUICE").ephemeris_settings = query.create_ephemeris_tabulated(
                frame_origin=global_frame_origin,
                frame_orientation=global_frame_orientation,
            )
        """
        if frame_orientation not in ["ECLIPJ2000", "J2000"]:
            raise ValueError(
                "refplane parameter must be one of: " + '"ECLIPJ2000", "J2000"'
            )
        vector = self.cartesian(
            frame_orientation=frame_orientation, aberations=aberations
        )

        table = {x[0]: x[1:7] for x in vector}

        return environment_setup.ephemeris.tabulated(
            body_state_history=table,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
        )

    def ephemerides(
        self,
        reference_system: str = "J2000",
        extra_precision: bool = False,
        *args,
        **kwargs,
    ) -> astropy.table.Table:
        """Implements the JPL Horizons ephemerides API and returns it in raw Astropy table format.
        Ephemerides API provides time-interpolated observer parameters such as right ascension and declination.
        Note that this means that values provided are not actual observations.

        A number of quantities are retrieved, their definitions can be found here:
        https://ssd.jpl.nasa.gov/horizons/manual.html#obsquan.
        By default all available quantities are retrieved.

        More parameters can be passed directly to the astroquery call. These can be passed as kwargs: kwargs=("refraction":True).
        Check the astroquery documentation for an overview:
        https://astroquery.readthedocs.io/en/latest/api/astroquery.jplhorizons.HorizonsClass.html#astroquery.jplhorizons.HorizonsClass.ephemerides


        Parameters
        ----------
        reference_system : str, optional
            Coordinate reference system, value must be one of `ICRF`/`J2000` or B1950, by default "J2000"
        extra_precision : bool, optional
            Enables extra precision in right ascension and declination values, by default False

        Returns
        -------
        astropy.table.Table
            Unprocessed output in astropy table format.

        Raises
        ------
        ValueError
            If time query has incorrect format or an incorrect reference system is chosen
        """
        if reference_system not in ["ICRF", "J2000", "B1950"]:
            raise ValueError(
                "`reference_system` must be one of: `J2000`(=`ICRF`), `ICRF`, `B1950`"
            )
        if reference_system == "J2000":
            reference_system = "ICRF"

        res_list = []
        for query in self.queries:
            # EPHEMERIDES ONLY TAKES UT TIME (UTC)
            # temporarily transform from TDB to UTC:
            # NOTE the conversion back and forth is accurate to 1e-5 seconds

            if isinstance(query.epochs, dict):
                query.epochs["start"] = Time(
                    query.epochs["start"], format="iso", scale="tdb"
                ).utc.iso
                query.epochs["stop"] = Time(
                    query.epochs["stop"], format="iso", scale="tdb"
                ).utc.iso
            elif isinstance(query.epochs, list) or isinstance(query.epochs, np.ndarray):
                query.epochs = [
                    (
                        Time(x, format="jd", scale="tdb").utc.jd1
                        + Time(x, format="jd", scale="tdb").utc.jd2
                    )
                    for x in query.epochs
                ]
            else:
                raise ValueError("query epoch has incorrect format")

            res = query.ephemerides(
                refsystem=reference_system,
                extra_precision=extra_precision,
                *args,
                **kwargs,
            )

            timeformatt = "%Y-%b-%d %H:%M:%S.%f"
            actual_time_tdb = Time.strptime(
                res["datetime_str"], format_string=timeformatt, scale="utc"
            ).tdb

            res["datetime_str"] = actual_time_tdb.strftime(timeformatt)

            new_actual_time_tdb = Time.strptime(
                res["datetime_str"], format_string=timeformatt
            )

            res["datetime_jd"] = new_actual_time_tdb.jd1 + new_actual_time_tdb.jd2
            res["epochJ2000secondsTDB"] = (
                (new_actual_time_tdb.jd1 - constants.JULIAN_DAY_ON_J2000)
                * constants.JULIAN_DAY
            ) + (new_actual_time_tdb.jd2 * constants.JULIAN_DAY)

            res_list.append(res)

            # convert back
            if isinstance(query.epochs, dict):
                query.epochs["start"] = Time(
                    query.epochs["start"], format="iso", scale="utc"
                ).tdb.iso
                query.epochs["stop"] = Time(
                    query.epochs["stop"], format="iso", scale="utc"
                ).tdb.iso
            elif isinstance(query.epochs, list):
                query.epochs = [
                    Time(x, format="jd", scale="utc").tdb.jd for x in query.epochs
                ]
            else:
                raise ValueError("query epoch has incorrect format")

        raw = vstack(res_list)

        return raw

    def interpolated_observations(
        self,
        degrees: bool = False,
        reference_system: str = "J2000",
        extra_precision: bool = True,
        *args,
        **kwargs,
    ) -> np.ndarray:
        """Retrieves interpolated Right Ascension and Declination from the Horizons ephemerides API.
        Note that these values are not real observations but instead interpolated
        values based on the Horizons ephemeris system.

        Parameters
        ----------
        degrees : bool, optional
            return values in degrees if True, radians if False, by default false
        reference_system : str, optional
            Coordinate reference system, value must be one of `ICRF`/`J2000` or B1950, by default "J2000"
        extra_precision : bool, optional
            Enables extra precision in Right Ascension and Declination values, by default False

        Returns
        -------
        np.ndarray
            Numpy array (N, 3) with time in seconds since J2000 TDB and the Right Ascension and Declination.

        Raises
        ------
        ValueError
            If time query has incorrect format or an incorrect reference system is chosen
        """

        kwargs["quantities"] = 1  # this gets only RA + DEC
        raw = self.ephemerides(
            reference_system=reference_system,
            extra_precision=extra_precision,
            *args,
            **kwargs,
        )

        res = raw.to_pandas().loc[:, ["epochJ2000secondsTDB", "RA", "DEC"]]

        if not degrees:
            res[["RA", "DEC"]] = res[["RA", "DEC"]].apply(np.radians)

        res = res.to_numpy()

        return res


class HorizonsBatch:
    def __init__(
        self,
        query_id_list: List[str],
        location: str,
        epoch_list: Union[list, None] = None,
        epoch_start: Union[datetime.datetime, float, None] = None,
        epoch_end: Union[datetime.datetime, float, None] = None,
        epoch_step: Union[str, None] = None,
        extended_query: bool = False,
    ) -> None:
        """Query object to retrieve Horizons data for multiple bodies.
        See documentation for `HorizonsQuery` for more extensive documentation.
        Note that the epochs requested must have data for all bodies queried.
        This class is useful for quickly creating ephemerides for many objects at once.
        For general purposes, use `HorizonsQuery instead`.

        Parameters
        ----------
        query_id : str
            List of query terms to retrieve data for,
            all queries behave as type default and the `query_type` behaviour
            can not be set for the `HorizonsBatch` class.
        location : str
            Coordinate centre for the data with syntax `site@body`.
        epoch_start : Union[datetime.datetime, float, None], optional
            Starting date to retrieve data for.
            Must be either a python datetime object or a float seconds since J2000 TDB.
            Combined with `epoch_end` and
            `epoch_step` can be used to retrieve a range of data.
            If `epoch_list` is used, value must be None, by default None
        epoch_end : Union[datetime.datetime, float, None], optional
            Final date to retrieve data for.
            Must be either a python datetime object or a float seconds since J2000 TDB.
            If the start stop, and step parameters dont result in an
            integer multiple of values,the final date is not used.
            Combined with `epoch_end` and
            `epoch_step` can be used to retrieve a range of data.
            If `epoch_list` is used, value must be None, by default None
        epoch_step : Union[str, None], optional
            Step for the range of epochs to retrieve data for.
            Can be be either a specific timestep such as `15m`,
            or a number of partitions `10` for 10 equidistant steps within the range.
            For the timestep, quantifiers `m` for minute, `h` for hour and `d` for day
            can be used. Seconds are not available, consider using the `epoch_list`
            parameter or a number of partitions instead.
            Month and Year are normally available in JPL Horizons but are restricted
            here because they may produce ambiguous results (leap years etc.).
            If `epoch_list` is used, value must be None, by default None
        epoch_list : Union[list, None], optional
            List of times in seconds since J2000 TDB. Can be used
            to retrieve specific times instead of a range.
            Must be None if start, end and step are set, by default None
        extended_query : bool, optional
            Enables the retrieval of larger collections of data, by default False.
        """
        self._query_objects = {}
        self._query_id_list = query_id_list
        self._names = []

        for query_id in query_id_list:
            if not isinstance(query_id, str):
                raise TypeError("Ids in query_id_list must be of type str")
            temp = HorizonsQuery(
                query_id=query_id,
                location=location,
                epoch_list=epoch_list,
                epoch_start=epoch_start,
                epoch_end=epoch_end,
                epoch_step=epoch_step,
                extended_query=extended_query,
            )

            self._query_objects[query_id] = temp

    @property
    def names(self) -> Union[None, List[str]]:
        """Retrieves a list of names of the query objects. Returns `None`
        if `add_batch_ephemerides` has not been run yet."""
        return self._names

    def add_batch_ephemerides(
        self,
        body_settings: environment_setup.BodyListSettings,
        frame_origin: str,
        frame_orientation: str = "ECLIPJ2000",
        aberations: str = "geometric",
    ) -> None:
        """Uses the data queried to add ephemerides of the bodies querried
        to the body_settings. The names of the bodies added can be retrieved
        using the names property.

        Parameters
        ----------
        body_settings : environment_setup.BodyListSettings
            Tudat body settings object.
        frame_origin : str
            Global frame origin, should match the queries' location parameter.
        frame_orientation : str, optional
            Reference Frame Orientation, equivalent to the astroquery refplane
            parameter. Options are 'J2000' and 'ECLIPJ2000', by default "ECLIPJ2000".
        aberations : str, optional
            Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
            'apparent', by default "geometric".

        See the Horizons System Manual for more info:

        https://ssd.jpl.nasa.gov/horizons/manual.html#output

        """

        names = []

        # add the body and its ephemerides to the body_settings
        for query in self._query_objects.values():
            eph = query.create_ephemeris_tabulated(
                frame_origin=frame_origin,
                frame_orientation=frame_orientation,
                aberations=aberations,
            )
            name = query.name
            names.append(name)

            body_settings.add_empty_settings(name)
            body_settings.get(name).ephemeris_settings = eph

        # retrieve the names to a list
        self._names = names
