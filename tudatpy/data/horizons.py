from astroquery.jplhorizons import Horizons
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.table import vstack

# import tqdm.tqdm as tqdm

import math
import numpy as np
import pandas as pd

from typing import Union, List, Tuple
import datetime

from ..kernel import constants
from ..kernel.numerical_simulation import environment_setup

import re


class HorizonsQuery:
    query_limit = 90024

    # 50 is chosen to be conservative,
    # best result i have seen is 77 (for ~18 digits per epoch).
    # 50 is chosen to prevent scenarios where higher precision is used,
    # to prevent undefined behaviour when very futuristic dates are used,
    # and because it is a nice human friendly number.
    query_limit_list = 50

    def __init__(
        self,
        query_id: str,
        location: str,
        query_type: str = "default",
        epoch_list: Union[list, None] = None,
        epoch_start: Union[datetime.datetime, float, None] = None,
        epoch_end: Union[datetime.datetime, float, None] = None,
        epoch_step: Union[str, None] = None,
        extended_query: bool = False,
    ) -> None:
        self.query_id = query_id
        self.query_type = None if query_type == "default" else query_type

        self._number = None
        self._name = None
        self._designation = None

        self.epoch_type = None

        self.queries = []
        self.query_lengths = []

        # ##############
        # error handling
        # ##############

        # epoch range format:
        # epoch_list is none rest is not none
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
        # start step end is none, list nnot none
        elif (epoch_start is None) and (epoch_end is None) and (epoch_step is None):
            if epoch_list is None:
                raise ValueError(
                    "Must specify either a list of times in sec since J2000 "
                    + "or a combined start, end and step parameters not both"
                )
            elif isinstance(epoch_list, list) or isinstance(epoch_list, tuple):
                self.epoch_type = "list"
                num_lines = len(epoch_list)

            else:
                raise ValueError(
                    "epoch_list must be a list of times in sec since J2000"
                )
        # no format
        else:
            raise ValueError(
                "Must specify either a list of times in sec since J2000 "
                + "or a combined start, end and step parameters"
            )

        # #####################
        # Check for query limit
        # #####################
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
                    id=self.query_id,
                    location=location,
                    id_type=self.query_type,
                    epochs=epoch_def,
                    # *args,
                    # **kwargs
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
                            id=self.query_id,
                            location=location,
                            id_type=self.query_type,
                            epochs=list(split),
                            # *args,
                            # **kwargs
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
                            id=self.query_id,
                            location=location,
                            id_type=self.query_type,
                            epochs=epoch_def,
                            # *args,
                            # **kwargs
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
                        id=self.query_id,
                        location=location,
                        id_type=self.query_type,
                        epochs=epoch_def,
                        # *args,
                        # **kwargs
                    )
                )

            else:
                raise RuntimeError("Undefined behaviour, (extended)")

        else:
            raise RuntimeError("Undefined behaviour, (unknown error in init)")

        self.num_queries = len(self.queries)

    @property
    def name(self):
        if self._name is None:
            return None
        elif "Spacecraft" in self._name:
            return self._number
        else:
            return self._name

    def _convert_time_to_astropy(self, time):
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

    def _validate_time_range(self, start, end):
        if not (
            isinstance(start, float)
            or isinstance(start, int)
            or isinstance(start, datetime.datetime)
        ):
            raise TypeError(
                "Incorrect start time given, must be datetime object or float seconds since J2000 TDB"
            )
        if not (
            isinstance(end, float)
            or isinstance(end, int)
            or isinstance(end, datetime.datetime)
        ):
            raise TypeError(
                "Incorrect start time given, must be datetime object or float seconds since J2000 TDB"
            )

    def _format_time_range(self, time: Union[float, datetime.datetime]):
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

    def _format_time_list(self, times):
        times = self._convert_time_to_astropy((np.array(times))).jd
        return times

    def _interpret_timestep(
        self,
        timestep: str,
        start: Union[datetime.datetime, float],
        end: Union[datetime.datetime, float],
    ):
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
        full_output: bool = False,
        refplane: str = "ecliptic",
        aberations: str = "geometric",
        *args,
        **kwargs,
    ):
        res_list = []
        for query in self.queries:
            res = query.vectors(
                refplane=refplane,
                aberrations=aberations,  # args=args, kwargs=kwargs
            )

            res_list.append(res)

        raw = vstack(res_list)

        if full_output:
            return raw
        else:
            return raw.to_pandas()

        # if len(raw["targetname"][0].split()) == 3:
        #     self._number, self._name, self._designation = raw["targetname"][0].split()
        # else:
        #     self._name = raw["targetname"][0]
        #     self._number = raw["targetname"][0]

        # if full_output:
        #     return raw

    # def vectors(
    #     self,
    #     full_output=False,
    #     refplane: str = "ecliptic",
    #     aberations: str = "geometric",
    #     *args,
    #     **kwargs,
    # ):
    #     raw = self.query.vectors(
    #         refplane=refplane,
    #         aberrations=aberations,  # args=args, kwargs=kwargs
    #     )

    #     if len(raw["targetname"][0].split()) == 3:
    #         self._number, self._name, self._designation = raw["targetname"][0].split()
    #     else:
    #         self._name = raw["targetname"][0]
    #         self._number = raw["targetname"][0]

    #     if full_output:
    #         return raw

    #     tab = (
    #         raw.to_pandas()
    #         .assign(
    #             epochJ2000secondsTDB=lambda x: (
    #                 x.datetime_jd - constants.JULIAN_DAY_ON_J2000
    #             )
    #             * constants.JULIAN_DAY
    #         )
    #         .assign(x=lambda i: i.x * constants.ASTRONOMICAL_UNIT)
    #         .assign(y=lambda i: i.y * constants.ASTRONOMICAL_UNIT)
    #         .assign(z=lambda i: i.z * constants.ASTRONOMICAL_UNIT)
    #         .assign(
    #             vx=lambda i: i.vx * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY
    #         )
    #         .assign(
    #             vy=lambda i: i.vy * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY
    #         )
    #         .assign(
    #             vz=lambda i: i.vz * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY
    #         )
    #         .loc[:, ["epochJ2000secondsTDB", "x", "y", "z", "vx", "vy", "vz"]]
    #     )

    #     return tab.to_numpy()

    # def create_ephemeris_tabulated(
    #     self,
    #     frame_origin,
    #     frame_orientation,
    #     refplane: str = "ecliptic",
    #     aberations: str = "geometric",
    #     *args,
    #     **kwargs,
    # ):
    #     vector = self.vectors(
    #         refplane=refplane, aberations=aberations, args=args, kwargs=kwargs
    #     )

    #     table = {x[0]: x[1:7] for x in vector}

    #     return environment_setup.ephemeris.tabulated(
    #         body_state_history=table,
    #         frame_origin=frame_origin,
    #         frame_orientation=frame_orientation,
    #     )
