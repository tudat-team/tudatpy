from astroquery.jplhorizons import Horizons
from typing import Union, List, Tuple
import datetime

from ..kernel import constants
from ..kernel.numerical_simulation import environment_setup


class HorizonsQuery:
    def __init__(
        self,
        query_id: str,
        location: str,
        query_type: str = "default",
        epoch_list: Union[float, list, None] = None,
        epoch_start: Union[float, datetime.datetime, None] = None,
        epoch_end: Union[float, datetime.datetime, None] = None,
        epoch_step: Union[float, str, None] = None,
    ) -> None:
        self.query_id = query_id
        self.query_type = None if query_type == "default" else query_type

        self._number = None
        self._name = None
        self._designation = None

        # error handling
        if epoch_list is None:
            if (epoch_start is None) or (epoch_end is None) or (epoch_step is None):
                raise ValueError(
                    "Must specify either a list of times in sec since J2000 \
                                 or a combined start, end and step parameters"
                )
            else:
                epoch_start = self._convert_time(epoch_start)
                epoch_end = self._convert_time(epoch_end)
                epoch_step = self._convert_timestep(epoch_step)

                epoch_def = dict(start=epoch_start, stop=epoch_end, step=epoch_step)

        elif (
            (epoch_start is not None)
            and (epoch_end is not None)
            and (epoch_step is not None)
        ):
            if epoch_list is not None:
                raise ValueError(
                    "Must specify either a list of times in sec since J2000 \
                                 or a combined start, end and step parameters not both"
                )
            if (
                isinstance(epoch_list, float)
                or isinstance(epoch_list, int)
                or isinstance(epoch_list, datetime.datetime)
            ):
                epoch_def = [self._convert_time(epoch_list)]
            elif isinstance(epoch_list, list) or isinstance(epoch_list, tuple):
                epoch_def = [self._convert_time(x) for x in epoch_list]
            else:
                raise ValueError(
                    "epoch_list must be a list of/ singular float int or datetime object"
                )
        else:
            raise ValueError(
                """Must specify either a list of times in sec since J2000 or a combined start, end and step parameters"""
            )

        self.body = Horizons(
            id=self.query_id,
            location=location,
            id_type=self.query_type,
            epochs=epoch_def,
            # *args,
            # **kwargs
        )

    @property
    def name(self):
        if self._name is None:
            return None
        elif "Spacecraft" in self._name:
            return self._number
        else:
            return self._name

    def _convert_time(self, time: Union[float, datetime.datetime]):
        # https://ssd.jpl.nasa.gov/horizons/manual.html#time
        if isinstance(time, float) or isinstance(time, int):
            # convert sec since J2000 TDB to JD TDB
            time = "JD" + str((time / 86400) + 2451545.0) + "TDB"
            return time
        elif isinstance(time, datetime.datetime):
            # using recommended 3 letter months may have compatibility issues with locales
            formatt = r"%Y-%m-%d %H:%M:%S.%f"
            return str(time.strftime(formatt))
        else:
            raise ValueError(
                "Incorrect time value given, must be datetime object or seconds since J2000"
            )

    def _convert_timestep(self, epoch_step):
        if isinstance(epoch_step, int) or isinstance(epoch_step, float):
            return str(epoch_step) + "s"
        elif isinstance(epoch_step, str):
            return epoch_step
        else:
            raise ValueError(
                "Incorrect time step value given, must be either float in seconds or str: '10h', '10m', '1d', etc."
            )

    def vectors(
        self,
        full_output=False,
        refplane: str = "ecliptic",
        aberations: str = "geometric",
        *args,
        **kwargs
    ):
        raw = self.body.vectors(
            refplane=refplane,
            aberrations=aberations,  # args=args, kwargs=kwargs
        )

        if len(raw["targetname"][0].split()) == 3:
            self._number, self._name, self._designation = raw["targetname"][0].split()
        else:
            self._name = raw["targetname"][0]
            self._number = raw["targetname"][0]

        if full_output:
            return raw

        tab = (
            raw.to_pandas()
            .assign(epochJ2000secondsTDB=lambda x: (x.datetime_jd - 2451545.0) * 86400)
            .assign(x=lambda i: i.x * constants.ASTRONOMICAL_UNIT)
            .assign(y=lambda i: i.y * constants.ASTRONOMICAL_UNIT)
            .assign(z=lambda i: i.z * constants.ASTRONOMICAL_UNIT)
            .assign(vx=lambda i: i.vx * constants.ASTRONOMICAL_UNIT / 86400)
            .assign(vy=lambda i: i.vy * constants.ASTRONOMICAL_UNIT / 86400)
            .assign(vz=lambda i: i.vz * constants.ASTRONOMICAL_UNIT / 86400)
            .loc[:, ["epochJ2000secondsTDB", "x", "y", "z", "vx", "vy", "vz"]]
        )

        return tab.to_numpy()

    def create_ephemeris_tabulated(
        self,
        frame_origin,
        frame_orientation,
        refplane: str = "ecliptic",
        aberations: str = "geometric",
        *args,
        **kwargs
    ):
        vector = self.vectors(
            refplane=refplane, aberations=aberations, args=args, kwargs=kwargs
        )

        table = {x[0]: x[1:7] for x in vector}

        return environment_setup.ephemeris.tabulated(
            body_state_history=table,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
        )
