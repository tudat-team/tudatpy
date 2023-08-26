# import os
import sys

sys.path.insert(0, r"/mnt/c/Users/Trez/Desktop/tudat-bundle/tudatpy/")
import pandas as pd
import numpy as np

# import astroquery.mpc.MPC as astroqueryMPC
from astroquery.mpc import MPC
import datetime
from typing import Union
import time

from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import estimation, estimation_setup
from tudatpy.kernel.numerical_simulation.estimation_setup import observation
from tudatpy.kernel.astro import element_conversion

from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from astropy import units as u

import matplotlib.pyplot as plt


class BatchMPC:
    def __init__(self) -> None:
        self._table: pd.DataFrame = pd.DataFrame()
        self.stations = []
        self.bands = []
        self.MPCcodes = []
        self.size = 0

        self.stationInfo = self._get_station_info()

        self._addObsCartPosition()

        self.epoch_start = 0
        self.epoch_end = 0

    # helper functions
    def _refresh_metadata(self):
        self.stations = list(self._table.observatory.unique())
        self.bands = list(self._table.band.unique())
        self.MPCcodes = list(self._table.number.unique())
        self.size = len(self._table)

        # NOTE consider making this the first and last value -> faster
        self.epoch_start = self._table.epochJ2000secondsTDB.min()
        self.epoch_end = self._table.epochJ2000secondsTDB.max()

    def _get_station_info(self):
        try:
            stationInfo = MPC.get_observatory_codes().to_pandas()
            return stationInfo
        except Exception as e:
            print("An error occured while retrieving observatory data")
            print(e)

            return None

    def _addObsCartPosition(self):
        temp = self.stationInfo

        # TODO replace with tudat constant:
        rConst = 6378137

        # https://github.com/matthewjohnpayne/MPCUtilities/blob/master/mpcutilities/obscode.py
        temp = (
            temp.assign(X=lambda x: x.cos * rConst * np.cos(np.radians(x.Longitude)))
            .assign(Y=lambda x: x.cos * rConst * np.sin(np.radians(x.Longitude)))
            .assign(Z=lambda x: x.sin * rConst)
        )
        self.stationInfo = temp

    # getters and setters for the table
    @property
    def table(self):
        return self._table

    @table.setter
    def table(self, value):
        self._table = value

    # methods for data retrievels
    def get_observations(self, MPCcodes: list):
        for code in MPCcodes:
            try:
                obs = MPC.get_observations(code)
                obs = obs.to_pandas()

                # convert JD to J2000
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

                obs["number"] = obs.number.astype(str)
                self._table = pd.concat([self._table, obs])

            except Exception as e:
                print(f"An error occured while retrieving {code}")
                print(e)

        self._refresh_metadata()

    def from_astroquery(
        self,
    ):
        # TODO add astroquery table manually
        # TODO validate to check if columns are present
        raise NotImplementedError()
        pass

    def from_pandas(
        self,
    ):
        # TODO same as above
        raise NotImplementedError()
        pass

    def filter(
        self,
        bands: Union[list, str] = None,
        stations: Union[list, str] = None,
        stations_exclude: Union[list, str] = None,
        epoch_start: Union[float, datetime.datetime] = None,
        epoch_end: Union[float, datetime.datetime] = None,
    ):
        assert not isinstance(
            stations, int
        ), "stations parameter must be of type 'str' or 'list'"

        if isinstance(bands, str):
            bands = [bands]
        if isinstance(stations, str):
            stations = [stations]
        if bands is not None:
            self._table = self._table.query("band == @bands")
        if stations is not None:
            self._table = self._table.query("observatory == @stations")
        if stations_exclude is not None:
            self._table = self._table.query("observatory != @stations_exclude")
        if epoch_start is not None:
            if isinstance(epoch_start, float) or isinstance(epoch_start, int):
                self._table = self._table.query("epochJ2000secondsTDB >= @epoch_start")
            elif isinstance(epoch_start, datetime.datetime):
                self._table = self._table.query("epochUTC >= @epoch_start")
        if epoch_end is not None:
            if isinstance(epoch_start, float) or isinstance(epoch_start, int):
                self._table = self._table.query("epochJ2000secondsTDB <= @epoch_start")
            elif isinstance(epoch_start, datetime.datetime):
                self._table = self._table.query("epochUTC <= @epoch_start")

        self._refresh_metadata()

    def to_tudat(self, bodies, station_body="Earth"):
        tempStations = self.stationInfo.query("Code == @self.stations").loc[
            :, ["Code", "X", "Y", "Z"]
        ]


        observations_table = pd.merge(
            left=self._table,
            right=tempStations,
            left_on="observatory",
            right_on="Code",
            how="left",
        )

        for body in self.MPCcodes:
            bodies.create_empty_body(str(body))

        for idx in range(len(tempStations)):
            ground_station_settings = environment_setup.ground_station.basic_station(
                station_name=tempStations.iloc[idx].Code,
                station_nominal_position=[
                    tempStations.iloc[idx].X,
                    tempStations.iloc[idx].Y,
                    tempStations.iloc[idx].Z,
                ],
            )

            # TODO CHECK IF THE BODy ALREADy EXISTS

            # Add the ground station to the environment
            environment_setup.add_ground_station(
                bodies.get_body(station_body), ground_station_settings
            )

        # TODO get unique combinations between links
        # Unique combinations
        unique_link_combos = (
            observations_table.loc[:, ["number", "observatory"]].drop_duplicates()
            # .assign(code=lambda x: x.number.astype(str) + "_" + x.observatory)
            # .set_index("code")
        ).values

        linksDictionary = {}
        observation_set_list = []
        for combo in unique_link_combos:
            MPC_number = combo[0]
            station_name = combo[1]
            # CREATE LINKS
            link_ends = dict()
            link_ends[observation.transmitter] = observation.body_origin_link_end_id(
                MPC_number
            )
            link_ends[
                observation.receiver
            ] = observation.body_reference_point_link_end_id(station_body, station_name)

            # TODO different link for satellite
            # TODO in this case also use actual satellite name and not the MPC code.
            # observation.body_origin_link_end_id(station_name)
            link_definition = observation.link_definition(link_ends)

            linksDictionary[f"{MPC_number}_{station_name}"] = link_definition

            observations_for_this_link = observations_table.query(
                "number == @MPC_number"
            ).query("observatory == @station_name")

            observation_angles = observations_for_this_link.loc[
                :, ["RA", "DEC"]
            ].to_numpy()


            observation_times = observations_for_this_link.loc[
                :, ["epochJ2000secondsTDB"]
            ].to_numpy()[:, 0]

            # TEST TRANSFORM TO GALACTIC?
            # observation_angles = SkyCoord(ra=observation_angles[:, 0], dec=observation_angles[:, 1], unit="rad",frame="fk4").transform_to("galactic").to_table().to_pandas().to_numpy()

            observation_set = estimation.single_observation_set(
                observation.angular_position_type,
                link_definition,
                observation_angles,
                observation_times,
                observation.receiver,
            )

            observation_set_list.append(observation_set)

        observation_collection = estimation.ObservationCollection(observation_set_list)
        return observation_collection, list(linksDictionary.values())
    
    def plotObservations(self):
        fig, ax = plt.subplots(1,1, subplot_kw={"projection":"aitoff"}, figsize =  (12, 6))

        a = plt.scatter((self._table.RA),( self._table.DEC), s =3, marker = "x", c=self._table.epochJ2000secondsTDB)
        plt.colorbar(mappable=a, ax=ax)
        ax.grid()
        return fig
    
    def summariseBatch(self):
        print("Batch Summary:")
        print(f"1. Batch includes {len(self.MPCcodes)} minor planets:")
        print(self.MPCcodes)
        # for u in self.MPCcodes:
        #     print(f"{u} | {self.table.query("num")}")
        # print(self.table)
        print(f"2. Batch includes {len(self.stations)} observatories:")
        # for u in self.stations:
            # print(f"{u} | {self.stationInfo.query("")}")
        print(self.stationInfo.query("Code == @self.stations").loc[:, ["Code", "Name", "X", "Y", "Z"]])
        print(f"3. Batch ranges from dates {self.table.epochUTC.min()} to {self.table.epochUTC.max()}")
        print(f"4. Batch includes {len(self.bands)} minor planets:")
        print(self.bands)
        print()

if __name__ == "__main__":
    pass
