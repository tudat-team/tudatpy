"""
Converters for processing ATDF data into structured tracking data.
"""

import abc
from datetime import datetime

import numpy as np
import pandas as pd
from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data import TrackingData


class AtdfConverter(abc.ABC):
    """Abstract base class for ATDF converters."""

    frequency_bands_mapping = {
        "S": "S-band",
        "X": "X-band",
        # "K": "Ku-band",
        "Ka": "Ka-band",
    }

    @abc.abstractmethod
    def extract(self, df: pd.DataFrame, sc_name: str) -> pd.DataFrame:
        """
        Extract data of interest from the decoded ATDF DataFrame.
        """
        pass

    @abc.abstractmethod
    def process(self, df: pd.DataFrame, earth_name: str = "Earth") -> list[TrackingData]:
        """
        Process an extracted DataFrame into Tudat structured format.
        """
        pass

    def build_link_ends(
        self, link_ends: tuple[str, str, str], earth_name: str
    ) -> list[tuple[tuple[str, str], str]]:
        """
        Construct a ``PlainLinkDefinition`` (list of ``((body, reference_point), link_end_role)``
        tuples) for TrackingData creation.
        """

        if link_ends[0] == "":
            return [
                ((link_ends[1], ""), "transmitter"),
                ((earth_name, link_ends[2]), "receiver"),
            ]

        return [
            ((earth_name, link_ends[0]), "transmitter"),
            ((link_ends[1], ""), "reflector_1"),
            ((earth_name, link_ends[2]), "receiver"),
        ]

    def atdf_station_to_tudat(self, station: str) -> str:
        return station.strip().replace(" ", "-")

    def get_link_ends(self, row: pd.Series, sc_name: str) -> tuple[str, str, str]:
        if row["Xmtr"] == "S/C":
            link_ends = ("", sc_name, self.atdf_station_to_tudat(row["Rcvr"]))

        else:
            link_ends = (
                self.atdf_station_to_tudat(row["Xmtr"]),
                sc_name,
                self.atdf_station_to_tudat(row["Rcvr"]),
            )

        return link_ends

    def get_link_end_delays(self, row: pd.Series) -> tuple[float, float, float]:

        xmtr_delay = row["XmtrDly (nsec)"] * 1e-9
        sc_delay = row["ScDly (nsec)"] * 1e-9
        rcvr_delay = row["RcvrDly (nsec)"] * 1e-9

        if row["Xmtr"] == "S/C":
            sc_delay += xmtr_delay

            xmtr_delay = 0.0

        return (xmtr_delay, sc_delay, rcvr_delay)

    def to_utc_epoch(self, datetime_utc: datetime) -> time_representation.Time:
        """
        Convert a datetime object in UTC into seconds since J2000 in UTC.

        Parameters
        ----------
        datetime_utc : datetime
            The datetime object to convert.

        Returns
        -------
        float
            The time in seconds since J2000 in UTC.
        """
        return time_representation.DateTime.from_python_datetime(
            datetime_utc
        ).to_epoch_time_object()


class AtdfRangeConverter(AtdfConverter):
    def extract(self, df: pd.DataFrame, sc_name: str) -> pd.DataFrame:

        df_nway_range = df[
            df["Data Type"].isin(
                [
                    "2-Way-Range",
                ]
            )
        ]

        data = {
            "epoch": df_nway_range["time_tag (DT)"].to_list(),
            "link_ends": [self.get_link_ends(row, sc_name) for _, row in df_nway_range.iterrows()],
            "band": [(row["UL"], row["DL"]) for _, row in df_nway_range.iterrows()],
            "link_end_delays": [
                self.get_link_end_delays(row) for _, row in df_nway_range.iterrows()
            ],
            "lowest_ranging_component": df_nway_range["Rng-LC"].to_list(),
            "obs": df_nway_range["Observed"].to_list(),
        }

        return pd.DataFrame(data)

    def process(self, df: pd.DataFrame, earth_name: str = "Earth") -> list[TrackingData]:

        tracking_data_list = []
        for link_end in df["link_ends"].unique():
            link_ends = self.build_link_ends(link_end, earth_name)
            df_le = df[df["link_ends"] == link_end]
            for band in df_le["band"].unique():
                df_band = df_le[df_le["band"] == band]
                for ttd in df_band["link_end_delays"].unique():
                    df_ttd = df_band[df_band["link_end_delays"] == ttd]
                    for lrc in df_ttd["lowest_ranging_component"].unique():
                        df_lrc = df_ttd[df_ttd["lowest_ranging_component"] == lrc]
                        # NOTE - we force to build a tracking data object per observable in order
                        # to store the conversion factor calculated in the
                        # dsnNWayRangeObservationModel.h. A faster approach would be to use the
                        # dataframe apply method as done in the doppler converter.
                        for _, row in df_lrc.iterrows():
                            obs_values = [np.array([row["obs"]], dtype=float).reshape((-1, 1))]
                            epoch_seconds = [self.to_utc_epoch(row["epoch"])]

                            tracking_data = TrackingData(
                                "DsnNWayRange",
                                link_ends,
                                obs_values,
                                epoch_seconds,
                                "receiver",
                                "UTC",
                            )
                            tracking_data.add_string_vector_ancillary_setting(
                                "frequency bands",
                                [
                                    self.frequency_bands_mapping[band[0]],
                                    self.frequency_bands_mapping[band[1]],
                                ],
                            )
                            tracking_data.add_double_ancillary_setting(
                                "DSN sequential range lowest ranging component",
                                float(lrc),
                            )
                            tracking_data.add_double_vector_ancillary_setting(
                                "link ends time delays", list(ttd)
                            )
                            tracking_data_list.append(tracking_data)

        return tracking_data_list


class AtdfNwayDopplerConverter(AtdfConverter):
    def extract(self, df: pd.DataFrame, sc_name: str) -> pd.DataFrame:

        df_nway_doppler = df[
            df["Data Type"].isin(
                [
                    "2-Way-Doppler",
                    "3-Way-Doppler",
                ]
            )
        ]

        data = {
            "epoch": df_nway_doppler["time_tag (DT)"].to_list(),
            "link_ends": [
                self.get_link_ends(row, sc_name) for _, row in df_nway_doppler.iterrows()
            ],
            "band": [(row["UL"], row["DL"]) for _, row in df_nway_doppler.iterrows()],
            "link_end_delays": [
                self.get_link_end_delays(row) for _, row in df_nway_doppler.iterrows()
            ],
            "count_time": df_nway_doppler["CT (sec)"].to_list(),
            "obs": df_nway_doppler["Observed"].to_list(),
            "ref_freq": df_nway_doppler["Ref-Freq (Hz)"].to_list(),
            "ex": df_nway_doppler["Ex"].to_list(),
        }

        return pd.DataFrame(data)

    def process(self, df: pd.DataFrame, earth_name: str = "Earth") -> list[TrackingData]:

        tracking_data_list = []

        for link_end in df["link_ends"].unique():
            link_ends = self.build_link_ends(link_end, earth_name)
            df_le = df[df["link_ends"] == link_end]
            for band in df_le["band"].unique():
                df_band = df_le[df_le["band"] == band]
                for ref_freq in df_band["ref_freq"].unique():
                    df_ref_freq = df_band[df_band["ref_freq"] == ref_freq]
                    for ttd in df_ref_freq["link_end_delays"].unique():
                        df_ttd = df_ref_freq[df_ref_freq["link_end_delays"] == ttd]
                        for ct in df_ttd["count_time"].unique():
                            df_ct = df_ttd[df_ttd["count_time"] == ct]
                            for ex in df_ct["ex"].unique():
                                df_ex = df_ct[df_ct["ex"] == ex]
                                obs_values = df_ex["obs"].to_numpy(dtype=float).reshape((-1, 1))
                                epoch_seconds = (
                                    df_ex["epoch"].apply(self.to_utc_epoch).to_numpy(dtype=float)
                                )

                                tracking_data = TrackingData(
                                    "DsnNWayAveragedDoppler",
                                    link_ends,
                                    obs_values,
                                    epoch_seconds,
                                    "receiver",
                                    "UTC",
                                )
                                tracking_data.add_string_vector_ancillary_setting(
                                    "frequency bands",
                                    [
                                        self.frequency_bands_mapping[band[0]],
                                        self.frequency_bands_mapping[band[1]],
                                    ],
                                )
                                tracking_data.add_string_ancillary_setting(
                                    "DSN reference frequency band at reception",
                                    self.frequency_bands_mapping[ex],
                                )
                                tracking_data.add_double_ancillary_setting(
                                    "DSN Doppler reference frequency", float(ref_freq)
                                )
                                tracking_data.add_double_ancillary_setting(
                                    "Doppler observable integration time", float(ct)
                                )
                                tracking_data.add_double_vector_ancillary_setting(
                                    "link ends time delays", list(ttd)
                                )
                                tracking_data_list.append(tracking_data)

        return tracking_data_list


class AtdfRampConverter:
    def extract(self, df: pd.DataFrame) -> pd.DataFrame:

        data = {
            "start_time": df["Start-Time (DT)"].to_list(),
            "end_time": df["End-Time (DT)"].to_list(),
            "station": df["Station"].apply(lambda x: x.strip().replace(" ", "-")).to_list(),
            "freq": df["Frequency (Hz)"].to_list(),
            "rate": df["Rate (Hz/sec)"].to_list(),
        }

        return pd.DataFrame(data)
