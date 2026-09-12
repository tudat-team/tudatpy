"""Legacy TNF converters, including their station-dependent time conversion.

Retained from develop at 6dfa06e689 for the deprecation period. The new
tracking-data converters remain in data_input.tracking_data.tnf.
"""

from datetime import datetime, timedelta

import numpy as _np
from pandas import DataFrame
from trk234 import bands, SFDU

from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data.tnf._converters import Converter
from tudatpy.dynamics.environment_setup.ground_station import (
    get_approximate_dsn_ground_station_positions,
)
from tudatpy.estimation.observable_models_setup import links
from tudatpy.estimation.observable_models_setup.links import link_definition, receiver
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.estimation.observations import create_single_observation_set, SingleObservationSet
from tudatpy.estimation.observations_setup.ancillary_settings import (
    FrequencyBands,
    dsn_n_way_doppler_ancillary_settings,
    dsn_n_way_range_ancillary_settings,
)


class RadioBase(Converter):

    trkModeDict = {
        0: "Unknown",
        1: "1W",
        2: "2W",
        3: "3W",
    }

    time_scale_converter = time_representation.default_time_scale_converter()

    frequencyBandsDict = {
        "S": FrequencyBands.s_band,
        "X": FrequencyBands.x_band,
        # "K": FrequencyBands.ku_band,
        "Ka": FrequencyBands.ka_band,
    }

    stationDict = get_approximate_dsn_ground_station_positions()

    def get_link_ends(self, sfdu: SFDU) -> tuple[str, str, str]:
        """
        Returns the uplink, spacecraft, and downlink IDs for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the link ends from.

        Returns
        -------
        tuple(str, str, str)
            A tuple containing the uplink, spacecraft, and downlink IDs.
            If the uplink is unknown or not valid, the uplink entry is a `"nan"` string.
        """
        upLink = (
            sfdu.sec_chdo.vld_ul_stn
            if sfdu.sec_chdo.vld_ul_stn != 0
            else sfdu.sec_chdo.ul_prdx_stn if sfdu.sec_chdo.ul_prdx_stn != 0 else "nan"
        )
        upLink = "DSS-" + str(upLink) if upLink != "nan" else upLink

        # Add minus sign to comply with NAIF convention
        scId = str(-sfdu.sec_chdo.scft_id)
        dlLink = "DSS-" + str(sfdu.sec_chdo.dl_dss_id)

        return (upLink, scId, dlLink)

    def get_band(self, sfdu: SFDU) -> tuple[str, str]:
        """
        Returns the uplink and downlink radio bands for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the radio bands from.

        Returns
        -------
        tuple(str, str)
            A tuple containing the uplink and downlink radio bands.
        """
        return (
            bands[sfdu.sec_chdo.ul_band_dl],
            bands[sfdu.sec_chdo.vld_dl_band],
        )

    def get_tracking_mode(self, sfdu: SFDU) -> str:
        """
        Returns the tracking mode for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the tracking mode from.

        Returns
        -------
        str
            The tracking mode of the SFDU record.
        """
        trkMode = (
            self.trkModeDict[sfdu.sec_chdo.vld_dop_mode]
            if sfdu.sec_chdo.vld_dop_mode != 0
            else sfdu.tracking_mode()
        )
        return trkMode

    def build_link_ends_dict(
        self, link_end_tuple: tuple[str, str, str], spacecraftName: str | None = None
    ) -> dict:
        """
        Construct a link ends dictionary for Doppler/Range observation creation.

        Parameters
        ----------
        link_end_tuple : tuple[str, str, str]
            A tuple containing the uplink, spacecraft, and downlink identifiers.
        spacecraftName : str, optional
            The name of the spacecraft to use in the simulation. If not provided, the spacecraft name
            is extracted from the TNF file.

        Returns
        -------
        dict
            A dictionary with link ends constructed according to the following logic:
            - If the uplink identifier (first element) is "nan", use the spacecraft as the transmitter
            and assign the downlink using Earth's reference point (from the third element).
            - Otherwise, assign the transmitter using Earth's reference point from the first element,
            the reflector as the spacecraft, and the receiver using Earth's reference point from the third element.
        """

        if len(link_end_tuple) != 3:
            raise ValueError(
                "Error when processing TNF file, building link ends dictionary: \n"
                + f"the link end tuple should contain exactly 3 elements: {link_end_tuple} provided."
            )

        # Set custom spacecraft name if provided
        if spacecraftName is not None:
            spacecraft = links.body_origin_link_end_id(spacecraftName)
        else:
            spacecraft = links.body_origin_link_end_id(link_end_tuple[1])

        if link_end_tuple[0] == "nan":
            return {
                links.transmitter: spacecraft,
                links.receiver: links.body_reference_point_link_end_id("Earth", link_end_tuple[2]),
            }
        else:
            return {
                links.transmitter: links.body_reference_point_link_end_id(
                    "Earth", link_end_tuple[0]
                ),
                links.reflector1: spacecraft,
                links.receiver: links.body_reference_point_link_end_id("Earth", link_end_tuple[2]),
            }

    def from_datetime_UTC_to_TDB(self, datetime_utc: datetime, station: str) -> float:
        """
        Convert a datetime object in UTC into seconds since J2000 in TDB.

        Parameters
        ----------
        datetime_utc : datetime
            The datetime object to convert.

        Returns
        -------
        float
            The time in seconds since J2000 in TDB.
        """
        if station not in self.stationDict:
            raise KeyError(
                "Error when processing TNF file, converting time from UTC to TDB: \n"
                + "the position of the ground station {} was not specified.".format(station)
            )

        epoch_utc = time_representation.DateTime.from_python_datetime(datetime_utc).to_epoch()
        epoch_tdb = self.time_scale_converter.convert_time(
            input_scale=time_representation.utc_scale,
            output_scale=time_representation.tdb_scale,
            input_value=epoch_utc,
            earth_fixed_position=self.stationDict[station],
        )

        return epoch_tdb


class DerivedDopplerConverter(RadioBase):
    def extract(self, sfdu_list: list[SFDU]) -> DataFrame:
        # Filter SFDU objects that represent derived Carrier Doppler data.
        # - Derived Carrier Doppler format_code == 16
        # - Only keep decoded ones
        # - Only keep 2W and 3W tracking mode
        doppler_sfdu = [sfdu for sfdu in sfdu_list if sfdu.pri_chdo.format_code == 16]
        for sfdu in doppler_sfdu:
            sfdu.decode(sfdu.binarydata, label=False, agg_chdo=False, pri_chdo=False)
        doppler_sfdu = [
            sfdu
            for sfdu in doppler_sfdu
            if sfdu.is_decoded and (self.get_tracking_mode(sfdu) not in ("None", "1W"))
        ]

        data = {
            "epoch": [sfdu.timestamp() for sfdu in doppler_sfdu],
            "link_ends": [self.get_link_ends(sfdu) for sfdu in doppler_sfdu],
            "band": [self.get_band(sfdu) for sfdu in doppler_sfdu],
            "tracking_mode": [self.get_tracking_mode(sfdu) for sfdu in doppler_sfdu],
            "link_delays": [self.get_link_delays(sfdu) for sfdu in doppler_sfdu],
            "count_time": [sfdu.trk_chdo.obs_cnt_time for sfdu in doppler_sfdu],
            "obs": [sfdu.trk_chdo.rcv_carr_obs[0] for sfdu in doppler_sfdu],
        }

        return DataFrame(data)

    def process(
        self, doppler_df: DataFrame, spacecraftName: str | None = None
    ) -> list[SingleObservationSet]:

        observation_set_list = []
        for link_end in doppler_df["link_ends"].unique():
            link_ends_dict = self.build_link_ends_dict(link_end, spacecraftName)
            link_def = link_definition(link_ends_dict)
            df_le = doppler_df[doppler_df["link_ends"] == link_end]
            for band in df_le["band"].unique():
                df_band = df_le[df_le["band"] == band]
                for ttd in df_band["link_delays"].unique():
                    df_ttd = df_band[df_band["link_delays"] == ttd]
                    for ct in df_ttd["count_time"].unique():
                        df_ct = df_ttd[df_ttd["count_time"] == ct]
                        ancillary_settings = dsn_n_way_doppler_ancillary_settings(
                            [
                                self.frequencyBandsDict[band[0]],
                                self.frequencyBandsDict[band[1]],
                            ],
                            self.frequencyBandsDict[band[1]],
                            0.0,
                            ct,
                            ttd,
                        )
                        obs_values = df_ct["obs"].to_numpy(dtype=float).reshape((-1, 1))
                        station = link_end[2] if len(link_end) == 3 else link_end[1]
                        epoch_seconds = (
                            df_ct["epoch"]
                            .apply(lambda t: self.from_datetime_UTC_to_TDB(t, station))
                            .tolist()
                        )
                        observation_set = create_single_observation_set(
                            ObservableType.dsn_n_way_averaged_doppler_type,
                            link_def.link_ends,
                            obs_values,
                            epoch_seconds,
                            receiver,
                            ancillary_settings,
                        )
                        observation_set_list.append(observation_set)

        return observation_set_list

    def get_link_delays(self, sfdu: SFDU) -> tuple[float, float, float]:
        """
        Returns the transmit time tag delay, spacecraft transmit delay, and receive time tag delay for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the time tag delays from.

        Returns
        -------
        tuple(float, float, float)
            A tuple containing the transmit time tag delay, spacecraft transmit delay, and receive time tag delay.
            If the values are not valid or not provided, the function sets the delays to 0.
        """
        uplinkDelay = 0.0
        uplinkDelay += (
            sfdu.sec_chdo.transmit_time_tag_delay
            if sfdu.sec_chdo.transmit_time_tag_delay != -1.0
            else 0.0
        )
        uplinkDelay += (
            sfdu.sec_chdo.ul_zheight_corr if sfdu.sec_chdo.ul_zheight_corr != -99.0 else 0.0
        )

        downlinkDelay = 0.0
        downlinkDelay += (
            sfdu.sec_chdo.rcv_time_tag_delay if sfdu.sec_chdo.rcv_time_tag_delay != -1.0 else 0.0
        )
        downlinkDelay += sfdu.sec_chdo.array_delay if sfdu.sec_chdo.array_flag != 0.0 else 0.0
        downlinkDelay += (
            sfdu.sec_chdo.dl_zheight_corr if sfdu.sec_chdo.dl_zheight_corr != -99.0 else 0.0
        )

        scft_transpd_delay = (
            sfdu.sec_chdo.scft_transpd_delay if sfdu.sec_chdo.scft_transpd_delay != -1.0 else 0.0
        )

        return (uplinkDelay, scft_transpd_delay, downlinkDelay)


class DerivedSraRangeConverter(RadioBase):
    def extract(self, sfdu_list: list[SFDU]) -> DataFrame:
        # Filter SFDU objects that represent SRA Range data.
        # - SRA Range format_code == 7
        # - Only keep decoded ones
        # - Skip records with invalid obs (obs == -1.0 and rng_type != 0)
        # - Skip if uplink frequency is zero
        # - Only keep 2W and 3W tracking mode
        range_sfdu_list = [sfdu for sfdu in sfdu_list if sfdu.pri_chdo.format_code == 7]
        for sfdu in range_sfdu_list:
            sfdu.decode(sfdu.binarydata, label=False, agg_chdo=False, pri_chdo=False)
        range_sfdu = [
            sfdu
            for sfdu in range_sfdu_list
            if sfdu.is_decoded
            and (self.get_tracking_mode(sfdu) == "2W" or self.get_tracking_mode(sfdu) == "3W")
            and not (sfdu.trk_chdo.meas_rng == -1.0 and sfdu.trk_chdo.rng_type != 0)
            and sfdu.trk_chdo.ul_freq != 0.0
        ]

        data = {
            "epoch": [sfdu.timestamp() for sfdu in range_sfdu],
            "link_ends": [self.get_link_ends(sfdu) for sfdu in range_sfdu],
            "band": [self.get_band(sfdu) for sfdu in range_sfdu],
            "tracking_mode": [self.get_tracking_mode(sfdu) for sfdu in range_sfdu],
            "link_delays": [self.get_link_delays(sfdu) for sfdu in range_sfdu],
            "obs": [
                _np.mod(sfdu.trk_chdo.meas_rng, sfdu.trk_chdo.rng_modulo) for sfdu in range_sfdu
            ],
            "zero_phase_times": [self.get_zero_phase_times(sfdu) for sfdu in range_sfdu],
            "lowest_ranging_component": [sfdu.trk_chdo.last_comp_num for sfdu in range_sfdu],
        }

        return DataFrame(data)

    def process(
        self, range_df: DataFrame, spacecraftName: str | None = None
    ) -> list[SingleObservationSet]:

        observation_set_list = []
        for link_end in range_df["link_ends"].unique():
            link_ends_dict = self.build_link_ends_dict(link_end, spacecraftName)
            link_def = link_definition(link_ends_dict)
            df_le = range_df[range_df["link_ends"] == link_end]
            for band in df_le["band"].unique():
                df_band = df_le[df_le["band"] == band]
                for ttd in df_band["link_delays"].unique():
                    df_ttd = df_band[df_band["link_delays"] == ttd]
                    for lrc in df_ttd["lowest_ranging_component"].unique():
                        df_lrc = df_ttd[df_ttd["lowest_ranging_component"] == lrc]
                        # NOTE - we force to build an observation set per observable in order to
                        # store conversion factor calculated in the dsnNWayRangeObservationModel.h
                        # A faster approach would be to use the dataframe apply method as done in
                        # the derivedDoppler converter.
                        for _, row in df_lrc.iterrows():
                            ancillary_settings = dsn_n_way_range_ancillary_settings(
                                [
                                    self.frequencyBandsDict[band[0]],
                                    self.frequencyBandsDict[band[1]],
                                ],
                                lrc,
                                ttd,
                            )
                            obs_values = [_np.array([row["obs"]], dtype=float).reshape((-1, 1))]
                            station = link_end[2] if len(link_end) == 3 else link_end[1]
                            epoch_seconds = [self.from_datetime_UTC_to_TDB(row["epoch"], station)]
                            observation_set = create_single_observation_set(
                                ObservableType.dsn_n_way_range_type,
                                link_def.link_ends,
                                obs_values,
                                epoch_seconds,
                                receiver,
                                ancillary_settings,
                            )
                            observation_set_list.append(observation_set)

        return observation_set_list

    def get_link_delays(self, sfdu: SFDU) -> tuple[float, float, float]:
        """
        Returns the transmit time tag delay, spacecraft transmit delay, and receive time tag delay for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function. If the ruToSeconds parameter is provided,

        There are several range calibration and delay fields in the records. The most crucial one is stn_cal, which represents
        the round-trip measured range calibration for a single pass. If the measured calibration has been split into two parts,
        it may appear in dl_stn_cal and ul_stn_cal.

        stn_cal contains the uplink, downlink, and array_delay calibrations. The array_delay is a delay that applies only to the downlink.
        If the calibration is split, it will be part of the dl_stn_cal field. The other calibration fields, rcv_timetag_delay and
        transmit_time_tag_delay, are nominal values and should be disregarded.

        Tudat uses range calibrations as time delays. The calibration is divided into rcvr (downlink) and xmtr (uplink) delays,
        which are set into each measurement.

        Since the stn_cal includes the uplink and downlink calibrations (which are identical for two-way data), we can remove
        the array_delay from stn_cal and divide it by two to obtain the individual leg calibration.
        Then, we can add the array_delay back into the downlink.

        The final step in the calibration is the z-height delays. These are the last leg of the electronics in the antenna and are not
        included in the calibration electronics. To obtain the final calibration values, we need to add the uplink and
        downlink z-height delays to each leg.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the time tag delays from.
        ruToSeconds : float, optional
            The conversion factor to convert the station calibration from range unit to seconds. Default is None.

        Returns
        -------
        tuple(float, float, float)
            A tuple containing the transmit time tag delay, spacecraft transmit delay, and receive time tag delay.
            If the values are not valid or not provided, the function sets the delays to 0.
        """
        ruToSeconds = (sfdu.trk_chdo.exc_scalar_den / sfdu.trk_chdo.exc_scalar_num) / (
            16 * sfdu.trk_chdo.ul_freq
        )

        uplinkDelay = 0.0
        uplinkDelay += (
            sfdu.trk_chdo.ul_stn_cal if sfdu.trk_chdo.ul_stn_cal != -1.0 else 0.0
        ) * ruToSeconds
        uplinkDelay += (
            sfdu.sec_chdo.ul_zheight_corr if sfdu.sec_chdo.ul_zheight_corr != -99.0 else 0.0
        )

        downlinkDelay = 0.0
        downlinkDelay += (
            sfdu.trk_chdo.dl_stn_cal if sfdu.trk_chdo.dl_stn_cal != -1.0 else 0.0
        ) * ruToSeconds
        downlinkDelay += (
            sfdu.sec_chdo.dl_zheight_corr if sfdu.sec_chdo.dl_zheight_corr != -99.0 else 0.0
        )

        scft_transpd_delay = (
            sfdu.sec_chdo.scft_transpd_delay if sfdu.sec_chdo.scft_transpd_delay != -1.0 else 0.0
        )

        return (uplinkDelay, scft_transpd_delay, downlinkDelay)

    def get_zero_phase_times(self, sfdu: SFDU) -> tuple[datetime, datetime]:
        """
        Get the zero phase times from the three way SRA Range.

        Parameters
        ----------
        sfdu : SFDU
            The SFDU object to extract the zero phase times from.

        Returns
        -------
        tuple
            A tuple containing the transmit and receive zero phase times.
        """

        return (
            sfdu.timestamp() + timedelta(seconds=sfdu.trk_chdo.transmit_inphs_time),
            sfdu.timestamp() + timedelta(seconds=sfdu.trk_chdo.rcv_inphs_time),
        )
