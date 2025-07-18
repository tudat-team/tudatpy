from tudatpy.numerical_simulation.estimation_setup.observation import (
    dsn_n_way_range_ancilliary_settings,
    link_definition,
    receiver,
    dsn_n_way_range,
)
from tudatpy.numerical_simulation.estimation import single_observation_set

from . import RadioBase
from pandas import DataFrame
import numpy as _np
from datetime import timedelta


class DerivedSraRangeConverter(RadioBase):
    def extract(self, sfdu_list):
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
            and (
                self.get_tracking_mode(sfdu) == "2W"
                or self.get_tracking_mode(sfdu) == "3W"
            )
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
                _np.mod(sfdu.trk_chdo.meas_rng, sfdu.trk_chdo.rng_modulo)
                for sfdu in range_sfdu
            ],
            "zero_phase_times": [
                self.get_zero_phase_times(sfdu) for sfdu in range_sfdu
            ],
            "lowest_ranging_component": [
                sfdu.trk_chdo.last_comp_num for sfdu in range_sfdu
            ],
        }

        return DataFrame(data)

    def process(self, range_df, spacecraftName):

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
                            ancillary_settings = dsn_n_way_range_ancilliary_settings(
                                [
                                    self.frequencyBandsDict[band[0]],
                                    self.frequencyBandsDict[band[1]],
                                ],
                                lrc,
                                ttd,
                            )
                            obs_values = [
                                _np.array([row["obs"]], dtype=float).reshape((-1, 1))
                            ]
                            station = link_end[2] if len(link_end) == 3 else link_end[1]
                            epoch_seconds = [
                                self.from_datetime_to_TBD(row["epoch"], station)
                            ]
                            observation_set = single_observation_set(
                                dsn_n_way_range,
                                link_def,
                                obs_values,
                                epoch_seconds,
                                receiver,
                                ancillary_settings,
                            )
                            observation_set_list.append(observation_set)

        return observation_set_list

    def get_link_delays(self, sfdu):
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
            sfdu.sec_chdo.ul_zheight_corr
            if sfdu.sec_chdo.ul_zheight_corr != -99.0
            else 0.0
        )

        downlinkDelay = 0.0
        downlinkDelay += (
            sfdu.trk_chdo.dl_stn_cal if sfdu.trk_chdo.dl_stn_cal != -1.0 else 0.0
        ) * ruToSeconds
        downlinkDelay += (
            sfdu.sec_chdo.dl_zheight_corr
            if sfdu.sec_chdo.dl_zheight_corr != -99.0
            else 0.0
        )

        scft_transpd_delay = (
            sfdu.sec_chdo.scft_transpd_delay
            if sfdu.sec_chdo.scft_transpd_delay != -1.0
            else 0.0
        )

        return (uplinkDelay, scft_transpd_delay, downlinkDelay)

    def get_zero_phase_times(self, sfdu):
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
