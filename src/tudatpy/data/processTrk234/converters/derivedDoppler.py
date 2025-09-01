from tudatpy.estimation.observations_setup.ancillary_settings import dsn_n_way_doppler_ancilliary_settings
from tudatpy.estimation.observable_models_setup.links import link_definition, receiver, reflector1
from tudatpy.estimation.observable_models_setup.model_settings import dsn_n_way_doppler_averaged
from tudatpy.estimation.observations import single_observation_set
from . import RadioBase
from pandas import DataFrame


class DerivedDopplerConverter(RadioBase):
    def extract(self, sfdu_list):
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

    def process(self, doppler_df, spacecraftName=None):

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
                        ancillary_settings = dsn_n_way_doppler_ancilliary_settings(
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
                            .apply(lambda t: self.from_datetime_to_TBD(t, station))
                            .tolist()
                        )
                        observation_set = single_observation_set(
                            dsn_n_way_averaged_doppler_type,
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
            sfdu.sec_chdo.ul_zheight_corr
            if sfdu.sec_chdo.ul_zheight_corr != -99.0
            else 0.0
        )

        downlinkDelay = 0.0
        downlinkDelay += (
            sfdu.sec_chdo.rcv_time_tag_delay
            if sfdu.sec_chdo.rcv_time_tag_delay != -1.0
            else 0.0
        )
        downlinkDelay += (
            sfdu.sec_chdo.array_delay if sfdu.sec_chdo.array_flag != 0.0 else 0.0
        )
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
