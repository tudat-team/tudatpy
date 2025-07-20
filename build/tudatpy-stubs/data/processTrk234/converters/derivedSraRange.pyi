import typing
from . import RadioBase as RadioBase
from tudatpy.estimation.observable_models_setup.links import link_definition as link_definition, receiver as receiver
from tudatpy.estimation.observable_models_setup.model_settings import dsn_n_way_range as dsn_n_way_range
from tudatpy.estimation.observations import single_observation_set as single_observation_set
from tudatpy.estimation.observations_setup.ancillary_settings import dsn_n_way_range_ancilliary_settings as dsn_n_way_range_ancilliary_settings

class DerivedSraRangeConverter(RadioBase):

    def extract(self, sfdu_list):
        ...

    def process(self, range_df, spacecraftName):
        ...

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