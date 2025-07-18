"""
Base class for radiometric converters sharing helper functions.
"""

from tudatpy.dynamics.environment_setup.ground_station import (
    # from tudatpy.dynamics.environment import (
    get_approximate_dsn_ground_station_positions,
)
from tudatpy.numerical_simulation.estimation_setup import observation  # type:ignore
from tudatpy.astro import time_conversion
from . import Converter
from trk234 import bands


class RadioBase(Converter):

    trkModeDict = {
        0: "Unknown",
        1: "1W",
        2: "2W",
        3: "3W",
    }

    time_scale_converter = time_conversion.default_time_scale_converter()

    frequencyBandsDict = {
        "S": observation.FrequencyBands.s_band,
        "X": observation.FrequencyBands.x_band,
        # "K": observation.FrequencyBands.ku_band,
        "Ka": observation.FrequencyBands.ka_band,
    }

    stationDict = get_approximate_dsn_ground_station_positions()

    def get_link_ends(self, sfdu):
        """
        Returns the uplink, spacecraft, and downlink IDs for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the link ends from.

        Returns
        -------
        tuple(int, int, int)
            A tuple containing the uplink, spacecraft, and downlink IDs.
            If the uplink or downlink are unknown or not valid, the function returns NaN.
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

    def get_band(self, sfdu):
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

    def get_tracking_mode(self, sfdu):
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

    def build_link_ends_dict(self, link_end_tuple, spacecraftName=None):
        """
        Construct a link ends dictionary for Doppler/Range observation creation.

        Parameters
        ----------
        link_end_tuple : tuple
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

        # Set custom spacecraft name if provided
        if spacecraftName is not None:
            spacecraft = observation.body_origin_link_end_id(spacecraftName)
        else:
            spacecraft = observation.body_origin_link_end_id(link_end_tuple[1])

        if link_end_tuple[0] == "nan" and len(link_end_tuple) == 2:
            return {
                observation.transmitter: spacecraft,
                observation.receiver: observation.body_reference_point_link_end_id(
                    "Earth", link_end_tuple[2]
                ),
            }
        else:
            return {
                observation.transmitter: observation.body_reference_point_link_end_id(
                    "Earth", link_end_tuple[0]
                ),
                observation.reflector1: spacecraft,
                observation.receiver: observation.body_reference_point_link_end_id(
                    "Earth", link_end_tuple[2]
                ),
            }

    def from_datetime_to_TBD(self, epoch, station):
        """
        Convert a datetime object in UTC into seconds since J2000 in TDB.

        Parameters
        ----------
        epoch : datetime
            The datetime object to convert.

        Returns
        -------
        float
            The time in seconds since J2000 in TDB.
        """
        if station not in self.stationDict:
            raise KeyError(
                "Error when processing TNF file, converting time from UTC to TDB: \n"
                + "the position of the ground station {} was not specified.".format(
                    station
                )
            )

        epoch_utc = time_conversion.datetime_to_tudat(epoch).epoch()
        epoch_tdb = self.time_scale_converter.convert_time(
            input_scale=time_conversion.utc_scale,
            output_scale=time_conversion.tdb_scale,
            input_value=epoch_utc,
            earth_fixed_position=self.stationDict[station],
        )

        return epoch_tdb
