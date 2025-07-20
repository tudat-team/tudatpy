import typing
from . import Converter as Converter
from _typeshed import Incomplete
from tudatpy.astro import time_representation as time_representation
from tudatpy.dynamics.environment_setup.ground_station import get_approximate_dsn_ground_station_positions as get_approximate_dsn_ground_station_positions
from tudatpy.estimation.observable_models_setup import links as links
from tudatpy.estimation.observations_setup.ancillary_settings import FrequencyBands as FrequencyBands

class RadioBase(Converter):
    trkModeDict: Incomplete
    time_scale_converter: Incomplete
    frequencyBandsDict: Incomplete
    stationDict: Incomplete

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

    def build_link_ends_dict(self, link_end_tuple, spacecraftName: Incomplete | None=None):
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
            and assign the downlink using Earth\'s reference point (from the third element).
            - Otherwise, assign the transmitter using Earth\'s reference point from the first element,
            the reflector as the spacecraft, and the receiver using Earth\'s reference point from the third element.
        """

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