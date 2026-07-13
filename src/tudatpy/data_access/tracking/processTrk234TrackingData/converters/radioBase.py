"""
Base class for radiometric converters sharing helper functions.
"""

from tudatpy.astro import time_representation
from . import Converter
from trk234 import bands, SFDU
from datetime import datetime


class RadioBase(Converter):

    trkModeDict = {
        0: "Unknown",
        1: "1W",
        2: "2W",
        3: "3W",
    }

    # Maps trk234 band letters to the frequency band string identifiers recognized by
    # tudat::observation_models::getFrequencyBandFromString.
    frequencyBandsDict = {
        "S": "S-band",
        "X": "X-band",
        # "K": "Ku-band",
        "Ka": "Ka-band",
    }
    frequencyBandIds = {
        "S": 0.0,
        "X": 1.0,
        "Ka": 2.0,
    }

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

    def build_link_ends(
        self, link_end_tuple: tuple[str, str, str], spacecraftName: str | None = None
    ) -> list[tuple[tuple[str, str], str]]:
        """
        Construct a ``PlainLinkDefinition`` (list of ``((body, reference_point), link_end_role)``
        tuples) for Doppler/Range TrackingData creation.

        Parameters
        ----------
        link_end_tuple : tuple
            A tuple containing the uplink, spacecraft, and downlink identifiers.
        spacecraftName : str, optional
            The name of the spacecraft to use in the simulation. If not provided, the spacecraft name
            is extracted from the TNF file.

        Returns
        -------
        list[tuple[tuple[str, str], str]]
            The link ends constructed according to the following logic:
            - If the uplink identifier (first element) is "nan", use the spacecraft as the transmitter
            and assign the downlink using Earth's reference point (from the third element).
            - Otherwise, assign the transmitter using Earth's reference point from the first element,
            the reflector as the spacecraft, and the receiver using Earth's reference point from the third element.
        """

        if len(link_end_tuple) != 3:
            raise ValueError(
                "Error when processing TNF file, building link ends: \n"
                + f"the link end tuple should contain exactly 3 elements: {link_end_tuple} provided."
            )

        # Set custom spacecraft name if provided
        spacecraft_name = spacecraftName if spacecraftName is not None else link_end_tuple[1]
        spacecraft_link_end = (spacecraft_name, "")

        if link_end_tuple[0] == "nan":
            return [
                (spacecraft_link_end, "transmitter"),
                (("Earth", link_end_tuple[2]), "receiver"),
            ]
        else:
            return [
                (("Earth", link_end_tuple[0]), "transmitter"),
                (spacecraft_link_end, "reflector_1"),
                (("Earth", link_end_tuple[2]), "receiver"),
            ]

    def to_utc_epoch(self, datetime_utc: datetime) -> float:
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
        return time_representation.DateTime.from_python_datetime(datetime_utc).to_epoch()
