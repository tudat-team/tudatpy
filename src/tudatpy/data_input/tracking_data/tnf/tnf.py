"""TNF/TRK-2-34 tracking-file loading.

This module reads NASA Deep Space Network (DSN) Tracking and Navigation Files
and converts supported Doppler, range, and ramp records to Tudat tracking-data
containers. The binary TNF format is described in ``820-013, TRK-2-34 DSN
Tracking System Data Archival Format``.
The bulk of the TNF/TRK-2-34 parsing is handled by the ``pytrk234`` Python
package; this module converts the decoded records to Tudat tracking-data and
supplementary-data objects.
"""

from tudatpy.data_input.tracking_data import TrackingData, TrackingSupplementaryData

from ._converters.ramp import OpenRampHandling
from ._processor import TnfTrackingDataProcessor


def read_tnf_data(
    file_names: list[str],
    requested_observable_types: list[str],
    spacecraft_name: str | None = None,
    open_ramp_handling: OpenRampHandling = OpenRampHandling.print_warning_once,
) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:
    """Read TNF/TRK-2-34 files into tracking-data containers.

    TNF/TRK-2-34 files are binary DSN Tracking and Navigation Files containing
    radiometric spacecraft tracking data. This reader extracts the requested
    observable groups and returns Tudat tracking-data objects together with the
    corresponding supplementary data.
    Internally, this function creates a :class:`TnfTrackingDataProcessor`, which
    uses the ``pytrk234`` Python package to decode the TNF/TRK-2-34 files,
    dispatches the decoded records to the corresponding Tudat converters,
    applies the selected :class:`OpenRampHandling` policy to ramp records, and
    returns the generated Tudat ``TrackingData`` and
    ``TrackingSupplementaryData`` objects.

    Parameters
    ----------
    file_names : list[str]
        Paths to TNF/TRK-2-34 tracking files.
    requested_observable_types : list[str]
        Observable groups to extract from the files (to be used for all files).
    spacecraft_name : str | None, default None
        Optional spacecraft body name used in generated link definitions (to be
        used for all files).
    open_ramp_handling : OpenRampHandling, default OpenRampHandling.print_warning_once
        Policy for handling ramp intervals that do not have an explicit end
        record (to be used for all files).

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    processor = TnfTrackingDataProcessor(
        file_names,
        requested_observable_types,
        spacecraft_name=spacecraft_name,
    )
    return processor.process(open_ramp_handling)
