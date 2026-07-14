"""TNF/TRK-2-34 tracking-file loading."""

from ._converters.ramp import OpenRampHandling
from ._processor import TnfTrackingDataProcessor


def read_tnf_files(
    file_names: list[str],
    requested_observable_types: list[str],
    spacecraft_name: str | None = None,
    open_ramp_handling: OpenRampHandling = OpenRampHandling.print_warning_once,
):
    """Read TNF/TRK-2-34 files into tracking-data containers.

    Parameters
    ----------
    file_names : list[str]
        Paths to TNF/TRK-2-34 tracking files.
    requested_observable_types : list[str]
        Observable groups to extract from the files.
    spacecraft_name : str | None, default None
        Optional spacecraft body name used in generated link definitions.
    open_ramp_handling : OpenRampHandling, default OpenRampHandling.print_warning_once
        Policy for handling ramp intervals that do not have an explicit end record.

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
