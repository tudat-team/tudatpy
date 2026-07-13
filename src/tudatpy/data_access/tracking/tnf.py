"""TNF/TRK-2-34 tracking-file loading."""

from tudatpy.data_access.tracking.processTrk234TrackingData import (
    OpenRampHandling,
    Trk234TrackingDataProcessor,
)


def read_tnf_files(
    file_names: list[str],
    requested_observable_types: list[str],
    spacecraft_name: str | None = None,
    open_ramp_handling: OpenRampHandling = OpenRampHandling.print_warning_once,
):
    processor = Trk234TrackingDataProcessor(
        file_names,
        requested_observable_types,
        spacecraft_name=spacecraft_name,
    )
    return processor.process(open_ramp_handling)


__all__ = [
    "OpenRampHandling",
    "Trk234TrackingDataProcessor",
    "read_tnf_files",
]
