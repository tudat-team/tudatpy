"""ATDF/TRK-2-25 tracking-file loading.

This module reads NASA Deep Space Network (DSN) closed-loop archival tracking
data files (ATDFs) and converts the supported Doppler, range, and ramp records
to Tudat tracking-data containers. ATDFs are the earliest closed-loop DSN
radio-science products, described in ``TRK-2-25``. The ``atdf2ascii`` Python
package is used to decode the binary ATDF files into intermediate ASCII
tables; this module reads those tables and converts them to Tudat
tracking-data and supplementary-data objects.
"""

from pathlib import Path

from ._processor import _DEFAULT_PROC_COUNT, AtdfTrackingDataProcessor
from tudatpy.data_input.tracking_data import TrackingData, TrackingSupplementaryData


def read_atdf_data(
    atdf_file_path: list[Path],
    spacecraft_name: str,
    output_dir: Path = Path("output"),
    count_time: list[float] | None = [60.0],
    proc_count: int = _DEFAULT_PROC_COUNT,
    doppler_one_way: bool = False,
    doppler_two_way: bool = True,
    doppler_three_way: bool = True,
    range_one_way: bool = False,
    range_two_way: bool = True,
) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:
    """Read ATDF/TRK-2-25 files into tracking-data containers.

    ATDF files are binary DSN closed-loop archival tracking data files
    containing radiometric spacecraft tracking data. This reader decodes the
    requested observable groups with ``atdf2ascii`` and returns Tudat tracking-
    data objects together with the corresponding supplementary data.
    Internally, this function creates an :class:`AtdfTrackingDataProcessor`,
    which uses the ``atdf2ascii`` Python package to decode the ATDF files into
    intermediate ``.msr``/``.ramp`` ASCII tables, dispatches the decoded
    records to the corresponding Tudat converters, and returns the generated
    Tudat ``TrackingData`` and ``TrackingSupplementaryData`` objects.

    Parameters
    ----------
    atdf_file_path : list[pathlib.Path]
        Paths to ATDF/TRK-2-25 tracking files.
    spacecraft_name : str
        Spacecraft body name used in generated link definitions (to be used
        for all files).
    output_dir : pathlib.Path, default Path("output")
        Directory in which ``atdf2ascii`` writes the intermediate ASCII
        tables.
    count_time : list[float] | None, default [60.0]
        Doppler count times (in seconds) to be decoded by ``atdf2ascii``.
    proc_count : int
        Number of processors to use for the ``atdf2ascii`` decoding step.
        Defaults to half the available cores.
    doppler_one_way, doppler_two_way, doppler_three_way, range_one_way, range_two_way : bool
        Observable groups to be decoded by ``atdf2ascii``. ``doppler_one_way``
        and ``range_one_way`` are reserved for future support and currently
        raise ``NotImplementedError`` if set to ``True``, since no converter
        exists yet for 1-way Doppler/range data.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    processor = AtdfTrackingDataProcessor(
        atdf_file_path=atdf_file_path,
        spacecraft_name=spacecraft_name,
        proc_count=proc_count,
        doppler_one_way=doppler_one_way,
        doppler_two_way=doppler_two_way,
        doppler_three_way=doppler_three_way,
        range_one_way=range_one_way,
        range_two_way=range_two_way,
    )
    processor._convert_atdf_to_ascii(output_dir, count_time=count_time)
    return processor.process(output_dir)
