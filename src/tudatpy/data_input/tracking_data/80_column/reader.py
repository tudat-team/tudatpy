"""MPC 80-column optical astrometry loading."""

from .parsers import parse_80cols_file
from tudatpy.data_input.tracking_data.optical_utilities import read_astropy_optical_data


def read_80_column_data(
    filename: str | list[str],
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str = "",
):
    """Read MPC 80-column files into TrackingData objects.

    Parameters
    ----------
    filename : str | list[str]
        Path or paths to MPC 80-column optical astrometry files.
    frame : str, default "J2000"
        Reference frame of the input observations.
    custom_name : str | None, default None
        Optional target name assigned to all observations.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.
    weighing_scheme : str, default ""
        Name of the weighing scheme assigned to the generated tracking data.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    return read_astropy_optical_data(
        parse_80cols_file(filename),
        in_degrees=False,
        frame=frame,
        custom_name=custom_name,
        add_weights=add_weights,
        add_star_catalog_corrections=add_star_catalog_corrections,
        add_ancillary_data=add_ancillary_data,
        weighing_scheme=weighing_scheme,
    )
