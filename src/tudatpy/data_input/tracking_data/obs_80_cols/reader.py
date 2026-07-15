"""MPC 80-column optical astrometry loading.

The MPC 80-column format is a fixed-width optical-observation format for minor
planets, comets, and natural satellites:
https://minorplanetcenter.net/iau/info/OpticalObs.html.
"""

from .parsers import parse_80cols_file
from tudatpy.data_input.tracking_data.optical_utilities import read_astropy_optical_data


def read_80_column_data(
    file_names: list[str],
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
):
    """Read MPC 80-column files into TrackingData objects.

    The input files are parsed as MPC fixed-width optical astrometry records.
    The :func:`parse_80cols_file` function reads the fixed-width records into
    an astropy table with standardized optical astrometry columns. This table
    is then passed to
    :func:`~tudatpy.data_input.tracking_data.optical_utilities.read_astropy_optical_data`,
    which creates the augmented optical table with
    :func:`~tudatpy.data_input.tracking_data.optical_utilities.create_augmented_optical_table`
    and converts the observations to Tudat optical tracking data through the
    common optical-data conversion pipeline.

    Parameters
    ----------
    file_names : list[str]
        Paths to MPC 80-column optical astrometry files.
    frame : str, default "J2000"
        Reference frame of the input observations, e.g. the frame in which the
        right ascension and declination are defined.
    custom_name : str | None, default None
        Optional target name assigned to all observations, to be used if the
        file-provided object name is different from the Tudat body name.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata. To be used only
        in special cases where additional metadata is to be associated with the
        observations.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    return read_astropy_optical_data(
        parse_80cols_file(file_names),
        in_degrees=False,
        frame=frame,
        custom_name=custom_name,
        add_weights=add_weights,
        add_star_catalog_corrections=add_star_catalog_corrections,
        add_ancillary_data=add_ancillary_data,
    )
