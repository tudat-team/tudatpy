"""
MPC ADES tracking-data loading.
"""

from .parsers import parse_ades_file
from tudatpy.data_input.tracking_data.optical_utilities import read_astropy_optical_data
from tudatpy.data_input.tracking_data.radar_utilities import (
    radar_data_from_table,
    radar_data_to_tracking_data,
)


def read_ades_data(
    file_path: str,
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    *,
    ades_weights: bool = False,
    vfcc17_weights: bool = False,
    nightly_deweighting: bool = False,
):
    """Function to read an MPC ADES file into TrackingData objects.

    Function to read the observations in an MPC ADES file (in XML or PSV format), and convert
    them to tracking data objects that can be used to create an observation collection. The file
    is first parsed into an astropy table using :func:`parse_ades_file`. Stellar occultation
    observations are converted to the right ascension and declination of the target, and are
    then processed in the same manner as optical observations, using
    :func:`~tudatpy.data_input.tracking_data.optical_utilities.read_astropy_optical_data`.
    Radar observations are converted using
    :func:`~tudatpy.data_input.tracking_data.radar_utilities.radar_data_to_tracking_data`.

    Parameters
    ----------
    file_path : str
        Path to the MPC ADES file (in XML or PSV format).
    frame : str, default "J2000"
        Reference frame of the input observations, e.g. the frame in which the
        right ascension and declination are defined.
    custom_name : str | None, default None
        Optional target name assigned to all observations, to be used if the
        file-provided object name is different from the Tudat body name.
    add_weights : bool, default False
        Deprecated, use ``vfcc17_weights`` instead.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.
    ades_weights : bool, default False
        Whether the rmsRA and rmsDec uncertainties of each observation (in arcseconds) are used as
        weights (see
        :meth:`~tudatpy.data_input.tracking_data.TrackingDataWeightSettings.set_observation_weight_settings`).
        Requires rmsRA and rmsDec for every observation.
    vfcc17_weights : bool, default False
        Whether the observations get the weights of :cite:p:`veres2017` (see
        :meth:`~tudatpy.data_input.tracking_data.TrackingDataWeightSettings.set_named_weights`).
        Cannot be combined with ``ades_weights``.
    nightly_deweighting : bool, default False
        Whether the ADES weights are divided by the number of observations of the same target,
        from the same observatory, in the same local night, as defined for the weights of
        :cite:p:`veres2017` (see
        :meth:`~tudatpy.data_input.tracking_data.TrackingDataWeightSettings.set_block_deweighting`).
        Requires ``ades_weights``.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.

    Notes
    -----
    Offset observations, radar observations that are not w.r.t. the center of mass of the
    target, and observations from Earth-based roving observers (with positions given in WGS84,
    ITRF or IAU) are not yet supported, and are skipped (with a warning).
    """
    parsed_table = parse_ades_file(file_path)
    optical_tracking_data, supplementary_data = [], []
    if len(parsed_table) > 0:
        optical_tracking_data, supplementary_data = read_astropy_optical_data(
            parsed_table,
            in_degrees=False,
            frame=frame,
            custom_name=custom_name,
            add_weights=add_weights,
            add_star_catalog_corrections=add_star_catalog_corrections,
            add_ancillary_data=add_ancillary_data,
            ades_weights=ades_weights,
            vfcc17_weights=vfcc17_weights,
            nightly_deweighting=nightly_deweighting,
        )
    radar_data = radar_data_from_table(parsed_table)
    if custom_name is not None:
        radar_data = radar_data.assign(target_body=str(custom_name))
    radar_tracking_data, radar_supplementary_data = radar_data_to_tracking_data(radar_data)
    return (
        optical_tracking_data + radar_tracking_data,
        supplementary_data + radar_supplementary_data,
    )
