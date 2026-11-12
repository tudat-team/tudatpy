"""
MPC ADES tracking-data loading.
"""

import numpy as np

from .parsers import parse_ades_file
from tudatpy.data_input.tracking_data.optical_utilities import read_astropy_optical_data
from tudatpy.data_input.tracking_data.radar_utilities import (
    radar_data_from_table,
    radar_data_to_tracking_data,
)

# ADES fields with the uncertainties of an observation as reported by the observer (observation fields), and the corresponding
# uncertainties adopted by an orbit computing center (residual fields), which also provide the adopted biases
OBSERVATION_UNCERTAINTY_FIELDS = ("rmsRA", "rmsDec", "rmsCorr", "rmsTime")
RESIDUAL_UNCERTAINTY_FIELDS = ("sigRA", "sigDec", "sigCorr", "sigTime")
RESIDUAL_BIAS_FIELDS = ("biasRA", "biasDec", "biasTime")


def _select_ades_observation_or_residual_fields(table, prefer_residual_fields: bool):
    """Select, for each observation, either the ADES observation fields or the residual fields.

    For each observation, all fields of the preferred kind are used if any of them is available,
    and otherwise all fields of the other kind; the two kinds are never combined for a single
    observation. The observation fields are rmsRA, rmsDec, rmsCorr and rmsTime; the residual
    fields are sigRA, sigDec, sigCorr and sigTime (uncertainties adopted by an orbit computing
    center) and biasRA, biasDec and biasTime (biases adopted by an orbit computing center).
    The selected uncertainties are stored in the rmsRA, rmsDec, rmsCorr and rmsTime columns
    (the sigRA, sigDec, sigCorr and sigTime columns are removed), and the biasRA, biasDec and
    biasTime columns are kept only for the observations that use the residual fields (NaN
    otherwise).

    Parameters
    ----------
    table : astropy.table.Table
        Table with the parsed ADES observations (see :func:`parse_ades_file`), modified in place.
    prefer_residual_fields : bool
        Whether the residual fields (instead of the observation fields) are used for an
        observation for which both are available.
    """
    values = {
        column: (
            np.ma.filled(np.ma.asarray(table[column], dtype=float), np.nan)
            if column in table.colnames
            else np.full(len(table), np.nan)
        )
        for column in OBSERVATION_UNCERTAINTY_FIELDS
        + RESIDUAL_UNCERTAINTY_FIELDS
        + RESIDUAL_BIAS_FIELDS
    }
    has_observation_fields = np.any(
        [~np.isnan(values[column]) for column in OBSERVATION_UNCERTAINTY_FIELDS], axis=0
    )
    has_residual_fields = np.any(
        [
            ~np.isnan(values[column])
            for column in RESIDUAL_UNCERTAINTY_FIELDS + RESIDUAL_BIAS_FIELDS
        ],
        axis=0,
    )
    use_residual_fields = has_residual_fields & (prefer_residual_fields | ~has_observation_fields)

    for observation_field, residual_field in zip(
        OBSERVATION_UNCERTAINTY_FIELDS, RESIDUAL_UNCERTAINTY_FIELDS
    ):
        table[observation_field] = np.where(
            use_residual_fields, values[residual_field], values[observation_field]
        )
    for column in RESIDUAL_BIAS_FIELDS:
        table[column] = np.where(use_residual_fields, values[column], np.nan)
    table.remove_columns(
        [column for column in RESIDUAL_UNCERTAINTY_FIELDS if column in table.colnames]
    )


def read_ades_data(
    file_path: str,
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str | None = "",
    minimum_ades_uncertainty: float | None = None,
    prefer_residual_fields: bool = False,
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
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.
    weighing_scheme : str, default ""
        Weighing scheme that is to be used for the observations: ``"VFCC17"`` (see
        :cite:p:`veres2017`), ``"ADES"`` (weights computed from the rmsRA and rmsDec
        uncertainties of each observation, and from their correlation rmsCorr if
        available) or ``"ADES_NIGHT"`` (as ``"ADES"``, but with the
        weights divided by the number of observations of the same target,
        from the same observatory, in the same local night, as defined for ``"VFCC17"``). In
        both ADES schemes, observations without rmsRA and rmsDec get ``"VFCC17"`` weights. The
        weights are computed when the observation collection is created. If empty, no weights
        are assigned. If weights are assigned, any rmsTime timing uncertainties are added to
        the observation covariance (see
        :meth:`~tudatpy.data_input.tracking_data.TrackingData.set_observation_time_uncertainties`).
    minimum_ades_uncertainty : float | None, default None
        Minimum rmsRA and rmsDec (in arcseconds) for which the ADES uncertainties of an
        observation are used (applied to sigRA and sigDec for observations that use the
        residual fields, see ``prefer_residual_fields``). For observations with a smaller value,
        these are removed, so that the ADES weighing schemes assign ``"VFCC17"`` weights instead.
    prefer_residual_fields : bool, default False
        Whether the uncertainties and biases adopted by an orbit computing center (the ADES
        residual fields sigRA, sigDec, sigCorr, sigTime, biasRA, biasDec and biasTime) are used
        instead of the uncertainties reported by the observer (the ADES observation fields
        rmsRA, rmsDec, rmsCorr and rmsTime), for observations that have both. Each observation
        uses either all fields of one kind or all fields of the other: those of the preferred
        kind if any of them is available, otherwise those of the other kind. The residual uncertainties then
        take the place of rmsRA, rmsDec, rmsCorr and rmsTime in the weighing schemes. The
        adopted biases (only available in the residual fields) are applied as follows, assuming
        that the reported values are the true values plus the bias: biasTime (in seconds) is
        set as the negative observation time bias (see
        :meth:`~tudatpy.data_input.tracking_data.TrackingData.set_observation_time_biases`), and
        biasRA (of RA cos(Dec)) and biasDec (in arcseconds) are subtracted as observation
        corrections, replacing the star-catalog corrections of these observations.

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
        _select_ades_observation_or_residual_fields(parsed_table, prefer_residual_fields)
        if minimum_ades_uncertainty is not None:
            rejected = (parsed_table["rmsRA"] < minimum_ades_uncertainty) | (
                parsed_table["rmsDec"] < minimum_ades_uncertainty
            )
            for column in ("rmsRA", "rmsDec"):
                parsed_table[column][rejected] = np.nan
        optical_tracking_data, supplementary_data = read_astropy_optical_data(
            parsed_table,
            in_degrees=False,
            frame=frame,
            custom_name=custom_name,
            add_weights=add_weights,
            add_star_catalog_corrections=add_star_catalog_corrections,
            add_ancillary_data=add_ancillary_data,
            weighing_scheme=weighing_scheme,
        )
    radar_data = radar_data_from_table(parsed_table)
    if custom_name is not None:
        radar_data = radar_data.assign(target_body=str(custom_name))
    radar_tracking_data, radar_supplementary_data = radar_data_to_tracking_data(radar_data)
    return (
        optical_tracking_data + radar_tracking_data,
        supplementary_data + radar_supplementary_data,
    )
