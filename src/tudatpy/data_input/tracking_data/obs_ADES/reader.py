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
    weighing_scheme: str | None = "",
):
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
            weighing_scheme=weighing_scheme,
        )
    return (optical_tracking_data, supplementary_data)
