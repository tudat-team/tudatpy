"""SLR CRD file loading."""

from tudatpy.kernel.data_input.tracking_data.slr import (
    CrdFullRateRecord,
    CrdMeteoRecord,
    CrdNormalPointRecord,
    CrdPass,
    CrdPassConfigurationData,
    CrdPassData,
    convert_crd_two_way_time_of_flight_to_slr_range,
    extract_full_rate_measurements,
    extract_full_rate_measurements_from_passes,
    extract_normal_point_measurements,
    extract_normal_point_measurements_from_passes,
    get_station_wavelengths,
    group_crd_data_per_station,
    group_crd_data_per_target,
    read_crd_file,
    read_crd_files,
)

__all__ = [
    "CrdFullRateRecord",
    "CrdMeteoRecord",
    "CrdNormalPointRecord",
    "CrdPass",
    "CrdPassConfigurationData",
    "CrdPassData",
    "convert_crd_two_way_time_of_flight_to_slr_range",
    "extract_full_rate_measurements",
    "extract_full_rate_measurements_from_passes",
    "extract_normal_point_measurements",
    "extract_normal_point_measurements_from_passes",
    "get_station_wavelengths",
    "group_crd_data_per_station",
    "group_crd_data_per_target",
    "read_crd_file",
    "read_crd_files",
]
