from ._compat import deprecated_dir, deprecated_getattr, warn_custom_deprecation

_MIGRATION_GUIDE_URL = (
    "https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html"
)

_ALIASES = {
    "get_resource_path": "tudatpy.data_input.resource_paths.get_resource_path",
    "get_ephemeris_path": "tudatpy.data_input.resource_paths.get_ephemeris_path",
    "get_earth_orientation_path": "tudatpy.data_input.resource_paths.get_earth_orientation_path",
    "get_quadrature_path": "tudatpy.data_input.resource_paths.get_quadrature_path",
    "get_spice_kernel_path": "tudatpy.data_input.resource_paths.get_spice_kernel_path",
    "get_atmosphere_tables_path": "tudatpy.data_input.resource_paths.get_atmosphere_tables_path",
    "get_gravity_models_path": "tudatpy.data_input.resource_paths.get_gravity_models_path",
    "get_space_weather_path": "tudatpy.data_input.resource_paths.get_space_weather_path",
    "save2txt": "tudatpy.util.save2txt",
    "save_time_history_to_file": "tudatpy.util.save_time_history_to_file",
    "read_vector_history_from_file": "tudatpy.util.read_vector_history_from_file",
    "read_matrix_history_from_file": "tudatpy.util.read_matrix_history_from_file",
    "CrdPassConfigurationData": "tudatpy.data_input.tracking_data.slr.CrdPassConfigurationData",
    "CrdNormalPointRecord": "tudatpy.data_input.tracking_data.slr.CrdNormalPointRecord",
    "CrdFullRateRecord": "tudatpy.data_input.tracking_data.slr.CrdFullRateRecord",
    "CrdMeteoRecord": "tudatpy.data_input.tracking_data.slr.CrdMeteoRecord",
    "CrdPassData": "tudatpy.data_input.tracking_data.slr.CrdPassData",
    "CrdPass": "tudatpy.data_input.tracking_data.slr.CrdPass",
    "convert_crd_two_way_time_of_flight_to_slr_range": "tudatpy.data_input.tracking_data.slr.convert_crd_two_way_time_of_flight_to_slr_range",
    "read_crd_files": "tudatpy.data_input.tracking_data.slr.read_slr_data",
    "group_crd_data_per_target": "tudatpy.data_input.tracking_data.slr.group_crd_data_per_target",
    "group_crd_data_per_station": "tudatpy.data_input.tracking_data.slr.group_crd_data_per_station",
    "extract_normal_point_measurements": "tudatpy.data_input.tracking_data.slr.extract_normal_point_measurements",
    "extract_normal_point_measurements_from_passes": "tudatpy.data_input.tracking_data.slr.extract_normal_point_measurements_from_passes",
    "extract_full_rate_measurements": "tudatpy.data_input.tracking_data.slr.extract_full_rate_measurements",
    "extract_full_rate_measurements_from_passes": "tudatpy.data_input.tracking_data.slr.extract_full_rate_measurements_from_passes",
    "get_station_wavelengths": "tudatpy.data_input.tracking_data.slr.get_station_wavelengths",
    "SinexStationState": "tudatpy.data_input.environment_data.ilrs.SinexStationState",
    "SinexStationEccentricity": "tudatpy.data_input.environment_data.ilrs.SinexStationEccentricity",
    "IlrsStationRegistryEntry": "tudatpy.data_input.environment_data.ilrs.IlrsStationRegistryEntry",
    "convert_sinex_datetime_to_seconds_since_epoch": "tudatpy.data_input.environment_data.ilrs.convert_sinex_datetime_to_seconds_since_epoch",
    "read_sinex_station_data": "tudatpy.data_input.environment_data.ilrs.read_sinex_station_data",
    "read_sinex_station_eccentricities": "tudatpy.data_input.environment_data.ilrs.read_sinex_station_eccentricities",
    "read_ilrs_station_registry_from_sinex_site_id": "tudatpy.data_input.environment_data.ilrs.read_ilrs_station_registry_from_sinex_site_id",
    "read_domes_id_numbers": "tudatpy.data_input.environment_data.ilrs.read_domes_id_numbers",
    "read_monument_numbers": "tudatpy.data_input.environment_data.ilrs.read_monument_numbers",
    "read_ground_station_names": "tudatpy.data_input.environment_data.ilrs.read_ground_station_names",
    "TrackingDataType": "tudatpy.data_input.tracking_data.generic_text_file.TrackingDataType",
    "TrackingTxtFileReadFilterType": "tudatpy.data_input.tracking_data.generic_text_file.TrackingTxtFileReadFilterType",
    "TrackingTxtFileContents": "tudatpy.data_input.tracking_data.generic_text_file.TrackingTxtFileContents",
    "FdetDateFormat": "tudatpy.data_input.tracking_data.fdets.FdetDateFormat",
    "SolarActivityData": "tudatpy.data_input.environment_data.space_weather.SolarActivityData",
    "SolarActivityContainer": "tudatpy.data_input.environment_data.space_weather.SolarActivityContainer",
    "read_solar_activity_data": "tudatpy.data_input.environment_data.space_weather.read_solar_activity_data",
    "OdfDataType": "tudatpy.data_input.tracking_data.odf.OdfDataType",
    "OdfCommonDataBlock": "tudatpy.data_input.tracking_data.odf.OdfCommonDataBlock",
    "OdfDataSpecificBlock": "tudatpy.data_input.tracking_data.odf.OdfDataSpecificBlock",
    "OdfDopplerDataBlock": "tudatpy.data_input.tracking_data.odf.OdfDopplerDataBlock",
    "OdfDataBlock": "tudatpy.data_input.tracking_data.odf.OdfDataBlock",
    "OdfRampBlock": "tudatpy.data_input.tracking_data.odf.OdfRampBlock",
    "OdfRawFileContents": "tudatpy.data_input.tracking_data.odf.RawOdfFileContents",
    "read_odf_file": "tudatpy.data_input.tracking_data.odf.read_raw_odf_file_contents",
    "grail_antenna_file_reader": "tudatpy.data_input.environment_data.missions.grail.grail_antenna_file_reader",
    "grail_mass_level_0_file_reader": "tudatpy.data_input.environment_data.missions.grail.grail_mass_level_0_file_reader",
    "grail_mass_level_1_file_reader": "tudatpy.data_input.environment_data.missions.grail.grail_mass_level_1_file_reader",
    "LoadPDS": "tudatpy.data_input.data_retrieval.missions.LoadPDS",
    "DownloadAtmosphericData": "tudatpy.data_input.data_retrieval.media_corrections.DownloadAtmosphericData",
    "Trk234Processor": "tudatpy.data_input.tracking_data.tnf.TnfTrackingDataProcessor",
    "IonexProduct": "tudatpy.data_input.data_retrieval.media_corrections.IonexProduct",
    "IonexResolution": "tudatpy.data_input.data_retrieval.media_corrections.IonexResolution",
    "VmfTechnique": "tudatpy.data_input.data_retrieval.media_corrections.VmfTechnique",
    "VmfProcessing": "tudatpy.data_input.data_retrieval.media_corrections.VmfProcessing",
    "DownloadResult": "tudatpy.data_input.data_retrieval.media_corrections.DownloadResult",
    "AncillaryDownloadError": "tudatpy.data_input.data_retrieval.media_corrections.AncillaryDownloadError",
    "AuthenticationError": "tudatpy.data_input.data_retrieval.media_corrections.AuthenticationError",
    "download_ionex": "tudatpy.data_input.data_retrieval.media_corrections.download_ionex",
    "download_vmf": "tudatpy.data_input.data_retrieval.media_corrections.download_vmf",
    "download_ancillary": "tudatpy.data_input.data_retrieval.media_corrections.download_ancillary",
}

_CUSTOM_DEPRECATION_NAMES = {
    "read_crd_file",
    "read_tracking_txt_file",
    "read_ifms_file",
    "read_fdets_file",
    "set_dsn_weather_data_in_ground_stations",
    "set_estrack_weather_data_in_ground_stations",
}

__all__ = sorted(set(_ALIASES) | _CUSTOM_DEPRECATION_NAMES)


def _warn_non_one_to_one(name, replacement, details):
    warn_custom_deprecation(
        __name__,
        name,
        (
            f"Use {replacement} in the new data_input workflow. "
            f"{details} See the TudatPy migration guide: {_MIGRATION_GUIDE_URL}"
        ),
    )


def read_crd_file(file_name):
    _warn_non_one_to_one(
        "read_crd_file",
        "tudatpy.data_input.tracking_data.slr.read_slr_data",
        "The new function takes a list of CRD file names.",
    )
    from tudatpy.data_input.tracking_data.slr import read_slr_data

    return read_slr_data([file_name])


def read_tracking_txt_file(
    file_name,
    column_types,
    comment_symbol="#",
    value_separators=",:\t ",
    ignore_omitted_columns=False,
    data_filter_method=None,
):
    _warn_non_one_to_one(
        "read_tracking_txt_file",
        "tudatpy.data_input.tracking_data.generic_text_file.read_generic_text_data",
        "The new function takes a list of file names and returns a list of parsed files.",
    )
    from tudatpy.data_input.tracking_data.generic_text_file import (
        TrackingTxtFileReadFilterType,
        read_generic_text_data,
    )

    if data_filter_method is None:
        data_filter_method = TrackingTxtFileReadFilterType.no_tracking_txt_file_filter

    return read_generic_text_data(
        [file_name],
        column_types,
        comment_symbol,
        value_separators,
        ignore_omitted_columns,
        data_filter_method,
    )[0]


def read_ifms_file(*args, **kwargs):
    _warn_non_one_to_one(
        "read_ifms_file",
        "tudatpy.data_input.tracking_data.ifms.read_ifms_data",
        (
            "The new function has a different signature, requires spacecraft and "
            "ground-station metadata, and returns tracking data and supplementary "
            "data instead of raw file contents."
        ),
    )
    raise NotImplementedError(
        "tudatpy.data.read_ifms_file cannot be mapped one-to-one to the new IFMS "
        "reader. Use tudatpy.data_input.tracking_data.ifms.read_ifms_data."
    )


def read_fdets_file(*args, **kwargs):
    _warn_non_one_to_one(
        "read_fdets_file",
        "tudatpy.data_input.tracking_data.fdets.read_fdets_data",
        (
            "The new function has a different signature, requires base frequencies "
            "and station metadata, and returns tracking data and supplementary data "
            "instead of raw file contents."
        ),
    )
    raise NotImplementedError(
        "tudatpy.data.read_fdets_file cannot be mapped one-to-one to the new FDETS "
        "reader. Use tudatpy.data_input.tracking_data.fdets.read_fdets_data."
    )


def set_dsn_weather_data_in_ground_stations(*args, **kwargs):
    _warn_non_one_to_one(
        "set_dsn_weather_data_in_ground_stations",
        "tudatpy.dynamics.environment_setup.ground_station.set_dsn_weather_data_in_ground_station_settings",
        "Weather data files now have to be assigned in ground-station settings before the system of bodies is created.",
    )
    raise NotImplementedError(
        "DSN weather data must now be configured through ground-station settings. "
        "Use tudatpy.dynamics.environment_setup.ground_station."
        "set_dsn_weather_data_in_ground_station_settings."
    )


def set_estrack_weather_data_in_ground_stations(*args, **kwargs):
    _warn_non_one_to_one(
        "set_estrack_weather_data_in_ground_stations",
        "tudatpy.dynamics.environment_setup.ground_station.set_estrack_weather_data_in_ground_station_settings",
        "Weather data files now have to be assigned in ground-station settings before the system of bodies is created.",
    )
    raise NotImplementedError(
        "ESTRACK weather data must now be configured through ground-station settings. "
        "Use tudatpy.dynamics.environment_setup.ground_station."
        "set_estrack_weather_data_in_ground_station_settings."
    )


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), __all__)
