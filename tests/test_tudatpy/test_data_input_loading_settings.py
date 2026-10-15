"""Cross-module binding checks; loading mathematics and parsers are tested in C++."""

from pathlib import Path

import numpy as np

from tudatpy.astro import frame_conversion
from tudatpy.data_input.environment_data import (
    atmospheric_loading,
    ocean_loading,
    station_loading,
    terrestrial_frames,
)
from tudatpy.dynamics.environment_setup import rotation_model, shape_deformation

TEST_DATA = Path(__file__).parents[1] / "test_tudat"


def test_station_loading_records_and_preloaded_settings_registration():
    metadata = station_loading.imls_merra2_atmospheric_loading_metadata("binding release")
    metadata.reference_frame = "ITRS"
    data = station_loading.read_ephedisp_station_loading_data(
        str(TEST_DATA / "src/io/station_loading_fixtures/imls_merra2_extract.eph"),
        metadata,
    )
    assert isinstance(data, station_loading.StationLoadingData)
    assert isinstance(data.epochs[0], station_loading.StationLoadingEpoch)
    site = data.sites["weuc"]
    assert isinstance(site.cartesian_position, list)
    assert isinstance(site.samples[0].displacement, list)
    assert len(site.samples[0].displacement) == 3
    assert data.metadata.release == "binding release"

    policy = shape_deformation.StationLoadingEvaluationSettings()
    policy.station_bindings = {"EnvironmentStation": "weuc"}
    policy.output_origin_definition = data.metadata.origin_definition
    assert policy.origin_translation is None
    translation = shape_deformation.StationLoadingOriginTranslation()
    translation.series = data
    policy.origin_translation = translation
    assert isinstance(policy.origin_translation.series, station_loading.StationLoadingData)
    policy.origin_translation = None
    settings = shape_deformation.station_loading(data, policy)
    assert isinstance(settings, shape_deformation.StationLoadingBodyDeformationSettings)
    assert isinstance(settings, shape_deformation.BodyDeformationSettings)


def test_tidal_loading_records_cross_into_environment_settings():
    atmosphere = atmospheric_loading.read_atmospheric_tide_station_data(
        str(TEST_DATA / "reference_files/atmospheric_tides/bernese_manual_stations.atl")
    )
    assert isinstance(atmosphere, atmospheric_loading.AtmosphericTideStationData)
    assert next(iter(atmosphere.stations.values())).coefficients_metres.shape == (3, 4)
    policy = shape_deformation.AtmosphericTideEvaluationSettings()
    policy.displacement_basis = (
        atmospheric_loading.AtmosphericTideDisplacementBasis.geocentric_radial_north_east
    )
    settings = shape_deformation.atmospheric_tide_from_stations(atmosphere, policy)
    assert isinstance(settings, shape_deformation.AtmosphericTideBodyDeformationSettings)

    ocean = ocean_loading.read_ocean_pole_tide_grid(
        str(TEST_DATA / "reference_files/ocean_pole_tide/iers_grid_extract.txt")
    )
    assert isinstance(ocean, ocean_loading.OceanPoleTideGrid)
    assert ocean.coefficients[0].shape == (3, 2)
    assert isinstance(
        shape_deformation.ocean_pole_tidal_from_data(ocean),
        shape_deformation.OceanPoleTideBodyDeformationSettings,
    )


# Raw frame records must cross into numerical frame operations and rotation settings.
def test_earth_frame_and_rotation_settings_registration():
    frames = terrestrial_frames.official_itrf2020_parameters()
    assert isinstance(frames, terrestrial_frames.TerrestrialFrameTransformationData)
    parameters = frames.get_parameters("ITRF2014")
    assert isinstance(parameters, terrestrial_frames.TerrestrialFrameTransformationParameters)
    assert parameters.source_frame.realization == "ITRF2020"
    assert parameters.translation_metres.shape == (3,)
    assert len(parameters.rotation_rate_radians_per_second) == 3

    # Exercise the Helmert bindings using an explicit decimal year and state-time rate.
    model = frame_conversion.HelmertTransformationModel(frames)
    state = np.array([1.0, 2.0, 3.0, 0.1, 0.2, 0.3])
    year, year_rate = 2025.0, 1.0 / 31557600.0
    affine = model.get_transformation("ITRF2020", "ITRF2014", year, year_rate)
    direct = frame_conversion.evaluate_linearized_helmert(parameters, year, year_rate)
    transformed = model.transform_state(state, "ITRF2020", "ITRF2014", year, year_rate)
    np.testing.assert_array_equal(transformed, direct.transform_state(state))
    np.testing.assert_array_equal(
        model.transform_position(state[:3], "ITRF2020", "ITRF2014", year),
        transformed[:3],
    )
    np.testing.assert_array_equal(affine.matrix, direct.matrix)
    assert affine.inverse().compose_after(affine).transform_position(state[:3]).shape == (3,)

    conventional = rotation_model.gcrs_to_itrs()
    assert conventional.ocean_tide_eop_model == rotation_model.OceanTideEopModel.iers_2010
    opted_in = rotation_model.gcrs_to_itrs(
        ocean_tide_eop_model=rotation_model.OceanTideEopModel.desai_sibois_2016
    )
    assert opted_in.ocean_tide_eop_model == rotation_model.OceanTideEopModel.desai_sibois_2016
