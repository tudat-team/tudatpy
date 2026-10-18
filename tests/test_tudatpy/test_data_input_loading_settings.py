"""Exercise loaded environmental records through their consuming Python factories."""

from pathlib import Path

import numpy as np

from tudatpy.astro import frame_conversion
from tudatpy.data_input import resource_paths
from tudatpy.data_input.environment_data import (
    atmospheric_loading,
    ocean_loading,
    station_loading,
    terrestrial_frames,
)
from tudatpy.dynamics.environment_setup import rotation_model, shape_deformation

TEST_DATA = Path(resource_paths.get_test_data_path()) / "pr1092"


def test_loaded_products_create_deformation_settings():
    """Pass real station, atmospheric and ocean-pole records into their deformation factories."""
    loading = station_loading.read_ephedisp_station_loading_data(
        str(TEST_DATA / "src/io/station_loading_fixtures/imls_merra2_extract.eph")
    )

    # Provider REN values stay in metres when the raw series enters the loading settings.
    np.testing.assert_array_equal(
        loading.sites["weuc"].displacements[0], [-0.00510, -0.00202, 0.00170]
    )
    assert isinstance(
        shape_deformation.station_loading(loading, {"EnvironmentStation": "weuc"}),
        shape_deformation.StationLoadingBodyDeformationSettings,
    )

    atmosphere = atmospheric_loading.read_atmospheric_tide_station_data(
        str(TEST_DATA / "reference_files/atmospheric_tides/bernese_manual_stations.atl")
    )
    policy = shape_deformation.AtmosphericTideEvaluationSettings()
    policy.displacement_basis = (
        atmospheric_loading.AtmosphericTideDisplacementBasis.geocentric_radial_north_east
    )

    # The signed provider station coefficient is converted to metres before construction.
    assert atmosphere.stations["ABPO"].coefficients_metres[0, 0] == -0.3414 * 1.0e-3
    assert isinstance(
        shape_deformation.atmospheric_tide_from_stations(atmosphere, policy),
        shape_deformation.AtmosphericTideBodyDeformationSettings,
    )

    ocean = ocean_loading.read_ocean_pole_tide_grid(
        str(TEST_DATA / "reference_files/ocean_pole_tide/iers_grid_extract.txt")
    )

    # Ocean-pole dimensionless response coefficients pass directly into the consuming factory.
    assert ocean.coefficients[0][0, 0] == 0.216192
    assert isinstance(
        shape_deformation.ocean_pole_tidal_from_data(ocean),
        shape_deformation.OceanPoleTideBodyDeformationSettings,
    )


def test_earth_frame_and_rotation_settings_registration():
    """Evaluate the public Helmert position/state operations and select both ocean-tide EOP models."""
    frames = terrestrial_frames.official_itrf2020_parameters()
    parameters = frames.get_parameters("ITRF2014")
    model = frame_conversion.HelmertTransformationModel(frames)
    state = np.array([1.0, 2.0, 3.0, 0.1, 0.2, 0.3])
    year, year_rate = 2025.0, 1.0 / 31557600.0
    transformation = model.get_transformation("ITRF2020", "ITRF2014", year, year_rate)
    direct = frame_conversion.evaluate_linearized_helmert(parameters, year, year_rate)
    transformed = model.transform_state(state, "ITRF2020", "ITRF2014", year, year_rate)

    # Catalogue and direct-parameter operations return the same position and full state.
    np.testing.assert_array_equal(transformed, direct.transform_state(state))
    np.testing.assert_array_equal(
        model.transform_position(state[:3], "ITRF2020", "ITRF2014", year), transformed[:3]
    )
    np.testing.assert_array_equal(transformation.matrix, direct.matrix)
    assert transformation.inverse().compose_after(transformation).transform_position(
        state[:3]
    ).shape == (3,)

    conventional = rotation_model.gcrs_to_itrs()
    opted_in = rotation_model.gcrs_to_itrs(
        ocean_tide_eop_model=rotation_model.OceanTideEopModel.desai_sibois_2016
    )

    # The default and explicitly selected short-period models retain distinct settings.
    assert conventional.ocean_tide_eop_model == rotation_model.OceanTideEopModel.iers_2010
    assert opted_in.ocean_tide_eop_model == rotation_model.OceanTideEopModel.desai_sibois_2016
