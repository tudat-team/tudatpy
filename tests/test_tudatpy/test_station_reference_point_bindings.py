"""Python ownership regression; the numerical model is tested natively."""

from pathlib import Path

import numpy as np

from tudatpy.data_input import resource_paths
from tudatpy.data_input.environment_data import station_history as sh
from tudatpy.data_input.environment_data import terrestrial_frames
from tudatpy.dynamics.environment_setup import ground_station as gs
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics import environment
from tudatpy.astro.time_representation import Time


def test_ground_station_time_overload_keeps_sub_float_jump_order():
    """Dispatch distinct precise epochs through the public station position and state overloads."""
    settings = environment_setup.BodyListSettings("SSB", "J2000")
    settings.add_empty_settings("Earth")
    settings.get("Earth").shape_settings = environment_setup.shape.spherical(6378137.0)
    bodies = environment_setup.create_system_of_bodies(settings)
    boundary = 1.0e9
    motion = gs.piecewise_constant_station_motion(
        {boundary - 1.0: np.zeros(3), boundary: np.array([0.1, 0.2, 0.3])}
    )
    nominal = np.array([6378137.0, 0.0, 0.0])
    environment_setup.add_ground_station(
        bodies.get("Earth"),
        gs.basic_station("TimeDispatch", nominal, station_motion_settings=[motion]),
    )
    station = bodies.get("Earth").get_ground_station("TimeDispatch").station_state
    left, right = Time(boundary) - 1.0e-10, Time(boundary) + 1.0e-10

    # Equal float values must still select opposite sides of the precise station-motion jump.
    assert left < right and float(left) == float(right)
    np.testing.assert_array_equal(station.get_cartesian_position(left), nominal)
    np.testing.assert_array_equal(station.get_cartesian_position(right), nominal + [0.1, 0.2, 0.3])
    np.testing.assert_array_equal(station.get_cartesian_state(left)[:3], nominal)
    np.testing.assert_array_equal(station.get_cartesian_state(right)[:3], nominal + [0.1, 0.2, 0.3])


# Mutating records returned by runtime models must leave their validated coefficients unchanged.
def test_returned_offsets_and_trajectory_terms_are_independent_snapshots():
    """Keep source and transformed coefficient records independent from returned Python snapshots."""
    fixtures = (
        Path(resource_paths.get_test_data_path()) / "pr1092" / "reference_files/station_history"
    )
    products = [
        sh.read_sinex_station_history(str(fixtures / filename))
        for filename in ("ab07_solutions.snx", "ab07_discontinuities.snx", "ab07_psd.snx")
    ]
    for product in products:
        product.frame = product.release = "ITRF2020"
        product.origin = "CM"
        product.tide_system = "explicit test tide convention"
    trajectory = gs.ign_station_history(
        products[0],
        products[1],
        sh.SinexSiteKey("AB07", "A"),
        post_seismic=products[2],
        enu_to_terrestrial=np.eye(3),
    )
    query = environment.ign_nominal_utc_station_time(2019, 200, 0)

    # Editing exported solution and PSD records must not alter the validated trajectory.
    expected = trajectory.position(query)
    trajectory.solutions[3].position_metres = np.zeros(3)
    trajectory.post_seismic_terms[0].amplitude_metres = np.zeros(3)
    np.testing.assert_array_equal(trajectory.position(query), expected)

    # The frame transformation returns a new history; offsets are added in the frame of that history.
    target_frame = terrestrial_frames.TerrestrialFrameIdentity()
    target_frame.realization = target_frame.release = "ITRF2014"
    transformed = gs.transform_station_history(
        trajectory,
        terrestrial_frames.official_itrf2020_parameters(),
        target_frame,
        trajectory.metadata.origin,
        trajectory.metadata.tide_system,
    )
    assert transformed.metadata.frame == "ITRF2014"
    assert trajectory.metadata.frame == "ITRF2020"

    # The station factory identifies the ARP in the frame where its offset is supplied.
    point_offset = environment.StationReferencePointOffset()
    point_offset.from_point = environment.station_reference_point_identity(transformed.metadata)
    point_offset.to_point = point_offset.from_point
    point_offset.to_point.kind = environment.StationReferencePointKind.antenna_reference_point
    point_offset.to_point.name = "AB07:ARP"
    point_offset.frame = transformed.metadata.frame
    point_offset.components_metres = np.array([0.1, 0.2, 0.3])
    point_offset.equipment_type = "antenna"
    point_offset.equipment_serial_number = "serial"
    point_offset.source = "explicit test reference vector"
    point_history = environment.StationReferencePointOffsetHistory([point_offset])
    station = gs.reference_point_station(transformed, [point_history], np.zeros(3))
    assert station.station_name == "AB07:ARP"

    # Returned coefficient and offset records remain copies.
    query = environment.ign_nominal_utc_station_time_from_tai(query.tai, 1.0)
    expected_point = transformed.position(query) + point_offset.components_metres
    np.testing.assert_array_equal(
        transformed.position(query) + point_history.position_correction(query), expected_point
    )
    transformed.solutions[3].position_metres = np.zeros(3)
    transformed.post_seismic_terms[0].amplitude_metres = np.zeros(3)
    point_history.offsets[0].components_metres = np.zeros(3)
    np.testing.assert_array_equal(
        transformed.position(query) + point_history.position_correction(query), expected_point
    )
