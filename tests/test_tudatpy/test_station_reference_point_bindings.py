"""Python ownership regression; the numerical model is tested natively."""

from pathlib import Path

import numpy as np

from tudatpy.data_input.environment_data import station_history as sh
from tudatpy.dynamics.environment_setup import ground_station as gs
from tudatpy.dynamics import environment_setup
from tudatpy.astro.time_representation import Time


def test_ground_station_time_overload_keeps_sub_float_jump_order():
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
    assert left < right and float(left) == float(right)
    np.testing.assert_array_equal(station.get_cartesian_position(left), nominal)
    np.testing.assert_array_equal(station.get_cartesian_position(right), nominal + [0.1, 0.2, 0.3])
    np.testing.assert_array_equal(station.get_cartesian_state(left)[:3], nominal)
    np.testing.assert_array_equal(station.get_cartesian_state(right)[:3], nominal + [0.1, 0.2, 0.3])


def test_returned_offsets_and_trajectory_terms_are_independent_snapshots():
    identity = gs.StationReferencePointIdentity()
    identity.site = sh.SinexSiteKey("TEST", "A")
    identity.domes_id = "00000M001"
    identity.name = "marker"
    target = gs.StationReferencePointIdentity()
    target.site = identity.site
    target.domes_id = identity.domes_id
    target.kind = gs.StationReferencePointKind.antenna_reference_point
    target.name = "ARP"
    offset = gs.StationReferencePointOffset()
    offset.from_point, offset.to_point = identity, target
    offset.frame = "ITRF2020"
    offset.equipment_type, offset.equipment_serial_number = "antenna", "serial"
    offset.source = "explicit test reference vector"
    offset.components_metres = np.array([0.1, 0.2, 0.3])
    history = gs.StationReferencePointOffsetHistory([offset])
    query = gs.StationHistoryTime()
    expected = history.position_correction(query)
    snapshot = history.offsets[0]
    snapshot.components_metres = np.ones(3)
    snapshot.equipment_serial_number = "different hardware"
    snapshot.from_point.name = "different point"
    np.testing.assert_array_equal(history.position_correction(query), expected)
    assert history.offsets[0].equipment_serial_number == "serial"

    fixtures = Path(__file__).parents[1] / "test_tudat/reference_files/station_history"
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
    query = gs.ign_nominal_utc_station_time(2019, 200, 0)
    expected = trajectory.position(query)
    trajectory.solutions[3].position_metres = np.zeros(3)
    trajectory.post_seismic_terms[0].amplitude_metres = np.zeros(3)
    np.testing.assert_array_equal(trajectory.position(query), expected)
