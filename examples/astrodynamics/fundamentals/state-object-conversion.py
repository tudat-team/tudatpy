"""
Example: State Object Conversions
---------------------------------
"""

###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro.two_body_dynamics import (
    CartesianState,
    KeplerianState)

spice_interface.load_standard_kernels()

# define reference frame
frame_origin = "SSB"
frame_orientation = "ECLIPJ2000"

# mission system parent body
parent_body = "Sun"
parent_body_mu = spice_interface.get_body_gravitational_parameter("Sun")

# mission planetary departure parameter definition
demo_date = "November 27, 2027"
demo_body = "Earth"

demo_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(demo_date)

demo_body_mu = spice_interface.get_body_gravitational_parameter(demo_body)
demo_body_state = spice_interface.get_body_cartesian_state_at_epoch(
    target_body_name=demo_body,
    observer_body_name=frame_origin,
    reference_frame_name=frame_orientation,
    aberration_corrections="none",
    ephemeris_time=demo_ephemeris_time)
demo_body_velocity = demo_body_state[3:]
demo_body_position = demo_body_state[:3]


def test_cartesian_state():
    test_cartesian_state = CartesianState(
        gravitational_parameter=demo_body_mu,
        position_vector=demo_body_position,
        velocity_vector=demo_body_velocity)

    test_keplerian_state = test_cartesian_state.to_keplerian()
    print(test_keplerian_state)

def test_keplerian_state():
    test_keplerian_state = KeplerianState(
        gravitational_parameter=demo_body_mu,
        semi_major_axis=1e5,
        eccentricity=0.0,
        inclination=0.0,
        longitude_ascending_node=0.0,
        argument_of_periapsis=0.0,
        true_anomaly=0.0)

    test_cartesian_state = test_keplerian_state.to_cartesian()
    print(test_cartesian_state)


if __name__ == "__main__":
    test_keplerian_state()
    test_cartesian_state()