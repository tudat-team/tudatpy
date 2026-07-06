"""
Functions to calculate photocenter corrections to observations
"""
import numpy as np
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.dynamics.environment import SystemOfBodies
from numpy.linalg import norm
from tudatpy.estimation.observable_models_setup.biases.observation_correction_utils import _offset_vector_to_corrections, _unit


def _photocenter_offset(diameter : float,
                        body_wrt_observer : np.ndarray,
                        body_wrt_sun : np.ndarray,) -> np.ndarray:
    """Helper function to calculate photocenter offset at one epoch"""
    # Get unit vectors
    body_wrt_observer_unit = _unit(body_wrt_observer)
    body_wrt_sun_unit = _unit(body_wrt_sun)

    # Angle spanned by observer, body and Sun:
    solar_phase_angle = np.arccos(np.dot(body_wrt_observer_unit, body_wrt_sun_unit))

    offset_dir = - (body_wrt_sun_unit - np.dot(body_wrt_sun_unit, body_wrt_observer_unit)
                    * body_wrt_observer_unit)
    # Calculate fraction of radius according to Fuentes-Munoz (2024):
    cot = lambda x: np.cos(x) / np.sin(x)
    num = 2 * (np.sin(solar_phase_angle) + (np.pi - solar_phase_angle) * np.cos(solar_phase_angle))
    denom = 3 * np.pi * (
            cot(solar_phase_angle / 2) - np.sin(solar_phase_angle / 2) * np.log(cot(solar_phase_angle / 4)))
    offset_ratio = num / denom # Fraction of body radius

    offset_mag = offset_ratio * (diameter/2) / norm(body_wrt_observer)

    return offset_dir * offset_mag


def photocenter_corrections_from_observations(
        observations: np.ndarray | ObservationCollection,
        diameter: float,
        body_name: str,
        observer_name : str,
        bodies: SystemOfBodies
) -> np.ndarray:
    """
    Calculate corrections to observations to account for the photocenter-barycenter offset. Should always be added
    to observations.

    Parameters
    ----------
    observations: np.ndarray | ObservationCollection
        Observations as [time, RA, DEC] with time in seconds since J2000, RA, DEC as radians or an ObservationCollection
    diameter : float
        Diameter of the body in meter
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_name : str
        Name of the observing body in the SystemOfBodies object
    bodies : SystemOfBodies
        The SystemOfBodies object that contains ephemerides of observer, body and sun in same reference frame


    Returns
    -------
    np.ndarray
        Corrections as [RA, DEC]
    """
    # Convert observation collection to an array
    if isinstance(observations, ObservationCollection):
        epochs = observations.get_concatenated_observation_times()
        observations = observations.get_concatenated_observations()
        ra, dec = observations[::2], observations[1::2]
        observations = np.column_stack((epochs, ra, dec))

    # Check inputs
    body_dne = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in (body_name, observer_name, 'Sun')]
    if any(body_dne):
        raise ValueError(f'Bodies {body_name}, {observer_name}, "Sun" not found in SystemOfBodies object, or their'
                         f'associated ephemerides do not exist')

    body_eph = bodies.get(body_name).ephemeris
    sun_eph = bodies.get('Sun').ephemeris
    observer_eph = bodies.get(observer_name).ephemeris

    frame_origins_match = (body_eph.frame_origin == sun_eph.frame_origin) and \
                          (observer_eph.frame_origin == sun_eph.frame_origin)
    frame_orientations_match = (body_eph.frame_orientation == sun_eph.frame_orientation) and \
                               (observer_eph.frame_orientation == sun_eph.frame_orientation)
    if not frame_orientations_match or not frame_origins_match:
        raise ValueError('Observer, body and sun ephemeris must be defined in the same reference frame')

    # Corrections to RA/DEC observations
    corrections = []

    for epoch, ra, dec in observations.T:

        # Calculate the offset as a vector in radians
        body_wrt_ssb = body_eph.cartesian_position(epoch)
        sun_wrt_ssb = sun_eph.cartesian_position(epoch)
        observer_wrt_ssb = observer_eph.cartesian_position(epoch)

        offset_vector = _photocenter_offset(diameter,
                                            (body_wrt_ssb - observer_wrt_ssb),
                                            (body_wrt_ssb - sun_wrt_ssb))

        # Calculate corrections to the astrometry
        corrections.append(
            _offset_vector_to_corrections(offset_vector, ra, dec)
        )

    corrections = np.array(corrections)

    return corrections

