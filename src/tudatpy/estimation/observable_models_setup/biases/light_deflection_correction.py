"""
Functions to calculate light deflection corrections to observations
"""
import numpy as np
from numpy.linalg import norm
import pandas as pd
from tudatpy.estimation.observable_models_setup.biases.observation_correction_utils import _offset_vector_to_corrections
from tudatpy.constants import SPEED_OF_LIGHT

def _calculate_light_deflection(gaia_state,
                                asteroid_state,
                                body_state,
                                mu_body):
    """
    Helper function to calculate post-Newtonian light-bending contribution from one body (delta k_pn in Klioner (2003))

    Args:
        gaia_state (np.ndarray): SSB position of Gaia at time of observation
        asteroid_state (np.ndarray): SSB position of asteroid at time of emission
        body_state (np.ndarray): SSB position of light-deflecting body at reference time
        mu_body (float): Gravitational parameter of body

    Returns:
        np.ndarray: Individual contribution to total light deflection
    """
    # Define vectors
    r_upper = gaia_state - asteroid_state  # capital R
    r_ea = asteroid_state - body_state
    r_oa = gaia_state - body_state

    # Calculate contribution of deflection from body
    factor = (2 * mu_body) / (SPEED_OF_LIGHT ** 2)
    num = np.cross(r_upper, np.cross(r_ea, r_oa))
    denom = norm(r_upper) * norm(r_oa) * (norm(r_ea) * norm(r_oa) + np.dot(r_oa, r_ea))
    offset_vec = factor * num / denom
    assert offset_vec.shape == (3,)

    return offset_vec


def relativistic_light_deflection(mpc_number,
                                  table: pd.DataFrame,
                                  bodies,
                                  bodies_to_include: list = ['Sun']):
    """
    Function that calculates offset due to relativistic bending of light around massive bodies.
    Uses pre-loaded ephemerides of asteroids and major bodies to calculate deflection

    Args:
        mpc_number (int): Asteroid to calculate corrections for
        table (pd.DataFrame): Observation table that contains asteroids observations
        bodies (SystemOfBodies): Tudat bodies object. Must have all appropriate ephemerides loaded
        bodies_to_include (list): A list of bodies that exert light bending

    Returns:
        np.ndarray: Corrections to RA and DEC, to be added to observations
    """
    table = table[table['number_mp'] == mpc_number]
    table = table.reset_index(drop=True)
    assert table['epoch'].is_monotonic_increasing
    assert not table.empty, 'Observation table contains no observations'

    ra_corrections = []
    dec_corrections = []

    # Retrieve asteroid ephemeris
    asteroid_ephemeris = bodies.get(str(mpc_number)).ephemeris

    # Loop over observations
    for row in table.itertuples():
        current_epoch = row.epoch

        # Calculate light-time from observer to asteroid (one iteration)
        gaia_state = np.array((row.x_gaia, row.y_gaia, row.z_gaia))
        # NOTE: possible small offset if geocentric ephemeris is used in estimation
        light_time_observer_asteroid = np.linalg.norm(
            gaia_state - asteroid_ephemeris.cartesian_position(current_epoch)
        )
        light_time_observer_asteroid /= SPEED_OF_LIGHT

        # Retrieve asteroid state at time of emission
        asteroid_state = asteroid_ephemeris.cartesian_position(current_epoch - light_time_observer_asteroid)

        total_offset_vec = [] # Total deflection from all bodies
        # Loop over light-deflecting bodies
        for body in bodies_to_include:

            # Retrieve body environment models
            body_ephemeris = bodies.get(body).ephemeris
            mu_body = bodies.get(body).gravitational_parameter

            # Calculate reference time for light-deflecting body (one iteration)
            light_time_observer_body = np.linalg.norm(
                gaia_state - body_ephemeris.cartesian_position(current_epoch)
            )
            light_time_observer_body /= SPEED_OF_LIGHT

            # Retrieve body state at reference time
            body_state = body_ephemeris.cartesian_position(current_epoch - light_time_observer_body)

            offset_from_body = _calculate_light_deflection(gaia_state, asteroid_state,body_state, mu_body,)
            total_offset_vec.append(offset_from_body)

        # Add contributions from all bodies
        # Note: we omit the minus sign from Klioner eq. 70 since this is just an odd convention in observation directions
        total_offset_vec = np.sum(total_offset_vec, axis=0)
        assert total_offset_vec.shape == (3,)

        ra_corr, dec_corr = _offset_vector_to_corrections(total_offset_vec, row.ra, row.dec)

        # Save differences
        ra_corrections.append(ra_corr)
        dec_corrections.append(dec_corr)

    corrections = np.column_stack((ra_corrections, dec_corrections))
    assert corrections.shape == (len(table),2)

    return corrections













