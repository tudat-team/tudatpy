"""
Functions to calculate photocenter corrections to observations
"""
import numpy as np
import pandas as pd
from tudatpy.data.sbdb import SBDBquery
from numpy.linalg import norm
from tudatpy.estimation.observable_models_setup.biases.observation_correction_utils import _offset_vector_to_corrections, _unit

def _retrieve_asteroid_diameter(mpc_number):
    """
    Retrieve asteroid diameter from JPL SBDB

    Args:
        mpc_number (int): MPC number of asteroid:

    Returns:
        float: Diameter in meters
    """
    # Query SBDB
    sbdb = SBDBquery(mpc_number)

    # Option 1: diameter is in SBDB
    try:
        diameter = sbdb.diameter  # Meters

    except ValueError:
        # Option 2: approximate from H magnitude
        try:
            magnitude_absolute = sbdb.query['phys_par']['H']
            diameter = 1329000 / np.sqrt(0.125) * 10 ** (-0.2 * magnitude_absolute)  # Meters
        # Assume zero diameter if no data available
        except ValueError:
            diameter = 0
            print(f'No diameter for object {mpc_number}')

    assert diameter >= 0
    return diameter

def _solar_phase_angle(asteroid_wrt_gaia_unit,
                       asteroid_wrt_sun_unit,):
    """
    Calculate solar phase angle in radians. Solar phase angle is angle between lines connecting observer, asteroid, sun

    Args:
        asteroid_wrt_gaia_unit: Unit vector of asteroid wrt. Gaia
        asteroid_wrt_sun_unit: Unit vector of asteroid wrt. sun

    Returns:
        float: Solar phase angle in radians
    """
    return np.arccos(np.dot(asteroid_wrt_gaia_unit, asteroid_wrt_sun_unit))


def _offset_magnitude(solar_phase_angle: float,
                      diameter: float,
                      asteroid_gaia_distance: float) -> float:
    """
    Calculate magnitude in radians of photocenter offset

    Args:
        solar_phase_angle (float): Solar phase angle in rad
        diameter (float): Diameter in meters
        asteroid_gaia_distance (float): Distance between asteroid and Gaia

    Returns:
        float: Magnitude in radians of photocenter offset
    """
    # Calculate fraction of radius according to Fuentes-Munoz (2024)
    cot = lambda x: np.cos(x) / np.sin(x)
    num = 2 * (np.sin(solar_phase_angle) + (np.pi - solar_phase_angle) * np.cos(solar_phase_angle))
    denom = 3 * np.pi * (
            cot(solar_phase_angle / 2) - np.sin(solar_phase_angle / 2) * np.log(cot(solar_phase_angle / 4)))
    offset_ratio = num / denom # Fraction of body radius
    assert 0 <= offset_ratio <= 1
    return offset_ratio * (diameter/2) / asteroid_gaia_distance


def _offset_direction(asteroid_wrt_gaia_unit: np.ndarray,
                      asteroid_wrt_sun_unit: np.ndarray,):
    """
    Calculate direction of the photocenter offset vector

    Args:
        asteroid_wrt_gaia_unit:  Unit vector of asteroid wrt. Gaia
        asteroid_wrt_sun_unit:  Unit vector of asteroid wrt. sun

    Returns:
        np.ndarray: Direction (unit) of photocenter offset vector
    """
    offset_dir = - (asteroid_wrt_sun_unit - np.dot(asteroid_wrt_sun_unit, asteroid_wrt_gaia_unit)
                    * asteroid_wrt_gaia_unit)
    return _unit(offset_dir)


def photocenter_offset_spherical(mpc_number: int,
                                 table: pd.DataFrame,
                                 bodies):
    """
    Estimate the photocenter offset assuming a spherical shape with isotropic scattering properties.

    Args:
        mpc_number (int): MPC number of asteroid to calculate corrections for
        table (pd.DataFrame): observation table dataframe containing observations for mpc_number
        bodies (SystemOfBodies): bodies object which contains Gaia, asteroid and solar ephemeris

    Returns:
        np.ndarray: RA/DEC corrections to be added to observations
    """
    # Take one observation per transit as geometry is approx constant over 39 sec
    table_reduced = table[table['number_mp'] == mpc_number].drop_duplicates(subset='transit_id', ignore_index=True)
    assert table_reduced['epoch'].is_monotonic_increasing # Ensure same order wrt original table
    assert not table_reduced.empty

    asteroid_ephemeris = bodies.get(str(mpc_number)).ephemeris # Should be wrt SSB
    sun_ephemeris = bodies.get('Sun').ephemeris # wrt SSB

    # Corrections to RA/DEC observations
    ra_corrections = []
    dec_corrections = []

    # Loop through 1 obs-per-transit table
    for row in table_reduced.itertuples():
        # Retrieve problem geometry
        asteroid_wrt_ssb = asteroid_ephemeris.cartesian_position(row.epoch)
        sun_wrt_ssb = sun_ephemeris.cartesian_position(row.epoch)
        gaia_wrt_ssb = np.array((row.x_gaia, row.y_gaia, row.z_gaia)) # Note: possible small inconsistency if geocentric gaia eph. is used

        asteroid_wrt_sun_unit = _unit(asteroid_wrt_ssb - sun_wrt_ssb)
        asteroid_wrt_gaia_unit = _unit(asteroid_wrt_ssb - gaia_wrt_ssb)
        asteroid_gaia_distance = norm(asteroid_wrt_ssb - gaia_wrt_ssb)

        # Calculate offset vector
        solar_phase_angle = _solar_phase_angle(asteroid_wrt_gaia_unit, asteroid_wrt_sun_unit)
        diameter = _retrieve_asteroid_diameter(mpc_number)
        offset_vec = (_offset_magnitude(solar_phase_angle, diameter,asteroid_gaia_distance)
                      * _offset_direction(asteroid_wrt_gaia_unit, asteroid_wrt_sun_unit))

        # Calculate corrections
        ra_corr, dec_corr = _offset_vector_to_corrections(offset_vec, row.ra, row.dec)

        # Apply corrections across all obs in a transit
        transit_length = np.count_nonzero(table['transit_id'] == row.transit_id)
        ra_corrections.extend([ra_corr] * transit_length)
        dec_corrections.extend([dec_corr] * transit_length)

    corrections = np.column_stack((ra_corrections, dec_corrections))
    assert corrections.shape == (np.count_nonzero(table['number_mp'] == mpc_number), 2)
    return corrections

