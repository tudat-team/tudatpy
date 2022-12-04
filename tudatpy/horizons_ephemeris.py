import os

if bool(os.getenv("READTHEDOCS")):
    os.environ['PROJ_LIB'] = "/home/docs/checkouts/readthedocs.org/user_builds/tudatpy/share/proj"
else:
    os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + '/share/proj'
import numpy as np
import astropy

from astropy.time import Time
from astroquery.jplhorizons import Horizons
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice

def create_ephemeris_dictionary_from_horizons_table(
        horizons_table ):

    """Creates a dictionary of times/states from a Horizons Table output (see https://astroquery.readthedocs.io/en/latest/api/astroquery.jplhorizons.HorizonsClass.html#astroquery.jplhorizons.HorizonsClass.vectors), to be used as input for tabulated ephemeris settings.

    Parameters
    ----------
    horizons_table : astropy.table.Table
        Table of states, as directly retried from JPL Horizons through astroquery (typically using the jpl_horizons_ephemeris_table function)

    Returns
    -------
    Dictionary (times as key; Cartesian states as value) storing a body state history, as retrieved from JPL Horizons

    """

    # Create dictionary, and iterate over horizons output table	
    ephemeris_table = dict( )
    for i in range(len(horizons_table)):

        # Conert time format to Tudat format
        current_time_jd_tdb = Time( horizons_table['datetime_jd'][ i ], format = 'jd', scale='utc')
        julian_days_since_j2000 = current_time_jd_tdb.jd1 - constants.JULIAN_DAY_ON_J2000
        current_time = ( julian_days_since_j2000 + current_time_jd_tdb.jd2 ) * constants.JULIAN_DAY
        
        # Convert state to SI units, and store in dictionary
        ephemeris_table[ current_time ] = np.array([
            horizons_table['x'][i] * constants.ASTRONOMICAL_UNIT,
            horizons_table['y'][ i ] * constants.ASTRONOMICAL_UNIT,
            horizons_table['z'][ i ] * constants.ASTRONOMICAL_UNIT,
            horizons_table['vx'][i] * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY_ON_J2000,
            horizons_table['vy'][ i ] * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY_ON_J2000,
            horizons_table['vz'][ i ] * constants.ASTRONOMICAL_UNIT / constants.JULIAN_DAY_ON_J2000
        ])

    return ephemeris_table


def jpl_horizons_ephemeris_table(
        body_identifier,
        central_body_identifier,
        epochs ):

    """Creates a Horizons Table (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html), and converts it to a dictionary of times/states, to be used as input for tabulated ephemeris settings.

    Parameters
    ----------
    body_identifier : str
        Horizons-compatible string identfying the body for which stats are to be generated
    central_body_identifier : str
	Horizons-compatible string identfying the central body (observer) w.r.t. which stats are to be generated
    epochs: list[float]
        List of seconds since J2000 (TDB) at which Horizons is to be queried. Note that the conversion to julian (as needed by Horizons) is done inside this funcion

    Returns
    -------
    Dictionary (times as key; Cartesian states as value) storing a body state history, as retrieved from JPL Horizons

    """

    # Convert time from seconds since J2000 to julian days
    epochs_julian_days = [ epoch / constants.JULIAN_DAY + constants.JULIAN_DAY_ON_J2000 for epoch in epochs]

    # Query Horizons for the required body and time interval
    obj = Horizons(id=body_identifier, location='@'+central_body_identifier,epochs=epochs_julian_days)

    # Extract the table of Cartesian states
    return create_ephemeris_dictionary_from_horizons_table(obj.vectors())



def jpl_horizons_tabulated_ephemeris_settings(
        body_identifier,
        central_body_identifier,
        epochs ):

    """Uses JPL Horizons (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html), to create tabulated ephemeris settings.

    Parameters
    ----------
    body_identifier : str
        Horizons-compatible string identfying the body for which stats are to be generated
    central_body_identifier : str
	Horizons-compatible string identfying the central body (observer) w.r.t. which stats are to be generated
    epochs: list[float]
        List of seconds since J2000 (TDB) at which Horizons is to be queried. Note that the conversion to julian (as needed by Horizons) is done inside this funcion


    Returns
    -------
    ``TabulatedEphemerisSettings`` containing data retrieved from JPL horizons, for the given bodies and times, for the creation of a tabulated ephemeris

    """
    # Create ephemeris settings from table
    ephemeris_settings = environment_setup.ephemeris.tabulated(
        jpl_horizons_ephemeris_table( body_identifier, central_body_identifier, epochs ),
        central_body_identifier )
    return ephemeris_settings

