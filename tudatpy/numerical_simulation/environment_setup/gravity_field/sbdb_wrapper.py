from tudatpy.data.sbdb import SBDBquery
from tudatpy.numerical_simulation.environment_setup import gravity_field
from typing import Union


def central_sbdb(MPCcode: Union[str, int]) -> gravity_field.GravityFieldSettings:
    """Factory function to create central gravity field settings from JPL's Small-Body Database (SBDB)

    JPL SBDB hosts information about small bodies such as asteroids, including the gravitational parameter.
    Please note that the gravitational parameter is not available for all objects and that the accuracy of this value varies per object

    This function is a wrapper for the tudatpy.data.sbdb functionality.
    That api is not available on the api documentation yet.
    For now, visit the HorizonsQuery souce code for extensive documentation: 
    https://github.com/tudat-team/tudatpy/blob/master/tudatpy/data/sbdb.py

    For more information on the JPL Small-Body Database, visit: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/

    
    Parameters
    ----------
    MPCcode : Union[str, int]
        MPC code for the object.

    Examples
    ----------
    Retrieve gravitational parameter settings for Eros

        >>> # create gravity field settings
        >>> body_settings.get( "Eros" ).gravity_field_settings = environment_setup.gravity_field.central_sbdb("433")

    """

    query = SBDBquery(MPCcode=MPCcode)
    settings = gravity_field.central(query.gravitational_parameter)
    return settings
