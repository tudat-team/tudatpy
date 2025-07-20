import typing
from .. import gravity_field as gravity_field
from tudatpy.data.sbdb import SBDBquery as SBDBquery

def central_sbdb(MPCcode: str | int) -> gravity_field.GravityFieldSettings:
    """Factory function to create central gravity field settings using a gravitational parameter retrieved from JPL\'s Small-Body Database (SBDB)
    
    JPL SBDB hosts information about small bodies such as asteroids, including the gravitational parameter for some objects.
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
        >>> body_settings.get( "Eros" ).gravity_field_settings = environment_setup.gravity_field.central_sbdb("433")"""

def central_sbdb_density(MPCcode: str | int, density: float) -> gravity_field.GravityFieldSettings:
    """Factory function to create a central gravity field settings using a diameter retrieved from JPL\'s Small-Body Database (SBDB) and an inputted density.
    A simplified gravitational parameter is calculated by assuming a spherical body and homogenous density.
    
    JPL SBDB hosts information about small bodies such as asteroids, including the diameter for some objects.
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
    density : float
        Mean density of the object in `kg m^-3`.
    
    Examples
    ----------
    Retrieve gravitational parameter settings for Eros
    
        >>> # create gravity field settings
        >>> body_settings.get( "Eros" ).gravity_field_settings = environment_setup.gravity_field.central_sbdb_density("433", 2670.0)"""