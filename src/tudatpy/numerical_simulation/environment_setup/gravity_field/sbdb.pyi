from astropy import units as u
import astroquery.jplsbdb.core
import math as math
from tudatpy.numerical_simulation.environment_setup.gravity_field import expose_gravity_field as gravity_field
import tudatpy.numerical_simulation.environment_setup.gravity_field.expose_gravity_field
import typing
from typing import Any
__all__ = ['Any', 'GRAVITATIONAL_CONSTANT', 'SBDBquery', 'astroquerySBDB', 'central_sbdb', 'central_sbdb_density', 'gravity_field', 'math', 'u']

class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body.
	Usefull for retrieving names and masses in conjunction with the MPC module.
	
	"""

    def __getitem__(self, index) -> typing.Any:
        ...

    def __init__(self, MPCcode: typing.Union[str, int], *args, **kwargs) -> None:
        """
        Create a Small-Body Database Query
        
                Additional parameters are available through args and kwards, see:
                https://astroquery.readthedocs.io/en/latest/jplsbdb/jplsbdb.html
        
                Parameters
                ----------
                MPCcode : Union[str, int]
                    MPC code for the object.
        
                
        """

    def __str__(self) -> str:
        ...

    def estimated_spherical_gravitational_parameter(self, density: float):
        """
        Calculate a very simple gravitational parameter by estimating the object's mass using a given density.
                Will raise an error if the body's diameter is not available on SBDB.
        
                Parameters
                ----------
                density : float
                    Density of the object in `kg m^-3`
        
                Returns
                -------
                float
                    Simplified estimation for the object's gravitational parameter
                
        """

    def estimated_spherical_mass(self, density: float):
        """
        Calculate a very simple mass by estimating the object's mass using a given density.
                Will raise an error if the body's diameter is not available on SBDB.
        
                Parameters
                ----------
                density : float
                    Density of the object in `kg m^-3`
        
                Returns
                -------
                float
                    Simplified estimation for the object's mass
                
        """

    @property
    def codes_300_spkid(self):
        """
        Returns spice kernel number for the codes_300ast_20100725.bsp spice kernel.
        
                Some objects may return a name instead of a number.
                These are objects specifically specified by name in the codes_300ast_20100725.bsp kernel.
        
                See https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/aa_summaries.txt for a list of exceptions
                
        """

    @property
    def diameter(self):
        """
        Returns diameter of the small body if available
        """

    @property
    def gravitational_parameter(self):
        """
        Returns the gravitational parameter for the small body if available
        """

    @property
    def name(self):
        """
        Short name in the format: `MPC NAME DESIGNATION`
        """

    @property
    def object_classification(self):
        """
        Returns the orbit class of the object for example Amor Group for Eros
        """

    @property
    def object_info(self):
        """
        Returns info about the object, including its designation and orbit class
        """

    @property
    def shortname(self):
        """
        Short name in the format: `MPC NAME`
        """

    @property
    def spkid(self):
        """
        Returns the JPL SPKID, the related codes_300_spkid method returns a modified ID for the Tudat standard kernel
        """

def central_sbdb(MPCcode: typing.Union[str, int]) -> tudatpy.numerical_simulation.environment_setup.gravity_field.expose_gravity_field.GravityFieldSettings:
    """Factory function to create central gravity field settings using a gravitational parameter retrieved from JPL's Small-Body Database (SBDB)
	
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
		>>> body_settings.get( "Eros" ).gravity_field_settings = environment_setup.gravity_field.central_sbdb("433")
	
	
	"""

def central_sbdb_density(MPCcode: typing.Union[str, int], density: float) -> tudatpy.numerical_simulation.environment_setup.gravity_field.expose_gravity_field.GravityFieldSettings:
    """Factory function to create a central gravity field settings using a diameter retrieved from JPL's Small-Body Database (SBDB) and an inputted density.
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
		>>> body_settings.get( "Eros" ).gravity_field_settings = environment_setup.gravity_field.central_sbdb_density("433", 2670.0)
	
	
	"""
GRAVITATIONAL_CONSTANT: float = 6.67259e-11
astroquerySBDB: astroquery.jplsbdb.core.SBDBClass