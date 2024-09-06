from astroquery.jplsbdb import SBDB as astroquerySBDB
from astropy import units as u
from typing import Any, Union
import math
from ....constants import GRAVITATIONAL_CONSTANT
from . import expose_gravity_field as gravity_field


class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body.
    Usefull for retrieving names and masses in conjunction with the MPC module.
    """

    def __init__(self, MPCcode: Union[str, int], *args, **kwargs) -> None:
        """Create a Small-Body Database Query

        Additional parameters are available through args and kwards, see:
        https://astroquery.readthedocs.io/en/latest/jplsbdb/jplsbdb.html

        Parameters
        ----------
        MPCcode : Union[str, int]
            MPC code for the object.

        """
        self.MPCcode = MPCcode
        self.query = astroquerySBDB.query(MPCcode, phys=True, *args, **kwargs)

    def __str__(self) -> str:
        return astroquerySBDB.schematic(self.query)

    def __getitem__(self, index) -> Any:
        # makes the instance directly subscriptable, i.e.: query["phys_par"]
        return dict(self.query[index])

    @property
    def name(self):
        """Short name in the format: `MPC NAME DESIGNATION`"""
        return self.query["object"]["fullname"]

    @property
    def shortname(self):
        """Short name in the format: `MPC NAME`"""
        return self.query["object"]["shortname"]

    @property
    def spkid(self):
        """Returns the JPL SPKID, the related codes_300_spkid method returns a modified ID for the Tudat standard kernel"""
        return self.query["object"]["spkid"]

    @property
    def codes_300_spkid(self):
        """Returns spice kernel number for the codes_300ast_20100725.bsp spice kernel.

        Some objects may return a name instead of a number.
        These are objects specifically specified by name in the codes_300ast_20100725.bsp kernel.

        See https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/aa_summaries.txt for a list of exceptions
        """
        spkid = self.spkid[0] + self.spkid[2:]
        if spkid == "2000001":
            return "Ceres"
        elif spkid == "2000004":
            return "Vesta"
        elif spkid == "2000021":
            return "Lutetia"
        elif spkid == "2000216":
            return "Kleopatra"
        elif spkid == "2000433":
            return "Eros"
        else:
            return spkid

    @property
    def gravitational_parameter(self):
        """Returns the gravitational parameter for the small body if available"""
        try:
            res = self.query["phys_par"]["GM"].to(u.meter**3 / u.second**2)
            return res.value
        except Exception as _:
            raise ValueError(
                f"Gravitational parameter is not available for object {self.name}"
            )

    @property
    def object_info(self):
        """Returns info about the object, including its designation and orbit class"""
        return self.query["object"]

    @property
    def object_classification(self):
        """Returns the orbit class of the object for example Amor Group for Eros"""
        return self.object_info["orbit_class"]["name"]

    @property
    def diameter(self):
        """Returns diameter of the small body if available"""
        try:
            res = self.query["phys_par"]["diameter"].to(u.meter)
            return res.value
        except Exception as _:
            raise ValueError(
                f"Gravitational parameter is not available for object {self.name}"
            )

    def estimated_spherical_mass(self, density: float):
        """Calculate a very simple mass by estimating the object's mass using a given density.
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
        volume = (math.pi / 6) * self.diameter**3
        mass = volume * density
        return mass

    def estimated_spherical_gravitational_parameter(self, density: float):
        """Calculate a very simple gravitational parameter by estimating the object's mass using a given density.
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
        return GRAVITATIONAL_CONSTANT * self.estimated_spherical_mass(density)


def central_sbdb(MPCcode: Union[str, int]) -> gravity_field.GravityFieldSettings:
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

    query = SBDBquery(MPCcode=MPCcode)
    settings = gravity_field.central(query.gravitational_parameter)
    return settings


def central_sbdb_density(
    MPCcode: Union[str, int], density: float
) -> gravity_field.GravityFieldSettings:
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

    query = SBDBquery(MPCcode=MPCcode)
    gm = query.estimated_spherical_gravitational_parameter(density)
    settings = gravity_field.central(gm)
    return settings
