import typing
from _typeshed import Incomplete
from tudatpy.constants import GRAVITATIONAL_CONSTANT as GRAVITATIONAL_CONSTANT
from typing import Any

class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body.
    Usefull for retrieving names and masses in conjunction with the MPC module."""
    MPCcode: Incomplete
    query: Incomplete

    def __init__(self, MPCcode: str | int, *args, **kwargs) -> None:
        """Create a Small-Body Database Query

        Additional parameters are available through args and kwards, see:
        https://astroquery.readthedocs.io/en/latest/jplsbdb/jplsbdb.html

        Parameters
        ----------
        MPCcode : Union[str, int]
            MPC code for the object.

        """

    def __getitem__(self, index) -> Any:
        ...

    @property
    def name(self):
        """Short name in the format: `MPC NAME DESIGNATION`"""

    @property
    def shortname(self):
        """Short name in the format: `MPC NAME`"""

    @property
    def spkid(self):
        """Returns the JPL SPKID, the related codes_300_spkid method returns a modified ID for the Tudat standard kernel"""

    @property
    def codes_300_spkid(self):
        """Returns spice kernel number for the codes_300ast_20100725.bsp spice kernel.

        Some objects may return a name instead of a number.
        These are objects specifically specified by name in the codes_300ast_20100725.bsp kernel.

        See https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/aa_summaries.txt for a list of exceptions
        """

    @property
    def gravitational_parameter(self):
        """Returns the gravitational parameter for the small body if available"""

    @property
    def object_info(self):
        """Returns info about the object, including its designation and orbit class"""

    @property
    def object_classification(self):
        """Returns the orbit class of the object for example Amor Group for Eros"""

    @property
    def diameter(self):
        """Returns diameter of the small body if available"""

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