from astroquery.jplsbdb import SBDB
from astropy import units as u
from typing import Union

class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body. 
    Usefull for retrieving names and masses in conjunction with the MPC module.
    """
    def __init__(self, MPCcode:Union[str, int]) -> None:
        """Create a Small-Body Database Query

        Parameters
        ----------
        MPCcode : Union[str, int]
            MPC code for the object.
        """
        self.query = SBDB.query(MPCcode)

    def __str__(self) -> str:
        return SBDB.schematic(self.query)
    
    @property
    def name(self):
        """Short name in the format: `MPC NAME DESIGNATION`"""
        return self.query["object"]["fullname"]
    
    @property
    def shortname(self):
        """Short name in the format: `MPC NAME`"""
        return self.query["object"]["shortname"]
    
    @property
    def spicename(self):
        """Returns spice kernel number for the codes_300ast_20100725.bsp spice kernel"""
        return self.query["object"]["spkid"]
    
    @property
    def gravitational_parameter(self):
        """Returns the gravitational parameter for the small body"""
        return self.query["phys_par"]["GM"].to(u.meter**3/u.second**2)