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
    def spkid(self):
        """Returns the JPL SPKID, the related codes_300_spkid method returns a modified ID for the Tudat standard kernel"""
        return self.query["object"]["spkid"]
    
    @property
    def codes_300_spkid(self):
        """Returns spice kernel number for the codes_300ast_20100725.bsp spice kernel. 

        Some Objects may return a name instead of a number.
        These are objects specifically specified by name in the codes_300ast_20100725.bsp kernel.
        
        See https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/aa_summaries.txt for a list of exceptions"""
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
        """Returns the gravitational parameter for the small body"""
        return self.query["phys_par"]["GM"].to(u.meter**3/u.second**2)