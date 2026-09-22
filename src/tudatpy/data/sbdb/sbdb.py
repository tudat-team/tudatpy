from astroquery.jplsbdb import SBDB as astroquerySBDB
from astropy import units as u
from astropy.time import Time
from typing import Any, Union
import math
import numpy as np
import datetime
from tudatpy.constants import GRAVITATIONAL_CONSTANT


class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body.
    Useful for retrieving names and masses in conjunction with the MPC module.
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
    def moid(self):
        """Returns the Small-Body's Minimum Orbital Intersection Distance with the Earth [m], if available"""
        return self.query["orbit"]["moid"].to(u.meter).value

    @property
    def moid_jup(self):
        """Returns the Small-Body's Minimum Orbital Intersection Distance with Jupiter [m], if available"""
        return self.query["orbit"]["moid_jup"].to(u.meter).value

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
            raise ValueError(f"Gravitational parameter is not available for object {self.name}")

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
            raise ValueError(f"Diameter is not available for object {self.name}")

    @property
    def nongrav_params(self):
        """Returns cometary non-gravitational model parameters for the small body in [m/s^2].
        If one or more parameter is unavailable a corresponding value of 0 is returned in the array
        """
        try:
            A1 = (
                (self.query["orbit"]["model_pars"]["A1"].value * u.au / u.day**2)
                .to(u.m / u.s**2)
                .value
            )
        except Exception as _:
            A1 = 0
        try:
            A2 = (
                (self.query["orbit"]["model_pars"]["A2"].value * u.au / u.day**2)
                .to(u.m / u.s**2)
                .value
            )
        except Exception as _:
            A2 = 0
        try:
            A3 = (
                (self.query["orbit"]["model_pars"]["A3"].value * u.au / u.day**2)
                .to(u.m / u.s**2)
                .value
            )
        except Exception as _:
            A3 = 0
        return np.array([A1, A2, A3])

    @property
    def Dt(self):
        """If available, returns asymmetric cometary non-gravitational model Dt (see Yeomans and Chodas, 1989) of the small body in [seconds]"""
        try:
            DT = self.query["orbit"]["model_pars"]["DT"].value
            return (
                (DT * u.day).to(u.s).value
            )  # a positive Dt will need to evaluate the position at (t-Dt)
        except Exception as _:
            raise ValueError(f"Asymmetry parameter DT is not available for object {self.name}")

    @property
    def first_obs(self):
        """If available, returns date of the first observation used for the orbit estimation of the small body in [JD]"""
        try:
            first_obs = [int(el) for el in self.query["orbit"]["first_obs"].split("-")]
            observation_start = datetime.datetime(first_obs[0], first_obs[1], first_obs[2])
            return Time(observation_start).jd
        except Exception as _:
            raise ValueError(f"Date of first observation is not available for object {self.name}")

    @property
    def last_obs(self):
        """If available, returns date of the last observation used for the orbit estimation of the small body in [JD]"""
        try:
            last_obs = [int(el) for el in self.query["orbit"]["last_obs"].split("-")]
            observation_end = datetime.datetime(last_obs[0], last_obs[1], last_obs[2])
            return Time(observation_end).jd
        except Exception as _:
            raise ValueError(f"Date of last observation is not available for object {self.name}")

    @property
    def perihelion(self):
        """If available, returns perihelion of the small body in [m]"""
        try:
            q = (self.query["orbit"]["elements"]["q"].value * u.au).to(u.m).value
            return q
        except Exception as _:
            raise ValueError(f"Perihelion distance is not available for object {self.name}")

    @property
    def time_perihelion(self):
        """If available, returns time of perihelion of the small body in [JD]"""
        try:
            tp = self.query["orbit"]["elements"]["tp"].value
            return tp
        except Exception as _:
            raise ValueError(f"Perihelion time is not available for object {self.name}")

    def estimated_spherical_mass(self, density: float) -> float:
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

    def estimated_spherical_gravitational_parameter(self, density: float) -> float:
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
