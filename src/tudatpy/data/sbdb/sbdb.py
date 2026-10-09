from astroquery.jplsbdb import SBDB as astroquerySBDB
from astropy import units as u
from typing import Any, Union
import math
import numpy as np
from datetime import datetime
from tudatpy.constants import GRAVITATIONAL_CONSTANT
from tudatpy.astro.time_representation import DateTime
from tudatpy.astro import time_representation


class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body.
    Useful for retrieving names and masses in conjunction with the MPC module.
    """

    def __init__(self, MPCcode: Union[str, int], *args, **kwargs) -> None:
        """Create a Small-Body Database Query

        Additional parameters are available through args and kwargs, see:
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
            A1 = (self.query["orbit"]["model_pars"]["A1"]).to(u.m / u.s**2).value
        except Exception as _:
            A1 = 0
        try:
            A2 = (self.query["orbit"]["model_pars"]["A2"]).to(u.m / u.s**2).value
        except Exception as _:
            A2 = 0
        try:
            A3 = (self.query["orbit"]["model_pars"]["A3"]).to(u.m / u.s**2).value
        except Exception as _:
            A3 = 0
        return np.array([A1, A2, A3])

    @property
    def Dt(self):
        """If available, returns asymmetric cometary non-gravitational model Dt (see Yeomans and Chodas, 1989) of the small body in [seconds]"""
        try:
            DT = (self.query["orbit"]["model_pars"]["DT"]).to(u.s).value
            return DT  # a positive Dt will need to evaluate the position at (t-Dt)
        except Exception as _:
            raise ValueError(f"Asymmetry parameter DT is not available for object {self.name}")

    @property
    def first_obs(self):
        """If available, returns date of the first observation used for the orbit estimation of the small body in [seconds since J2000]"""
        try:
            first_obs = self.query["orbit"]["first_obs"]
            observation_start = datetime.strptime(first_obs, "%Y-%m-%d")
            return DateTime.from_python_datetime(observation_start).to_epoch()
        except Exception as _:
            raise ValueError(f"Date of first observation is not available for object {self.name}")

    @property
    def last_obs(self):
        """If available, returns date of the last observation used for the orbit estimation of the small body in [seconds since J2000]"""
        try:
            last_obs = self.query["orbit"]["last_obs"]
            observation_end = datetime.strptime(last_obs, "%Y-%m-%d")
            return DateTime.from_python_datetime(observation_end).to_epoch()
        except Exception as _:
            raise ValueError(f"Date of last observation is not available for object {self.name}")

    @property
    def perihelion(self):
        """If available, returns perihelion of the small body in [m]"""
        try:
            q = (self.query["orbit"]["elements"]["q"]).to(u.m).value
            return q
        except Exception as _:
            raise ValueError(f"Perihelion distance is not available for object {self.name}")

    @property
    def time_perihelion(self):
        """If available, returns time of perihelion of the small body in [seconds since J2000]"""
        try:
            tp = self.query["orbit"]["elements"]["tp"].value
            tp = DateTime.from_julian_day(tp).to_epoch()
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

    def get_sbdb_close_approaches_with_bodies(
        self,
        body: list[str],
        start_time: DateTime | None = None,
        end_time: DateTime | None = None,
        maximum_distance: float | None = None,
    ) -> dict:
        """Retrieve the close approaches of the small body to the requested bodies, as listed by the SBDB.
        Requires the SBDB query to be created with `close_approach=True`.

        The quantities are converted from the SBDB units to the tudat conventions (SI units, epochs in seconds
        since J2000 TDB).

        Parameters
        ----------
        body : list[str]
            Names of the bodies to consider (e.g. ``["Earth", "Moon"]``)
        start_time : DateTime, optional
            Only keep close approaches at or after this epoch (TDB)
        end_time : DateTime, optional
            Only keep close approaches at or before this epoch (TDB)
        maximum_distance : float, optional
            Only keep close approaches whose nominal closest approach distance (``dist``) is at most this value, in [m]

        Returns
        -------
        dict
            Close approach data, with one array per quantity and one entry per close approach that satisfies all
            the filters:

            - ``epoch``: epoch of the closest approach, in [s] since J2000 TDB
            - ``jd``: epoch of the closest approach, as Julian day in TDB (as given by the SBDB)
            - ``cd``: epoch of the closest approach, as calendar date string in TDB (as given by the SBDB)
            - ``dist``, ``dist_min``, ``dist_max``: nominal, minimum and maximum (3-sigma) closest approach
              distance, in [m]
            - ``v_inf``, ``v_rel``: velocity relative to the body at infinity and at the closest approach, in [m/s]
            - ``sigma_t``: 3-sigma uncertainty of the closest approach epoch, in [s]
            - ``sigma_tf``: same as ``sigma_t``, as formatted string (as given by the SBDB)
            - ``unc_major``, ``unc_minor``: semi-major and semi-minor axes of the 3-sigma uncertainty ellipse in the
              target plane, in [m]
            - ``unc_angle``: orientation angle of the uncertainty ellipse in the target plane, in [rad]
            - ``orbit_ref``: reference of the orbit solution used by the SBDB
            - ``body``: name of the body
        """

        if "ca_data" not in self.query.keys():
            raise ValueError(
                "Query does not contain close approach data. \n To request close approaches in the query, you must type:\n"
                "'sbdb_query = SBDBquery(str(asteroid_number), full_precision=True, close_approach = True)'\n"
            )

        all_close_approach_data = self.query["ca_data"]

        # Convert SBDB units (au, km, km/s degrees) to Tudat units
        unit_conversions = {
            "dist": (u.au, u.m),
            "dist_min": (u.au, u.m),
            "dist_max": (u.au, u.m),
            "v_inf": (u.km / u.s, u.m / u.s),
            "v_rel": (u.km / u.s, u.m / u.s),
            "sigma_t": (u.min, u.s),
            "unc_major": (u.km, u.m),
            "unc_minor": (u.km, u.m),
            "unc_angle": (u.deg, u.rad),
        }

        close_approach_data = dict()
        for key, values in all_close_approach_data.items():
            if key in unit_conversions:
                sbdb_unit, tudat_unit = unit_conversions[key]
                if not isinstance(values, u.Quantity):
                    values = np.asarray(values, dtype=float) * sbdb_unit
                close_approach_data[key] = values.to(tudat_unit).value
            elif key == "jd":
                close_approach_data[key] = (
                    values.to(u.d).value
                    if isinstance(values, u.Quantity)
                    else np.asarray(values, dtype=float)
                )
            else:
                close_approach_data[key] = np.asarray(values)

        # Epochs of the closest approaches in seconds since J2000 TDB (the SBDB epochs are in TDB)
        close_approach_data["epoch"] = np.array(
            [
                time_representation.julian_day_to_seconds_since_epoch(jd)
                for jd in close_approach_data["jd"]
            ]
        )

        # Select encounters involving any requested body
        mask = np.isin(close_approach_data["body"], body)

        # Filter by start time
        if start_time is not None:
            mask &= close_approach_data["epoch"] >= start_time.to_epoch()

        # Filter by end time
        if end_time is not None:
            mask &= close_approach_data["epoch"] <= end_time.to_epoch()

        # Filter by closest approach distance
        if maximum_distance is not None:
            mask &= close_approach_data["dist"] <= maximum_distance

        wanted_approaches = {key: values[mask] for key, values in close_approach_data.items()}

        return wanted_approaches
