from astroquery.jplsbdb import SBDB as astroquerySBDB
from astropy import units as u
from os import PathLike
from typing import Any, Iterable, Optional, Union
from datetime import datetime
import math
import numpy as np
import pandas as pd
import requests
from tudatpy.astro.time_representation import DateTime
from tudatpy.constants import GRAVITATIONAL_CONSTANT


class SBDBbatch:
    """Access the JPL Small-Body Database catalogue as a table.

    Load a saved CSV file, or fetch all catalogue rows from JPL through
    requests of at most ``page_size`` rows each. Use :meth:`get`
    to select rows by primary MPC designation (the ``pdes`` field).

    Parameters
    ----------
    csv_path : str or os.PathLike, optional
        Path to a previously saved CSV file. When provided, the catalogue is
        read from this file instead of queried from JPL.
    fields : iterable of str, optional
        SBDB fields to request when downloading. If omitted, all fields
        advertised by the API are requested. Include ``"pdes"`` to use
        :meth:`get`.
    page_size : int, optional
        Maximum number of rows requested per JPL request. Default is 100,000.
    timeout : float, optional
        HTTP timeout in seconds for each request. Default is 120.

    Attributes
    ----------
    dataframe : pandas.DataFrame
        The downloaded or loaded catalogue.
    """

    api_url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"

    def __init__(
        self,
        csv_path: Optional[Union[str, PathLike]] = None,
        fields: Optional[Iterable[str]] = None,
        page_size: int = 100_000,
        timeout: float = 120,
    ) -> None:
        """Load the catalogue from JPL or a local CSV file."""
        self.dataframe = (
            pd.read_csv(csv_path, dtype={"pdes": str, "spkid": str})
            if csv_path is not None
            else self._download(fields, page_size, timeout)
        )

    @classmethod
    def _download(
        cls, fields: Optional[Iterable[str]], page_size: int, timeout: float
    ) -> pd.DataFrame:
        """Fetch all catalogue rows through successive API requests.

        Parameters
        ----------
        fields : iterable of str or None
            SBDB field names. If None, request all fields advertised by the API.
        page_size : int
            Maximum number of rows requested per API call.
        timeout : float
            HTTP timeout in seconds for each request.

        Returns
        -------
        pandas.DataFrame
            Catalogue rows with the columns returned by JPL.
        """
        if fields is None:
            response = requests.get(cls.api_url, params={"info": "field"}, timeout=timeout)
            response.raise_for_status()
            fields = [
                field["name"]
                for group in response.json()["info"]["field"].values()
                for field in group["list"]
            ]
        fields = list(fields)

        rows, offset, columns = [], 0, []
        while True:
            response = requests.get(
                cls.api_url,
                params={"fields": ",".join(fields), "limit": page_size, "limit-from": offset},
                timeout=timeout,
            )
            response.raise_for_status()
            result = response.json()
            columns = result["fields"]
            rows.extend(result.get("data", []))
            offset += len(result.get("data", []))
            if offset >= int(result["count"]) or not result.get("data"):
                break
        return pd.DataFrame(rows, columns=columns)

    def get(self, mpc_codes: Union[str, int, Iterable[Union[str, int]]]) -> pd.DataFrame:
        """Select catalogue rows by primary MPC designation.

        Parameters
        ----------
        mpc_codes : str, int, or iterable of str or int
            One or more designations to match against the ``pdes`` column.

        Returns
        -------
        pandas.DataFrame
            Matching rows, or an empty table if no designation matches.
        """
        codes = [mpc_codes] if isinstance(mpc_codes, (str, int)) else mpc_codes
        return self.dataframe[self.dataframe["pdes"].astype(str).isin(map(str, codes))]

    def to_csv(self, path: Union[str, PathLike]) -> None:
        """Save the loaded catalogue to a CSV file for offline reuse.

        Parameters
        ----------
        path : str or os.PathLike
            Output CSV path. The DataFrame index is not written.
        """
        self.dataframe.to_csv(path, index=False)


class SBDBquery:
    """Small-Body Database Query for retrieving various properties of a small body.
    Useful for retrieving names and masses in conjunction with the MPC module.
    """

    _AU_TO_M = (u.au, u.m)
    _KM_TO_M = (u.km, u.m)
    _KM_S_TO_M_S = (u.km / u.s, u.m / u.s)
    _AU_D2_TO_M_S2 = (u.au / u.d**2, u.m / u.s**2)

    # SBDB field -> (unit of the field in the SBDB API, tudat SI unit)
    # The SBDB units are those of the SBDB API (https://ssd-api.jpl.nasa.gov/doc/sbdb.html): documented there for
    # moid, moid_jup, q and the close approach fields, and given in the units entry of each API response for GM,
    # diameter, A1, A2, A3 and DT
    field_units = {
        "moid": _AU_TO_M,
        "moid_jup": _AU_TO_M,
        "q": _AU_TO_M,
        "dist": _AU_TO_M,
        "dist_min": _AU_TO_M,
        "dist_max": _AU_TO_M,
        "diameter": _KM_TO_M,
        "unc_major": _KM_TO_M,
        "unc_minor": _KM_TO_M,
        "GM": (u.km**3 / u.s**2, u.m**3 / u.s**2),
        "v_inf": _KM_S_TO_M_S,
        "v_rel": _KM_S_TO_M_S,
        "sigma_t": (u.min, u.s),
        "DT": (u.d, u.s),
        "unc_angle": (u.deg, u.rad),
        "A1": _AU_D2_TO_M_S2,
        "A2": _AU_D2_TO_M_S2,
        "A3": _AU_D2_TO_M_S2,
    }

    # Full body name (lower case) -> abbreviated name used by the SBDB close approach data
    close_approach_body_names = {
        "mercury": "Merc",
        "venus": "Venus",
        "earth": "Earth",
        "moon": "Moon",
        "mars": "Mars",
        "jupiter": "Juptr",
        "saturn": "Satrn",
        "uranus": "Urnus",
        "neptune": "Neptn",
        "pluto": "Pluto",
    }

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
        """**read-only**

        Short name in the format: `MPC NAME DESIGNATION`
        """
        return self.query["object"]["fullname"]

    @property
    def shortname(self):
        """**read-only**

        Short name in the format: `MPC NAME`
        """
        return self.query["object"]["shortname"]

    @property
    def spkid(self):
        """**read-only**

        Returns the JPL SPKID, the related codes_300_spkid method returns a modified ID for the Tudat standard kernel
        """
        return self.query["object"]["spkid"]

    @property
    def moid(self):
        """**read-only**

        Small body's minimum orbital intersection distance with Earth, in metres.
        """
        return self._to_si("moid", self.query["orbit"]["moid"])

    @property
    def moid_jup(self):
        """**read-only**

        Small body's minimum orbital intersection distance with Jupiter, in metres.
        """
        return self._to_si("moid_jup", self.query["orbit"]["moid_jup"])

    @property
    def codes_300_spkid(self):
        """**read-only**

        Returns spice kernel number for the codes_300ast_20100725.bsp spice kernel.

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
        """**read-only**

        Returns the gravitational parameter for the small body if available
        """
        try:
            return self._to_si("GM", self.query["phys_par"]["GM"])
        except Exception as _:
            raise ValueError(f"Gravitational parameter is not available for object {self.name}")

    @property
    def object_info(self):
        """**read-only**

        Returns info about the object, including its designation and orbit class
        """
        return self.query["object"]

    @property
    def object_classification(self):
        """**read-only**

        Returns the orbit class of the object for example Amor Group for Eros
        """
        return self.object_info["orbit_class"]["name"]

    @property
    def diameter(self):
        """**read-only**

        Returns diameter of the small body if available
        """
        try:
            return self._to_si("diameter", self.query["phys_par"]["diameter"])
        except Exception as _:
            raise ValueError(f"Diameter is not available for object {self.name}")

    @property
    def nongrav_params(self):
        """**read-only**

        Cometary non-gravitational model parameters ``A1``, ``A2`` and ``A3``
        in m/s². Missing parameters are returned as zero.
        """
        parameters = []
        for name in ("A1", "A2", "A3"):
            try:
                value = self._to_si(name, self.query["orbit"]["model_pars"][name])
            except Exception:
                value = 0.0
            parameters.append(value)
        return np.array(parameters)

    @property
    def Dt(self):
        """**read-only**

        Asymmetric cometary non-gravitational model parameter ``DT`` in seconds.
        """
        try:
            return self._to_si("DT", self.query["orbit"]["model_pars"]["DT"])
        except Exception as exception:
            raise ValueError(
                f"Asymmetry parameter DT is not available for object {self.name}"
            ) from exception

    @property
    def first_obs(self):
        """**read-only**

        Epoch of the first observation used in the orbit solution, in seconds since J2000.
        """
        try:
            observation_start = datetime.strptime(self.query["orbit"]["first_obs"], "%Y-%m-%d")
            return DateTime.from_python_datetime(observation_start).to_epoch()
        except Exception as exception:
            raise ValueError(
                f"Date of first observation is not available for object {self.name}"
            ) from exception

    @property
    def last_obs(self):
        """**read-only**

        Epoch of the last observation used in the orbit solution, in seconds since J2000.
        """
        try:
            observation_end = datetime.strptime(self.query["orbit"]["last_obs"], "%Y-%m-%d")
            return DateTime.from_python_datetime(observation_end).to_epoch()
        except Exception as exception:
            raise ValueError(
                f"Date of last observation is not available for object {self.name}"
            ) from exception

    @property
    def perihelion(self):
        """**read-only**

        Perihelion distance in metres.
        """
        try:
            return self._to_si("q", self.query["orbit"]["elements"]["q"])
        except Exception as exception:
            raise ValueError(
                f"Perihelion distance is not available for object {self.name}"
            ) from exception

    @property
    def time_perihelion(self):
        """**read-only**

        Perihelion epoch in seconds since J2000.
        """
        try:
            epoch_julian_day = self.query["orbit"]["elements"]["tp"].value
            return DateTime.from_julian_day(epoch_julian_day).to_epoch()
        except Exception as exception:
            raise ValueError(
                f"Perihelion time is not available for object {self.name}"
            ) from exception

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
        body: str | list[str],
        start_time: DateTime | None = None,
        end_time: DateTime | None = None,
        maximum_distance: float | None = None,
    ) -> dict:
        """Retrieve the close approaches of the small body to the requested bodies, as listed by the SBDB.
        Requires the SBDB query to be created with ``close_approach=True``.

        The quantities are converted from the SBDB units to the tudat conventions (SI units, epochs in seconds
        since J2000 TDB).

        Parameters
        ----------
        body : str | list[str]
            Name(s) of the bodies to consider (e.g. ``["Earth", "Jupiter"]``). Full planet names are mapped to the
            abbreviations used by the SBDB (e.g. ``"Jupiter"`` to ``"Juptr"``), other names are used as given
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
            - ``cd``: epoch of the closest approach, as calendar date string in TDB (as given by the SBDB, for display)
            - ``dist``, ``dist_min``, ``dist_max``: nominal, minimum and maximum (3-sigma) closest approach
              distance, in [m]
            - ``v_inf``, ``v_rel``: velocity relative to the body at infinity and at the closest approach, in [m/s]
            - ``sigma_t``: 3-sigma uncertainty of the closest approach epoch, in [s]
            - ``sigma_tf``: same as ``sigma_t``, as formatted string (as given by the SBDB)
            - ``unc_major``, ``unc_minor``: semi-major and semi-minor axes of the 3-sigma uncertainty ellipse in the
              target plane, in [m]
            - ``unc_angle``: orientation angle of the uncertainty ellipse in the target plane, in [rad]
            - ``orbit_ref``: reference of the orbit solution used by the SBDB, as string
            - ``body``: name of the body, as abbreviated by the SBDB
        """
        if "ca_data" not in self.query.keys():
            raise ValueError(
                "Query does not contain close approach data.\nTo request close approaches in the query, you must type:\n"
                "'sbdb_query = SBDBquery(str(asteroid_number), full_precision=True, close_approach=True)'\n"
            )

        close_approach_data = dict()
        for key, values in self.query["ca_data"].items():
            if key == "jd":
                # Epoch, not a duration: the origin is moved to J2000 as well (the SBDB epochs are in TDB)
                julian_days = np.atleast_1d(u.Quantity(values, u.d).value)
                close_approach_data["epoch"] = np.array(
                    [DateTime.from_julian_day(julian_day).to_epoch() for julian_day in julian_days]
                )
            elif key == "orbit_ref":
                # astroquery parses the orbit solution reference as a number (e.g. 220.0)
                close_approach_data[key] = np.array(
                    [
                        (
                            str(int(reference))
                            if isinstance(reference, float) and reference.is_integer()
                            else str(reference)
                        )
                        for reference in np.atleast_1d(values)
                    ]
                )
            elif key in self.field_units:
                close_approach_data[key] = np.atleast_1d(self._to_si(key, values))
            else:
                close_approach_data[key] = np.atleast_1d(np.asarray(values))

        # Select encounters involving any requested body
        bodies = [body] if isinstance(body, str) else body
        sbdb_bodies = [self.close_approach_body_names.get(name.lower(), name) for name in bodies]
        mask = np.isin(close_approach_data["body"], sbdb_bodies)

        # Filter by start time
        if start_time is not None:
            mask &= close_approach_data["epoch"] >= start_time.to_epoch()

        # Filter by end time
        if end_time is not None:
            mask &= close_approach_data["epoch"] <= end_time.to_epoch()

        # Filter by closest approach distance
        if maximum_distance is not None:
            mask &= close_approach_data["dist"] <= maximum_distance

        return {key: values[mask] for key, values in close_approach_data.items()}

    def _to_si(self, name, value):
        """Convert an SBDB field to SI units.

        astroquery attaches units only to some fields: ``dist`` and ``unc_angle`` always arrive as plain floats, and
        so does any close approach column with a missing entry. For these, the unit the SBDB API uses for the field
        is taken from ``field_units``. When astroquery does attach a unit, the conversion uses it, and fails if it
        is incompatible with the SI unit.
        """
        sbdb_unit, si_unit = self.field_units[name]
        if not isinstance(value, u.Quantity):
            value = np.asarray(value, dtype=float) * sbdb_unit
        return value.to(si_unit).value
