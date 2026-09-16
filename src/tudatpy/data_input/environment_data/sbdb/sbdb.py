from astroquery.jplsbdb import SBDB as astroquerySBDB
from astropy import units as u
from os import PathLike
from typing import Any, Iterable, Optional, Union
import math
import pandas as pd
import requests
from tudatpy.constants import GRAVITATIONAL_CONSTANT


class SBDBbatch:
    """Access the JPL Small-Body Database catalogue as a table.

    Download the catalogue in pages or load it from a CSV file. Use :meth:`get`
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
        Maximum number of rows requested per API page. Default is 100,000.
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
        """Download the requested catalogue fields in pages.

        Parameters
        ----------
        fields : iterable of str or None
            SBDB field names. If None, request all fields advertised by the API.
        page_size : int
            Maximum number of rows requested per page.
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
        return self.query["orbit"]["moid"].to(u.meter).value

    @property
    def moid_jup(self):
        """**read-only**

        Small body's minimum orbital intersection distance with Jupiter, in metres.
        """
        return self.query["orbit"]["moid_jup"].to(u.meter).value

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
            res = self.query["phys_par"]["GM"].to(u.meter**3 / u.second**2)
            return res.value
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
            res = self.query["phys_par"]["diameter"].to(u.meter)
            return res.value
        except Exception as _:
            raise ValueError(f"Gravitational parameter is not available for object {self.name}")

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
