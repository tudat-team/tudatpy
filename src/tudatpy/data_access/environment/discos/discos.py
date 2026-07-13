import requests
from typing import Any


class DiscosQuery:
    """
    A class to query the ESA DISCOS (Database and Information System Characterising Objects in Space) API.
    """

    def __init__(
        self, token: str, timeout: int | None = 10, url: str = "https://discosweb.esoc.esa.int"
    ):
        """
        Initialize the DiscosQuery client.

        Parameters
        ----------
        token : str
            API authentication token from DISCOS.
        timeout : int, optional
            Request timeout in seconds. Defaults to 10.
        url : str, optional
            Base URL for the DISCOS API. Defaults to "https://discosweb.esoc.esa.int".
        """
        self.token = token
        self.url = url
        self.api_version = "2"
        self.headers = {
            "Authorization": f"Bearer {self.token}",
            "DiscosWeb-Api-Version": self.api_version,
        }
        self.timeout = timeout

    def query_object(
        self, sat_id: str | int, is_discos_id: bool = False, verbose: bool = True
    ) -> dict[str, Any] | None:
        """
        Query object information from DISCOS using either a NORAD ID or a DISCOS ID.

        Parameters
        ----------
        sat_id : str | int
            The satellite identifier (NORAD ID or DISCOS ID).
        is_discos_id : bool, optional
            If True, treats sat_id as a DISCOS internal ID. If False, treats it as a NORAD ID. Defaults to False.
        verbose : bool, optional
            If True, prints status and error messages to the console. Defaults to True.

        Returns
        -------
        dict[str, Any] | None
            A dictionary containing object attributes if successful, None otherwise.
        """
        if is_discos_id:
            query_url = f"{self.url}/api/objects/{sat_id}"
        else:
            query_url = f"{self.url}/api/objects?filter=eq(satno,{sat_id})"

        try:
            response = requests.get(query_url, headers=self.headers, timeout=self.timeout)
            response.raise_for_status()
            data = response.json().get("data")

            if not data:
                if verbose:
                    print(f"No object found for satellite {sat_id}.")
                return None

            # Handle response structure (List vs Single Object)
            attributes = data[0]["attributes"] if isinstance(data, list) else data["attributes"]

            if verbose:
                print(f"Successfully retrieved data for {sat_id}")

            return attributes

        except requests.exceptions.Timeout:
            if verbose:
                print(f"Timeout Error: The request for {sat_id} took too long.")
            return None
        except requests.exceptions.RequestException as e:
            # Catches ConnectionErrors, HTTP Errors, etc.
            if verbose:
                print(f"Network Error for {sat_id}: {e}")
            return None
        except (KeyError, TypeError) as e:
            if verbose:
                print(f"Data Parsing Error: API response structure changed. {e}")
            return None
