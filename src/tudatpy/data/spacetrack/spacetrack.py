import json
import getpass
import os
import requests
import requests_cache
from collections import defaultdict
import numpy as np
from tudatpy.dynamics import environment, environment_setup, propagation_setup
from datetime import datetime, timedelta
import math

class SpaceTrackQuery:

    def __init__(self, username: str | None = None, password: str | None = None) -> None:
        """
        Initializes the query client.

        This will:
        1. Set up a 1-hour cache for all API requests.
        2. Create a single, persistent session.
        3. Log that session in one time.
        """

        if not username:
            username = input('space-track username: ')
        if not password:
            password = getpass.getpass('space-track password: ')

        self.username: str = username
        """ Space-Track.org username """

        self.password: str = password
        """ Space-Track.org username """

        self.spacetrack_url: str = 'https://www.space-track.org'
        """ Space-Track.org url """

        # 1. Install a cache. All requests will be cached for 1 hour (3600s).
        # This respects Space-Track's TLE guidelines perfectly.
        requests_cache.install_cache('spacetrack_cache', expire_after=3600, backend='sqlite')

        # 2. Create a single, persistent session
        self.session: requests.Session = requests.Session()
        """ Space-Track.org Single Session """

        self._login()
        """ Space-Track.org login call """

        self.OMMUtils: SpaceTrackQuery.OMMUtils = self.OMMUtils(self)
        """ Tudatpy's OMM-Utils functions """

    def _login(self) -> None:
        """
        Logs into Space-Track.org using the provided username and password.
        Handles potential request exceptions and prints login status.
        """
        print("Logging into Space-Track...")
        try:
            response = self.session.post(
                self.spacetrack_url + "/ajaxauth/login",
                json={'identity': self.username, 'password': self.password}
            )
            response.raise_for_status()  # Raise an error for bad responses
            if response.status_code == 200:
                print("Login successful.")
            else:
                print(f"Login failed with status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f"Login failed: {e}")
            raise

    def _get_json_and_save(self, url: str, json_name: str) -> dict | list | None:
        """
        Retrieves JSON data from a given URL, saves it to a file, and returns the parsed JSON.

        This method handles API requests, caching, error checking, and file saving.

        Parameters
        ----------
        url : str
            The URL to fetch the JSON data from.
        json_name : str
            The name of the file to save the JSON data to.

        Returns
        -------
        dict or list or None
            The parsed JSON data if successful, otherwise None.
        """
        try:
            response = self.session.get(url)
            response.raise_for_status()  # Check for HTTP errors
            json_data = response.json()

            # Check if the request came from the cache
            if getattr(response, 'from_cache', False):
                print(f"Loaded data for {json_name} from cache.")
            else:
                print(f"Fetched new data for {json_name} from API.")

            # Save the file
            if os.path.exists(json_name):
                os.remove(json_name)
            with open(json_name, "w") as json_file:
                json.dump(json_data, json_file, indent=4)
                print(f"Data saved to {json_name} file.")

            return json_data

        except requests.exceptions.RequestException as e:
            print(f"Failed to retrieve data from the API. Error: {e}")
            return None
        except json.JSONDecodeError:
            print(f"Failed to decode JSON response. Content: {response.text[:200]}...")
            return None

    #####################################
    # --- AVAILABLE DOWNLOAD METHODS ---#
    #####################################

    def latest_on_orbit(self) -> dict | list | None:
        """
        Retrieves the newest propagable element set for all on-orbit payloads.

        Returns
        -------
        dict or list or None
            The parsed JSON data if successful, otherwise None.

        See Also
        --------
        _get_json_and_save : The internal method used to fetch and save the data.
        """
        json_name = "latest_on_orbit.json"
        url = os.path.join(
            self.spacetrack_url,
            'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/decay_date/null-val/epoch/>now-30/orderby/norad_cat_id/format/json'
        )
        return self._get_json_and_save(url, json_name)

    def descending_epoch(self, N: int | None = None) -> dict | list | None:
        """
        Retrieves GP data for payloads, ordered by most recent epoch.

        Parameters
        ----------
        N : int, optional
            The maximum number of results to return. If None, all available data is returned.

        Returns
        -------
        dict or list or None
            The parsed JSON data if successful, otherwise None.
        """
        if N is not None:
            json_name = f"gp_descending_{N}.json"
            url = os.path.join(
                self.spacetrack_url,
                f'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/orderby/epoch desc/limit/{N}/format/json'
            )
        else:
            json_name = "gp_descending.json"
            url = os.path.join(
                self.spacetrack_url,
                'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/orderby/epoch desc/format/json'
            )

        json_data = self._get_json_and_save(url, json_name)

        if json_data and N != 1:
            # Call the filtering function
            filtered_values = self.OMMUtils.filter_tles_keep_latest_creation_from_json(json_name)

            # Get the count of filtered TLEs
            filtered_count = len(list(filtered_values))

            # Rename the file based on the count of *filtered* values
            parts_json_name = json_name.split('.')[0].split('_')
            parts_json_name[-1] = f'{filtered_count}'
            new_name = '_'.join(parts_json_name) + '.json'

            # Rename the temp file to the new, descriptive name
            os.system(f'mv temp_filtered_tles.json {new_name}')

            if new_name != json_name:
                # Remove the original, unfiltered download
                os.system(f'rm {json_name}')

            print(f"Filtered TLEs downloaded to {new_name} file.")

        return json_data

    def get_tles_by_norad_ids(
            self, norad_ids: int | list[int], history: bool = False,
            orderby: str = 'epoch desc', limit_per_object: int = 1
    ) -> dict | list | None:
        """
        Retrieves TLEs for a list of NORAD IDs in a single batch query.

        Parameters
        ----------
        norad_ids : int or list[int]
            A single NORAD ID or a list of NORAD IDs (e.g., [25544, 25338]).
        history : bool, optional
            If True, queries `gp_history` to retrieve more than one TLE per object.
            Defaults to False.
        orderby : str, optional
            Specifies how to sort the results. Defaults to 'epoch desc'.
        limit_per_object : int, optional
            The maximum number of TLEs to retrieve per object. Defaults to 1.

        Returns
        -------
        dict or list or None
            The parsed JSON data if successful, otherwise None.
        """
        if not isinstance(norad_ids, (list, tuple, set)):
            norad_ids = [norad_ids] # Allow a single int/str to be passed

        # Convert list of IDs [25544, 25338] to string "25544,25338"
        id_string = ",".join(map(str, norad_ids))

        tle_class = 'gp'
        # By default, use the 'orderby' parameter
        use_orderby = True

        if history or limit_per_object > 1:
            tle_class = 'gp_history'
            # If we are querying gp_history for MULTIPLE objects,
            # the 'orderby' clause breaks the 'limit per object' logic.
            if len(norad_ids) > 1:
                print("Note: Removing 'orderby' for batch history query to ensure 'limit' applies per-object.")
                use_orderby = False

        if len(norad_ids) == 1:
            single_id = norad_ids[0]
            json_name = f"{tle_class}_{single_id}_{limit_per_object}.json"
        else:
            json_name = f"{tle_class}_{len(norad_ids)}_ids.json"

        # Dynamically Build URL
        url_parts = [
            'basicspacedata/query',
            f'class/{tle_class}',
            f'NORAD_CAT_ID/{id_string}'
        ]

        if use_orderby:
            url_parts.append(f'orderby/{orderby}')

        url_parts.extend([
            f'limit/{limit_per_object}',
            'format/json'
        ])

        url = os.path.join(self.spacetrack_url, "/".join(url_parts))
        json_data = self._get_json_and_save(url, json_name)

        # Run filtering logic if we queried history or asked for more than 1 TLE
        if json_data and (history or limit_per_object > 1):

            # Call the filtering function
            filtered_values = self.OMMUtils.filter_tles_keep_latest_creation_from_json(json_name)
            filtered_count = len(list(filtered_values)) # Get the count of filtered TLEs

            # Rename the file based on the count of *filtered* values
            parts_json_name = json_name.split('.')[0].split('_')
            parts_json_name[-1] = f'{filtered_count}'
            new_name = '_'.join(parts_json_name) + '.json'

            # Rename the temp file to a descriptive name
            os.system(f'mv temp_filtered_tles.json {new_name}')

            if new_name != json_name:
                # Remove the original, unfiltered download
                os.system(f'rm {json_name}')
            print(f"Filtered TLEs downloaded to {new_name} file.")

        return json_data

    def filtered_by_oe_dict(
            self, filter_oe_dict: dict[str, tuple[float | None, float | None]],
            limit: int = 100, output_file: str = 'filtered_results.json'
    ) -> dict | list | None:
        """
        Retrieves GP data for payloads, filtered by a dictionary of orbital elements.

        This method constructs a Space-Track.org query to filter TLEs based on
        specified ranges for orbital elements (e.g., semimajor axis, eccentricity).
        It supports filtering by minimum, maximum, or both for each element.

        Parameters
        ----------
        filter_oe_dict : dict
            A dictionary where keys are orbital element names (e.g., 'SEMIMAJOR_AXIS',
            'ECCENTRICITY', 'INCLINATION') and values are tuples or lists representing
            the bounds (min, max). If only one value is provided, it's treated as a minimum
            if it's the first element, or a maximum if it's the second.
            Example: {'SEMIMAJOR_AXIS': (7000, 8000), 'ECCENTRICITY': (None, 0.1)}
        limit : int, optional
            The maximum number of results to return. Defaults to 100.
        output_file : str, optional
            The name of the file to save the JSON data to. Defaults to 'filtered_results.json'.
        """

        base_url_parts = ['basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD']

        # Build the filter string
        for oe, bounds in filter_oe_dict.items():
            min_val = bounds[0] if len(bounds) > 0 else None
            max_val = bounds[1] if len(bounds) > 1 else None

            # Use the Space-Track range operator: 'min--max'
            if min_val is not None and max_val is not None:
                range_str = f"{min_val}--{max_val}"
                base_url_parts.extend([oe.upper(), range_str])
            elif min_val is not None:
                base_url_parts.extend([oe.upper(), f">{min_val}"])
            elif max_val is not None:
                base_url_parts.extend([oe.upper(), f"<{max_val}"])

        # Add sorting, limit, and format
        base_url_parts.extend([f"orderby/epoch desc/limit/{limit}/format/json"])
        # Join all parts with '/'
        url = os.path.join(self.spacetrack_url, "/".join(base_url_parts))
        print(f"Built single filter query: {url}")
        # This now makes ONE query instead of 2*N queries
        return self._get_json_and_save(url, output_file)

    #################
    # OMM Utilities #
    #################
    class OMMUtils:
        def __init__(self, parent: 'SpaceTrackQuery') -> None:
            #Make spacetrackquery parent class available to OMMUtils
            self.parent: SpaceTrackQuery = parent
        def save_batch_to_individual_files(
                self, json_data: list | dict, limit_per_object: int = 1
        ) -> list[str] | None:
            """
            Takes a list of TLEs from a batch query and saves them
            into individual files named in the 'gp_NORAD_CAT_ID_limit_per_object.json' format.

            This is particularly useful when `get_tles_by_norad_ids` is called with multiple
            NORAD IDs and `limit_per_object` > 1, as Space-Track returns all TLEs in a single
            JSON array, making it difficult to distinguish TLEs belonging to different objects
            without parsing the data.

            Parameters
            ----------
            json_data : list | dict
                The JSON data containing TLEs, typically obtained from a Space-Track query.
            limit_per_object : int, optional
                The limit of TLEs per object that was used in the original query. Defaults to 1.

            Returns
            -------
            list[str] or None
                A list of filenames where the individual TLEs were saved, or None if no data.
            """
            if not json_data:
                print("No data to split.")
                return None

            # Group TLEs by NORAD ID
            grouped_tles: dict[str, list] = defaultdict(list)
            for tle in json_data:
                grouped_tles[tle['NORAD_CAT_ID']].append(tle)

            # Save each group to its own file
            saved_files: list[str] = []
            for norad_id, tle_list in grouped_tles.items():
                # Use 'gp' for consistency, as that's the class
                filename = f"gp_{norad_id}_{limit_per_object}.json"

                # Save the list of TLEs for this ID
                with open(filename, "w") as f:
                    json.dump(tle_list, f, indent=4)
                saved_files.append(filename)

            print(f"Split batch data into {len(saved_files)} individual files.")
            return saved_files

        def get_norad_id_name_map(self, limit: int | None = None) -> dict[int, str]:
            """
            This function gives a norad_cat_id - object_name mapping.

            It queries the 'satcat' class from Space-Track.org to retrieve satellite catalog
            information and creates a dictionary mapping NORAD Catalog IDs to their corresponding
            object names. The data is cached and saved to a JSON file.

            Parameters
            ----------
            limit : int, optional
                The maximum number of SATCAT entries to retrieve. If None, all available data is returned.

            Returns
            -------
            dict[int, str]
                A dictionary where keys are NORAD Catalog IDs (int) and values are object names (str).
            Refactored to use the parent's session.
            """
            query_parts = ['basicspacedata/query/class/satcat/orderby/NORAD_CAT_ID/format/json']
            if limit:
                query_parts.append(f"limit/{limit}")

            url = os.path.join(self.parent.spacetrack_url, "/".join(query_parts))
            json_name = 'norad_id_to_name.json'

            # Use the parent's _get_json_and_save method for caching
            data = self.parent._get_json_and_save(url, json_name)

            if not data:
                raise RuntimeError("Failed to fetch SATCAT data.")

            map_data: dict[int, str] = {int(obj['NORAD_CAT_ID']): obj['OBJECT_NAME'] for obj in data if obj['OBJECT_NAME']}

            # The _get_json_and_save already saves, but this is for the map
            map_json_name = 'norad_id_to_name_map.json'
            with open(map_json_name, "w") as f:
                json.dump(map_data, f, indent=4)
                print(f"Saved {len(map_data)} entries to {map_json_name}.")

            return map_data

        def filter_tles_keep_latest_creation_from_json(self, json_filename: str) -> any:
            """
            Filters a list of TLEs from a JSON file, keeping only the latest TLE for each
            unique epoch based on the 'CREATION_DATE'.

            This function is useful when a Space-Track query returns multiple TLEs for the
            same epoch (e.g., due to updates or re-releases). It ensures that only the most
            recently created TLE for a given epoch is retained. The filtered TLEs are saved
            to a temporary file named "temp_filtered_tles.json".

            Parameters
            ----------
            json_filename : str
                The path to the JSON file containing the list of TLE objects.

            Returns
            -------
            any
                A view object containing the values (filtered TLE dictionaries) from the internal dictionary.
            # ... (this method is fine, no changes needed) ...
            """
            with open(json_filename, "r") as f:
                objects_list = json.load(f)

            filtered: dict[str, dict] = {}
            if objects_list and isinstance(objects_list[0], list):
                objects_list = objects_list[0]

            for i, object_ in enumerate(objects_list):
                epoch = object_['EPOCH']
                creation_date = datetime.fromisoformat(object_['CREATION_DATE'])

                if epoch not in filtered:
                    filtered[epoch] = object_
                else:
                    existing_creation_date = datetime.fromisoformat(filtered[epoch]['CREATION_DATE'])
                    if creation_date > existing_creation_date:
                        filtered[epoch] = object_

            with open("temp_filtered_tles.json", "w") as f_out:
                json.dump(list(filtered.values()), f_out, indent=4)

            return filtered.values()

        def get_tles(self, json_dict: list[dict[str, any]] | dict[str, any]) -> dict[str, tuple[str, str]]:
            """"
            Extracts TLE lines (TLE_LINE1 and TLE_LINE2) from a JSON dictionary or list
            of TLE objects and organizes them into a dictionary keyed by NORAD_CAT_ID.

            Parameters
            ----------
            json_dict : list | dict
                A JSON object (list of dictionaries or a single dictionary) containing TLE data.
                Each dictionary is expected to have 'NORAD_CAT_ID', 'TLE_LINE1', and 'TLE_LINE2' keys.

            Returns
            -------
            dict[str, tuple[str, str]]
                A dictionary where keys are NORAD Catalog IDs (str) and values are tuples of (TLE_LINE1, TLE_LINE2).
            """
            tle_dict: dict[str, tuple[str, str]] = {}
            if type(json_dict) is list:
                # Handle case where json_dict is a list of lists
                if json_dict and isinstance(json_dict[0], list):
                    json_dict = json_dict[0]
                for json_entry in json_dict:
                    norad_cat_id = json_entry['NORAD_CAT_ID']
                    tle_line_1 = json_entry['TLE_LINE1']
                    tle_line_2 = json_entry['TLE_LINE2']
                    tle_dict[norad_cat_id] = (tle_line_1, tle_line_2)
            return tle_dict

        def get_tudat_keplerian_element_set(self, json_dict: list | dict) \
                -> tuple[float | None, float | None, float | None, float | None, float | None, float | None]:
            """
            Extracts Keplerian orbital elements from a JSON dictionary or list of TLE objects
            and converts them into a format suitable for Tudat.

            The function assumes the input JSON contains at least one TLE object with keys
            like 'SEMIMAJOR_AXIS', 'ECCENTRICITY', 'INCLINATION', etc. It converts units
            (e.g., km to meters, degrees to radians) and calculates true anomaly from mean anomaly.

            Parameters
            ----------
            json_dict : list | dict
                A JSON object (list of dictionaries or a single dictionary) containing TLE data.
                Expected to have keys for semimajor axis, eccentricity, inclination,
                argument of pericenter, right ascension of ascending node, and mean anomaly.

            Returns
            -------
            tuple[float | None, float | None, float | None, float | None, float | None, float | None]
                A tuple containing (semimajor_axis, eccentricity, inclination, argument_of_pericenter,
                right_ascension_of_ascending_node, true_anomaly). All angular values are in radians,
                and semimajor axis is in meters. Returns a tuple of Nones if `json_dict` is empty.
            """
            if not json_dict:
                return (None,) * 6 # Return Nones if dict is empty

            # Handle list of lists
            if json_dict and isinstance(json_dict[0], list):
                json_dict = json_dict[0]

            first_obj = json_dict[0]

            a = float(first_obj['SEMIMAJOR_AXIS'])*1e3 # [meters]
            e =float(first_obj['ECCENTRICITY'])
            i = float(first_obj['INCLINATION']) * math.pi/180 # [rad]
            omega = float(first_obj['ARG_OF_PERICENTER']) * math.pi/180 # [rad]
            raan = float(first_obj['RA_OF_ASC_NODE']) * math.pi/180 # [rad]
            mo = float(first_obj['MEAN_ANOMALY']) * math.pi/180 # [rad]

            true_anomaly = self.mean_to_true_anomaly(mo,e)
            return a,e,i,omega,raan,true_anomaly

        def tle_to_TleEphemeris_object(self, tle_line_1: str, tle_line_2: str) -> environment.TleEphemeris:
            """
            Converts two TLE lines into a Tudat `TleEphemeris` object.

            This utility function takes the two-line element set strings and creates
            a `Tle` object, which is then used to instantiate a `TleEphemeris` object
            for use in Tudat simulations.

            Parameters
            ----------
            tle_line_1 : str
                The first line of the Two-Line Element set.
            tle_line_2 : str
                The second line of the Two-Line Element set.

            Returns
            -------
            environment.TleEphemeris
                A Tudat `TleEphemeris` object representing the ephemeris of the satellite.
            """
            object_tle = environment.Tle(tle_line_1, tle_line_2)
            ephemeris_object = environment.TleEphemeris("Earth", "J2000", object_tle, False)
            return ephemeris_object

        def plot_earth(self, ax: any, radius: float = 6378, color: str = 'lightblue', alpha: float = 0.5, resolution: int = 50) -> None:
            """
            Plots a spherical representation of Earth on a given Matplotlib 3D axes.

            This function is intended for visualization purposes, typically in conjunction
            with plotting satellite orbits.

            Parameters
            ----------
            ax : any
                The Matplotlib 3D axes object (`Axes3D`) to plot on.
            radius : float, optional
                The radius of the Earth in kilometers. Defaults to 6378 km.
            color : str, optional
                The color of the Earth sphere. Defaults to 'lightblue'.
            alpha : float, optional
                The transparency of the Earth sphere (0.0 to 1.0). Defaults to 0.5.
            resolution : int, optional
                The number of points used to generate the sphere, affecting its smoothness.
                Defaults to 50.
            """
            u = np.linspace(0, 2 * np.pi, resolution)
            v = np.linspace(0, np.pi, resolution)
            x = radius * np.outer(np.cos(u), np.sin(v))
            y = radius * np.outer(np.sin(u), np.sin(v))
            z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, rstride=1, cstride=1, color=color, alpha=alpha, edgecolor='none')

        def get_tle_reference_epoch(self, tle_line_1: str) -> datetime:
            """
            Extracts the reference epoch (time of applicability) from the first TLE line.

            The epoch is encoded in the TLE as a two-digit year and a day of the year
            (including fractional part). This function parses these values and converts
            them into a `datetime` object.

            Parameters
            ----------
            tle_line_1 : str
                The first line of the Two-Line Element set.

            Returns
            -------
            datetime
                A `datetime` object representing the reference epoch of the TLE.
            """
            year_str = tle_line_1[18:20]
            day_str = tle_line_1[20:32]
            year = int(year_str)
            year += 2000 if year < 57 else 1900
            day_of_year = float(day_str)
            epoch = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
            return epoch

        def mean_to_true_anomaly(self, mo: float, e: float, tol: float = 1e-8, max_iter: int = 100) -> float:
            """
            Converts mean anomaly to true anomaly using an iterative solver for Kepler's equation.

            This function solves Kepler's equation (M = E - e sin(E)) for the eccentric anomaly (E)
            given the mean anomaly (M) and eccentricity (e). It then uses the eccentric anomaly
            to calculate the true anomaly (ν).

            Parameters
            ----------
            mo : float
                Mean anomaly in radians.
            e : float
                Eccentricity (dimensionless).
            tol : float, optional
                The tolerance for convergence of Kepler's equation. Defaults to 1e-8.
            max_iter : int, optional
                The maximum number of iterations for solving Kepler's equation. Defaults to 100.

            Returns
            -------
            float
                True anomaly in radians.
            """

            # solve Kepler's equation
            E = mo if e < 0.8 else np.pi
            for _ in range(max_iter):
                f = E - e * np.sin(E) - mo
                f_prime = 1 - e * np.cos(E)
                E_new = E - f / f_prime
                if abs(E_new - E) < tol:
                    break
                E = E_new
            else:
                raise RuntimeError("Kepler's equation did not converge.")

            # compute true anomaly
            v = 2 * np.arctan2(
                np.sqrt(1 + e) * np.sin(E / 2),
                np.sqrt(1 - e) * np.cos(E / 2)
            )
            return v