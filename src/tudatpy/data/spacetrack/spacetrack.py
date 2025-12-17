import json
import getpass
import os
import requests
from collections import defaultdict
import numpy as np
from tudatpy.dynamics import environment, environment_setup, propagation_setup
from datetime import datetime, timedelta
import math
from tudatpy.data import get_resource_path

class SpaceTrackQuery:
    """
    A class to handle queries to Space-Track.org API for retrieving TLEs and other space data.
    It manages authentication, session persistence, and local caching of TLE files to minimize API usage.
    """

    def __init__(self, username: str | None = None, password: str | None = None, tle_data_folder: str = get_resource_path() + "/tle_data") -> None:
        """
        Initializes the query client.

        This will:
        1. Create a single, persistent session.
        2. Log that session in one time.
        3. Initialize local TLE storage folder.

        Parameters
        ----------
        username : str | None, optional
            Space-Track.org username. If None, prompts user for input.
        password : str | None, optional
            Space-Track.org password. If None, prompts user for input.
        tle_data_folder : str, optional
            Path to the local folder where TLE files will be stored. Defaults to "tle_data".
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

        # Setup local TLE data folder
        self.tle_data_folder = tle_data_folder
        if not os.path.exists(self.tle_data_folder):
            os.makedirs(self.tle_data_folder)

        # 1. Create a single, persistent session
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

        Raises
        ------
        requests.exceptions.RequestException
            If the login request fails due to network issues or invalid credentials.
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
        Retrieves JSON data. Checks local folder FIRST.
        If missing, fetches from API, saves to file, and returns parsed JSON.

        Parameters
        ----------
        url : str
            The URL to fetch the JSON data from.
        json_name : str
            The filename (not path) to check/save in the local TLE folder.

        Returns
        -------
        dict | list | None
            The parsed JSON data if successful, otherwise None.
        """
        filepath = os.path.join(self.tle_data_folder, json_name)

        # 1. Check Local Cache
        if os.path.exists(filepath):
            print(f"Found local file: {filepath}. Skipping API call.")
            try:
                with open(filepath, 'r') as f:
                    return json.load(f)
            except json.JSONDecodeError:
                print(f"Local file {filepath} is corrupted. Re-downloading.")

        # 2. Fetch from API
        try:
            print(f"Fetching new data for {json_name} from API...")
            response = self.session.get(url)
            response.raise_for_status()
            json_data = response.json()

            # Save the file to the cache folder
            with open(filepath, "w") as json_file:
                json.dump(json_data, json_file, indent=4)
                print(f"Data saved to {filepath}.")

            return json_data

        except requests.exceptions.RequestException as e:
            print(f"Failed to retrieve data from the API. Error: {e}")
            return None
        except json.JSONDecodeError:
            print(f"Failed to decode JSON response. Content: {response.text[:200]}...")
            return None

    def _retrieve_local_tle_data(self, norad_cat_id: str | int, start_date: str, end_date: str) -> list[dict] | None:
        """
        Checks for a file covering the specific ID and Date Range.
        Pattern: tle_{ID}_{START}_{END}.json

        Parameters
        ----------
        norad_cat_id : str | int
            The NORAD ID of the object.
        start_date : str
            The start date of the requested range (YYYY-MM-DD).
        end_date : str
            The end date of the requested range (YYYY-MM-DD).

        Returns
        -------
        list[dict] | None
            A list of filtered TLE dictionaries if a covering file is found locally, otherwise None.
        """
        norad_cat_id = str(norad_cat_id)
        tle_file_name_pattern = f"tle_{norad_cat_id}_"

        try:
            req_start = datetime.strptime(start_date, "%Y-%m-%d")
            req_end = datetime.strptime(end_date, "%Y-%m-%d")
        except ValueError as e:
            print(f"Date parsing error: {e}")
            return None

        if not os.path.exists(self.tle_data_folder):
            return None

        for file_name in os.listdir(self.tle_data_folder):
            if file_name.startswith(tle_file_name_pattern) and file_name.endswith(".json"):
                try:
                    # Parse: tle_12345_2023-01-01_2023-01-10.json
                    parts = file_name.split('_')
                    if len(parts) < 4: continue

                    f_start_str = parts[2]
                    f_end_str = parts[3].split('.')[0]

                    f_start = datetime.strptime(f_start_str, "%Y-%m-%d")
                    f_end = datetime.strptime(f_end_str, "%Y-%m-%d")

                    # Check if file fully covers the requested range
                    if f_start <= req_start and f_end >= req_end:
                        filepath = os.path.join(self.tle_data_folder, file_name)
                        print(f"Found valid local coverage: {file_name}")
                        with open(filepath, 'r') as f:
                            data = json.load(f)
                            # Filter explicitly for the requested range
                            filtered = [
                                d for d in data
                                if start_date <= d['EPOCH'][:10] <= end_date
                            ]
                            return filtered
                except (ValueError, IndexError):
                    continue

        return None

    def get_tles_for_date_range(self, norad_id: int | str, start_date: str, end_date: str) -> list[dict] | None:
        """
        Primary method for simulation workflows.
        Checks local file cache for date coverage before hitting API.

        Parameters
        ----------
        norad_id : int | str
            The NORAD ID of the object.
        start_date : str
            The start date of the range (YYYY-MM-DD).
        end_date : str
            The end date of the range (YYYY-MM-DD).

        Returns
        -------
        list[dict] | None
            A list of TLE dictionaries covering the requested range, or None if retrieval fails.
        """
        # 1. Smart Range Check
        local_tles = self._retrieve_local_tle_data(norad_id, start_date, end_date)
        if local_tles is not None:
            return local_tles

        # 2. Query Space-Track
        print(f"No local file covers range {start_date} to {end_date} for {norad_id}. Querying API...")

        url = os.path.join(
            self.spacetrack_url,
            f"basicspacedata/query/class/gp/NORAD_CAT_ID/{norad_id}/EPOCH/{start_date}--{end_date}/orderby/EPOCH asc/format/json"
        )

        filename = f"tle_{norad_id}_{start_date}_{end_date}.json"

        # Use centralized saver (which handles saving to folder)
        # Note: We pass the URL. Since we know the file didn't exist (from step 1 check),
        # _get_json_and_save will fetch and save it.
        return self._get_json_and_save(url, filename)

    #####################################
    # --- AVAILABLE DOWNLOAD METHODS ---#
    #####################################

    def latest_on_orbit(self) -> dict | list | None:
        """
        Retrieves the newest propagable element set for all on-orbit payloads.
        Checks 'latest_on_orbit.json' in local folder first.

        Returns
        -------
        dict | list | None
            The parsed JSON data containing the latest TLEs, or None if retrieval fails.
        """
        json_name = "latest_on_orbit.json"
        url = os.path.join(
            self.spacetrack_url,
            'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/decay_date/null-val/epoch/>now-30/orderby/norad_cat_id/format/json'
        )
        return self._get_json_and_save(url, json_name)

    def descending_epoch(self, N: int | None = None) -> dict | list | None:
        """
        Retrieves GP data ordered by epoch. Checks local folder first.

        Parameters
        ----------
        N : int | None, optional
            Limit the number of results returned. If None, no limit is applied.

        Returns
        -------
        dict | list | None
            The parsed JSON data containing TLEs ordered by descending epoch, or None if retrieval fails.
        """
        if N is not None:
            json_name = f"gp_descending_limit_{N}.json"
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

        # Optional: Filter duplicates if new data was fetched
        if json_data and N != 1:
            filepath = os.path.join(self.tle_data_folder, json_name)

            # Perform filtering
            filtered_values = list(self.OMMUtils.filter_tles_keep_latest_creation_from_json(filepath))

            # Save FILTERED data back to the SAME filename to ensure cache stability
            with open(filepath, "w") as f:
                json.dump(filtered_values, f, indent=4)
            print(f"Filtered and updated local file: {json_name}")

        return json_data

    def get_tles_by_norad_ids(
            self, norad_ids: int | list[int], history: bool = False,
            orderby: str = 'epoch desc', limit_per_object: int = 1
    ) -> dict | list | None:
        """
        Retrieves TLEs for specific IDs. Checks local folder first.

        Parameters
        ----------
        norad_ids : int | list[int]
            A single NORAD ID or a list of NORAD IDs to query.
        history : bool, optional
            If True, queries the 'gp_history' class instead of 'gp'. Defaults to False.
        orderby : str, optional
            Ordering criteria for the query results. Defaults to 'epoch desc'.
        limit_per_object : int, optional
            Limit the number of TLEs returned per object. Defaults to 1.

        Returns
        -------
        dict | list | None
            The parsed JSON data containing the requested TLEs, or None if retrieval fails.
        """
        if not isinstance(norad_ids, (list, tuple, set)):
            norad_ids = [norad_ids]

        id_string = ",".join(map(str, norad_ids))
        tle_class = 'gp'
        use_orderby = True

        if history or limit_per_object > 1:
            tle_class = 'gp_history'
            if len(norad_ids) > 1:
                print("Note: Removing 'orderby' for batch history query.")
                use_orderby = False

        # Stable Filename Generation
        if len(norad_ids) == 1:
            single_id = norad_ids[0]
            # Use query params in filename to make it unique
            json_name = f"{tle_class}_{single_id}_limit_{limit_per_object}.json"
        else:
            # Hash or summary for multiple IDs could be used, but simplified here:
            json_name = f"{tle_class}_batch_{len(norad_ids)}_ids_limit_{limit_per_object}.json"

        # Build URL
        url_parts = ['basicspacedata/query', f'class/{tle_class}', f'NORAD_CAT_ID/{id_string}']
        if use_orderby:
            url_parts.append(f'orderby/{orderby}')
        url_parts.extend([f'limit/{limit_per_object}', 'format/json'])

        url = os.path.join(self.spacetrack_url, "/".join(url_parts))

        # 1. Fetch (or Load Local)
        json_data = self._get_json_and_save(url, json_name)

        # 2. Post-Process (Filter Duplicates)
        # Only needed if we suspect duplicates and we just downloaded it.
        # But for consistency, we can re-filter or ensure the saved file is clean.
        if json_data and (history or limit_per_object > 1):
            filepath = os.path.join(self.tle_data_folder, json_name)
            filtered_values = list(self.OMMUtils.filter_tles_keep_latest_creation_from_json(filepath))

            # Save FILTERED data back to SAME filename
            with open(filepath, "w") as f:
                json.dump(filtered_values, f, indent=4)
            print(f"Filtered local file content: {json_name}")

        return json_data

    def filtered_by_oe_dict(
            self, filter_oe_dict: dict[str, tuple[float | None, float | None]],
            limit: int = 100, output_file: str = 'filtered_results.json'
    ) -> dict | list | None:
        """
        Retrieves TLEs by orbital elements. Checks local folder first.

        Parameters
        ----------
        filter_oe_dict : dict[str, tuple[float | None, float | None]]
            A dictionary where keys are orbital element names (e.g., 'SEMIMAJOR_AXIS') and values are tuples of (min, max).
            Use None for unbound limits.
        limit : int, optional
            Limit the number of results returned. Defaults to 100.
        output_file : str, optional
            Filename to save the results to locally. Defaults to 'filtered_results.json'.

        Returns
        -------
        dict | list | None
            The parsed JSON data containing the filtered TLEs, or None if retrieval fails.
        """
        base_url_parts = ['basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD']

        for oe, bounds in filter_oe_dict.items():
            min_val = bounds[0] if len(bounds) > 0 else None
            max_val = bounds[1] if len(bounds) > 1 else None

            if min_val is not None and max_val is not None:
                range_str = f"{min_val}--{max_val}"
                base_url_parts.extend([oe.upper(), range_str])
            elif min_val is not None:
                base_url_parts.extend([oe.upper(), f">{min_val}"])
            elif max_val is not None:
                base_url_parts.extend([oe.upper(), f"<{max_val}"])

        base_url_parts.extend([f"orderby/epoch desc/limit/{limit}/format/json"])
        url = os.path.join(self.spacetrack_url, "/".join(base_url_parts))

        return self._get_json_and_save(url, output_file)

    #################
    # OMM Utilities #
    #################
    class OMMUtils:
        """
        Helper class for handling OMM (Orbit Mean-Elements Message) data and TLE manipulations.
        """
        def __init__(self, parent: 'SpaceTrackQuery') -> None:
            """
            Initializes the OMMUtils class.

            Parameters
            ----------
            parent : SpaceTrackQuery
                The parent SpaceTrackQuery instance.
            """
            self.parent: SpaceTrackQuery = parent

        def save_batch_to_individual_files(
                self, json_data: list | dict, limit_per_object: int = 1
        ) -> list[str] | None:
            """
            Splits batch TLE data into individual files per NORAD ID.

            Parameters
            ----------
            json_data : list | dict
                The batch JSON data retrieved from Space-Track.
            limit_per_object : int, optional
                The limit used in the query, for filename generation. Defaults to 1.

            Returns
            -------
            list[str] | None
                A list of created filenames, or None if no data provided.
            """
            if not json_data:
                print("No data to split.")
                return None

            grouped_tles: dict[str, list] = defaultdict(list)
            for tle in json_data:
                grouped_tles[tle['NORAD_CAT_ID']].append(tle)

            saved_files: list[str] = []
            for norad_id, tle_list in grouped_tles.items():
                filename = f"gp_{norad_id}_limit_{limit_per_object}.json"
                filepath = os.path.join(self.parent.tle_data_folder, filename)

                with open(filepath, "w") as f:
                    json.dump(tle_list, f, indent=4)
                saved_files.append(filename)

            print(f"Split batch data into {len(saved_files)} individual files in {self.parent.tle_data_folder}.")
            return saved_files

        def get_norad_id_name_map(self, limit: int | None = None) -> dict[int, str]:
            """
            Retrieves a mapping of NORAD IDs to object names from the satellite catalog.

            Parameters
            ----------
            limit : int | None, optional
                Limit the number of catalog entries to retrieve. If None, retrieves all.

            Returns
            -------
            dict[int, str]
                A dictionary mapping NORAD IDs (int) to object names (str).

            Raises
            ------
            RuntimeError
                If the SATCAT data fails to fetch.
            """
            query_parts = ['basicspacedata/query/class/satcat/orderby/NORAD_CAT_ID/format/json']
            if limit:
                query_parts.append(f"limit/{limit}")

            url = os.path.join(self.parent.spacetrack_url, "/".join(query_parts))
            json_name = 'norad_id_to_name.json'

            data = self.parent._get_json_and_save(url, json_name)

            if not data:
                raise RuntimeError("Failed to fetch SATCAT data.")

            map_data: dict[int, str] = {int(obj['NORAD_CAT_ID']): obj['OBJECT_NAME'] for obj in data if obj['OBJECT_NAME']}

            map_json_name = 'norad_id_to_name_map.json'
            map_path = os.path.join(self.parent.tle_data_folder, map_json_name)

            with open(map_path, "w") as f:
                json.dump(map_data, f, indent=4)
                print(f"Saved map to {map_path}.")

            return map_data

        def filter_tles_keep_latest_creation_from_json(self, full_filepath: str) -> any:
            """
            Filters duplicates. Requires full path.

            Parameters
            ----------
            full_filepath : str
                Full path to the JSON file containing TLEs.

            Returns
            -------
            any
                An iterable of filtered TLE dictionaries (latest creation date per epoch).
            """
            with open(full_filepath, "r") as f:
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
                        
            return filtered.values()

        def get_tles(self, json_dict: list[dict[str, any]] | dict[str, any]) -> dict[str, tuple[str, str]]:
            """
            Extracts TLE lines (Line 1 and Line 2) from JSON data.

            Parameters
            ----------
            json_dict : list[dict[str, any]] | dict[str, any]
                The JSON data containing TLE information.

            Returns
            -------
            dict[str, tuple[str, str]]
                A dictionary mapping NORAD IDs (str) to tuples of (Line 1, Line 2).
            """
            tle_dict: dict[str, tuple[str, str]] = {}
            if type(json_dict) is list:
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
            Extracts and converts Keplerian elements from TLE JSON data for Tudat usage.

            Parameters
            ----------
            json_dict : list | dict
                The JSON data containing TLE information.

            Returns
            -------
            tuple[float | None, float | None, float | None, float | None, float | None, float | None]
                A tuple of (a, e, i, omega, RAAN, true_anomaly) in SI units (meters, radians).
                Returns a tuple of Nones if input is empty.
            """
            if not json_dict:
                return (None,) * 6

            if json_dict and isinstance(json_dict[0], list):
                json_dict = json_dict[0]

            first_obj = json_dict[0]

            a = float(first_obj['SEMIMAJOR_AXIS'])*1e3
            e =float(first_obj['ECCENTRICITY'])
            i = float(first_obj['INCLINATION']) * math.pi/180
            omega = float(first_obj['ARG_OF_PERICENTER']) * math.pi/180
            raan = float(first_obj['RA_OF_ASC_NODE']) * math.pi/180
            mo = float(first_obj['MEAN_ANOMALY']) * math.pi/180

            true_anomaly = self.mean_to_true_anomaly(mo,e)
            return a,e,i,omega,raan,true_anomaly

        def tle_to_TleEphemeris_object(self, tle_line_1: str, tle_line_2: str) -> environment.TleEphemeris:
            """
            Converts TLE lines into a Tudat TleEphemeris object.

            Parameters
            ----------
            tle_line_1 : str
                The first line of the TLE.
            tle_line_2 : str
                The second line of the TLE.

            Returns
            -------
            environment.TleEphemeris
                The configured TleEphemeris object.
            """
            object_tle = environment.Tle(tle_line_1, tle_line_2)
            ephemeris_object = environment.TleEphemeris("Earth", "J2000", object_tle, False)
            return ephemeris_object

        def plot_earth(self, ax: any, radius: float = 6378, color: str = 'lightblue', alpha: float = 0.5, resolution: int = 50) -> None:
            """
            Plots a 3D sphere representing Earth on the given axes.

            Parameters
            ----------
            ax : any
                The matplotlib 3D axes to plot on.
            radius : float, optional
                Radius of the sphere. Defaults to 6378.
            color : str, optional
                Color of the sphere. Defaults to 'lightblue'.
            alpha : float, optional
                Transparency of the sphere. Defaults to 0.5.
            resolution : int, optional
                Mesh resolution of the sphere. Defaults to 50.
            """
            u = np.linspace(0, 2 * np.pi, resolution)
            v = np.linspace(0, np.pi, resolution)
            x = radius * np.outer(np.cos(u), np.sin(v))
            y = radius * np.outer(np.sin(u), np.sin(v))
            z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, rstride=1, cstride=1, color=color, alpha=alpha, edgecolor='none')

        def get_tle_reference_epoch(self, tle_line_1: str) -> datetime:
            """
            Parses the reference epoch from the first line of a TLE.

            Parameters
            ----------
            tle_line_1 : str
                The first line of the TLE.

            Returns
            -------
            datetime
                The reference epoch as a datetime object.
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
            Converts mean anomaly to true anomaly using Newton-Raphson iteration on Kepler's Equation.

            Parameters
            ----------
            mo : float
                Mean anomaly in radians.
            e : float
                Eccentricity.
            tol : float, optional
                Convergence tolerance. Defaults to 1e-8.
            max_iter : int, optional
                Maximum number of iterations. Defaults to 100.

            Returns
            -------
            float
                The true anomaly in radians.

            Raises
            ------
            RuntimeError
                If the iteration fails to converge.
            """
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

            v = 2 * np.arctan2(
                np.sqrt(1 + e) * np.sin(E / 2),
                np.sqrt(1 - e) * np.cos(E / 2)
            )
            return v
