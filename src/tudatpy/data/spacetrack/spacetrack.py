import json
import getpass
import os
import requests
from collections import defaultdict
import numpy as np
from mypyc.ir.rtypes import none_rprimitive
from tudatpy.dynamics import environment, environment_setup, propagation_setup
from datetime import datetime, timedelta, date
import math
from tudatpy.data import get_resource_path

class SpaceTrackQuery:
    """
    A class to handle queries to Space-Track.org API for retrieving TLEs and other space data.
    It manages authentication, session persistence, and local caching of TLE files to minimize API usage.
    """

    def __init__(self, username: str | None = None, password: str | None = None, spacetrack_url = 'https://www.space-track.org', tle_data_folder: str = get_resource_path() + "/tle_data") -> None:
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

        self.spacetrack_url: str = spacetrack_url
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

    def _get_json_and_save(self, url: str, json_name: str, merge: bool = False) -> dict | list | None:
        """
        Retrieves JSON data.

        Parameters
        ----------
        url : str
            API Endpoint.
        json_name : str
            Filename.
        merge : bool
            If True, merges new data with existing file (Append Mode).
            If False, overwrites the file with new data (Replace Mode).
        """
        filepath = os.path.join(self.tle_data_folder, json_name)

        # 1. Fetch from API
        try:
            print(f"Fetching data for {json_name}...")
            response = self.session.get(url)
            response.raise_for_status()
            new_data = response.json()

            # 2. Save Logic
            if merge and os.path.exists(filepath):
                print(f"Merge requested. Updating {json_name}...")
                # Use our robust merger to combine old + new
                self._save_unique_sorted(filepath, new_data)

                # Reload the full merged set to return it
                with open(filepath, 'r') as f:
                    final_data = json.load(f)
                    # Handle new dict structure if present
                    if isinstance(final_data, dict) and 'data' in final_data:
                        return final_data['data']
                    return final_data

            else:
                # Overwrite Mode (Default)
                print(f"Overwrite requested (or file missing). Saving {json_name}...")

                # _save_unique_sorted reads the file if it exists.
                # To force overwrite, we should just delete the file first if it exists.
                if os.path.exists(filepath):
                    os.remove(filepath)

                self._save_unique_sorted(filepath, new_data)
                return new_data

        except requests.exceptions.RequestException as e:
            print(f"API Error: {e}")
            return None
    def get_tles_for_date_range(self, norad_id: int | str, start_date: str, end_date: str, override_last_api_hit: bool = False) -> list[dict] | None:
        """
        Retrieves TLEs. Respects a 2-hour API cooldown to prevent spamming queries
        for missing data (e.g., if start_date is before the satellite was launched).
        """
        norad_id = str(norad_id)
        filename = f"tle_{norad_id}.json"
        filepath = os.path.join(self.tle_data_folder, filename)

        # 1. Parse dates
        req_start = datetime.strptime(start_date, "%Y-%m-%d")
        req_end = datetime.strptime(end_date, "%Y-%m-%d") + timedelta(days=1, microseconds=-1)

        # 2. Load Local Data & Metadata
        local_data = []
        last_hit = None

        if os.path.exists(filepath):
            try:
                with open(filepath, 'r') as f:
                    content = json.load(f)

                    if isinstance(content, dict):
                        local_data = content.get('data', [])
                        if 'last_api_hit' in content:
                            last_hit = datetime.fromisoformat(content['last_api_hit'])
                    elif isinstance(content, list):
                        # Legacy format support
                        local_data = content
            except json.JSONDecodeError:
                pass

        # 3. Check Cooldown Logic
        needs_fetch = True

        # Check 1: Do we have full coverage?
        if local_data:
            local_min = datetime.fromisoformat(local_data[0]['EPOCH'])
            local_max = datetime.fromisoformat(local_data[-1]['EPOCH'])

            if req_start >= local_min and req_end <= local_max:
                print(f"Local cache fully covers request. Skipping API.")
                needs_fetch = False

        # Check 2: Are we in the cooldown period?
        if needs_fetch and last_hit:
            time_since_last_hit = datetime.now() - last_hit
            cooldown_period = timedelta(hours=1.5)

            if not override_last_api_hit and time_since_last_hit < cooldown_period:
                print(f"Cache is recent ({time_since_last_hit} ago). Skipping API despite gaps.")
                needs_fetch = False
            elif override_last_api_hit and time_since_last_hit < cooldown_period:
                needs_fetch = True
        # 4. Fetch from API (Only if needed AND allowed)
        if needs_fetch:
            print(f"Fetching {norad_id} from API ({start_date} -- {end_date})...")
            url = os.path.join(
                self.spacetrack_url,
                f"basicspacedata/query/class/gp/NORAD_CAT_ID/{norad_id}/EPOCH/{start_date}--{end_date}/orderby/EPOCH asc/format/json"
            )

            # Even if API returns empty list, we MUST save to update the timestamp!
            fetched_data = self._fetch_json(url)

            # This updates 'last_api_hit' and merges any new data
            self._save_unique_sorted(filepath, fetched_data if fetched_data else [])

            # Reload fresh data
            with open(filepath, 'r') as f:
                content = json.load(f)
                local_data = content.get('data', [])

        # 5. Filter and Return
        filtered_omm = []
        for omm in local_data:
            tle_epoch = datetime.fromisoformat(omm['EPOCH'])
            if req_start <= tle_epoch <= req_end:
                filtered_omm.append(omm)

        return filtered_omm
    def _fetch_json(self, url: str) -> list | None:
        """Fetches and returns list."""
        response = self.session.get(url)
        response.raise_for_status()
        data = response.json()
        return data if isinstance(data, list) else [data]

    def _save_unique_sorted(self, filepath: str, new_data: list[dict]):
        """
        Merges new data into the file and updates the 'last_api_hit' timestamp.
        File Structure: {"last_api_hit": "ISO_STR", "data": [ ... ]}
        """
        # 1. Load existing content
        existing_data = []
        if os.path.exists(filepath):
            try:
                with open(filepath, 'r') as f:
                    content = json.load(f)

                    # Handle migration from old list-only format to new dict format
                    if isinstance(content, list):
                        existing_data = content
                    elif isinstance(content, dict):
                        existing_data = content.get('data', [])
            except json.JSONDecodeError:
                existing_data = []

        # This prevents Sat A from overwriting Sat B if they share an Epoch
        unique_map = {}
        for item in existing_data:
            key = f"{item['NORAD_CAT_ID']}_{item['EPOCH']}"
            unique_map[key] = item

        for new_item in new_data:
            key = f"{new_item['NORAD_CAT_ID']}_{new_item['EPOCH']}"

            if key not in unique_map:
                unique_map[key] = new_item
            else:
                # Revision check (check creation date)
                if new_item.get('CREATION_DATE', '') > unique_map[key].get('CREATION_DATE', ''):
                    unique_map[key] = new_item

        # Sort by EPOCH (and ID for cleanliness)
        sorted_data = sorted(unique_map.values(), key=lambda x: (x['EPOCH'], x['NORAD_CAT_ID']))

        file_content = {
            "last_api_hit": datetime.now().isoformat(),
            "data": sorted_data
        }

        with open(filepath, 'w') as f:
            json.dump(file_content, f, indent=4)

    def clean_tle_file(self, filepath: str) -> None:
        """
        Removes duplicate TLE entries from a JSON file.
        Supports both legacy (list) and new (dict with metadata) formats.
        """
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            return

        try:
            with open(filepath, 'r') as f:
                content = json.load(f)

            # 1. Detect Format and Extract Data
            is_dict_format = False
            metadata = {}

            if isinstance(content, list):
                data = content
            elif isinstance(content, dict) and 'data' in content:
                is_dict_format = True
                data = content['data']
                # Preserve metadata (like last_api_hit)
                metadata = {k: v for k, v in content.items() if k != 'data'}
            else:
                print(f"Skipping {filepath}: Unknown format.")
                return

            # 2. Deduplicate
            initial_count = len(data)
            unique_map = {}

            for entry in data:
                epoch = entry.get('EPOCH')
                if not epoch: continue

                if epoch not in unique_map:
                    unique_map[epoch] = entry
                else:
                    # Conflict resolution: keep newest CREATION_DATE
                    existing = unique_map[epoch]
                    if entry.get('CREATION_DATE', '') > existing.get('CREATION_DATE', ''):
                        unique_map[epoch] = entry

            cleaned_data = sorted(unique_map.values(), key=lambda x: x['EPOCH'])
            final_count = len(cleaned_data)

            # 3. Save back in the same format it was found
            if is_dict_format:
                output = metadata
                output['data'] = cleaned_data
            else:
                output = cleaned_data

            with open(filepath, 'w') as f:
                json.dump(output, f, indent=4)

            print(f"Successfully cleaned {filepath}. Removed {initial_count - final_count} duplicates.")

        except json.JSONDecodeError:
            print(f"Error: {filepath} contains invalid JSON.")

    #####################################
    # --- AVAILABLE DOWNLOAD METHODS ---#
    #####################################
    def latest_on_orbit(self, update_existing: bool = False, filename: str | None = None) -> dict | list | None:
        """
        Retrieves the newest propagable element set for all on-orbit payloads.

        Parameters
        ----------
        update_existing : bool
            If True, merges the new snapshot into the target file.
            If False (default), overwrites the target file with the current snapshot.
        filename : str | None
            Optional override for the filename.
            If None, defaults to 'latest_on_orbit.json'.
        """
        # 1. Determine Filename
        if filename:
            json_name = filename
        else:
            json_name = "latest_on_orbit.json"

        # Note: We use 'epoch/>now-30' to get data from the last 30 days.
        url = os.path.join(
            self.spacetrack_url,
            'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/decay_date/null-val/epoch/>now-30/orderby/norad_cat_id/format/json'
        )

        return self._get_json_and_save(url, json_name, merge=update_existing)

    def descending_epoch(self, N: int | None = None, update_existing: bool = False, filename: str | None = None) -> dict | list | None:
        """
        Retrieves GP data ordered by epoch.

        Parameters
        ----------
        N : int | None
            Limit number of results.
        update_existing : bool
            If True, merges results into the local file.
        filename : str | None
            Optional override for the filename.
            Useful if you want to update a 'master' file with a small N query.
        """
        # 1. Determine Filename
        if filename:
            json_name = filename
        elif N is not None:
            json_name = f"gp_descending_limit_{N}.json"
        else:
            json_name = "gp_descending.json"

        url_parts = ['basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/orderby/epoch desc']

        if N is not None:
            url_parts.append(f'limit/{N}')

        url_parts.append('format/json')
        url = os.path.join(self.spacetrack_url, "/".join(url_parts))

        return self._get_json_and_save(url, json_name, merge=update_existing)

    def get_tles_by_norad_ids(
            self, norad_ids: int | list[int], history: bool = False,
            orderby: str = 'epoch desc', limit_per_object: int = 1,
            update_existing: bool = False, filename: str | None = None
    ) -> dict | list | None:
        """
        Retrieves TLEs for specific IDs.

        Parameters
        ----------
        norad_ids : int | list[int]
            A single NORAD ID or a list of NORAD IDs to query.
        history : bool
            If True, queries 'gp_history'.
        orderby : str
            Ordering criteria.
        limit_per_object : int
            Limit TLEs per object.
        update_existing : bool
            If True, merges results into the target file.
        filename : str | None
            Force a specific filename. Essential if building a custom list
            (e.g., 'my_favorites.json') from multiple queries.
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

        # 1. Determine Filename
        if filename:
            json_name = filename
        elif len(norad_ids) == 1:
            single_id = norad_ids[0]
            json_name = f"{tle_class}_{single_id}_limit_{limit_per_object}.json"
        else:
            json_name = f"{tle_class}_batch_{len(norad_ids)}_ids_limit_{limit_per_object}.json"

        # Build URL
        url_parts = ['basicspacedata/query', f'class/{tle_class}', f'NORAD_CAT_ID/{id_string}']
        if use_orderby:
            url_parts.append(f'orderby/{orderby}')
        url_parts.extend([f'limit/{limit_per_object}', 'format/json'])

        url = os.path.join(self.spacetrack_url, "/".join(url_parts))

        # 2. Fetch using the new Merge logic
        json_data = self._get_json_and_save(url, json_name, merge=update_existing)

        # 3. Post-Process (Filter Duplicates)
        # Since _get_json_and_save now handles deduplication via _save_unique_sorted
        # internally if merge=True, we strictly only need this extra step
        # if you want to ensure the file is perfectly clean based on *Creation Date* specifically
        # or if _save_unique_sorted wasn't called (e.g. merge=False).
        # However, keeping it is safe as it re-verifies the file.
        if json_data and (history or limit_per_object > 1):
            filepath = os.path.join(self.tle_data_folder, json_name)

            # Note: We read from the file to ensure we are filtering the *merged* result
            filtered_values = list(self.OMMUtils.filter_tles_keep_latest_creation_from_json(filepath))

            # Since filter_tles_keep_latest... returns a list, we need to handle
            # the dictionary structure if we want to keep metadata,
            # OR just save the list back if we don't care about 'last_api_hit' here.
            # For simplicity in this specific helper method, we just dump the list.
            with open(filepath, "w") as f:
                json.dump(filtered_values, f, indent=4)
            print(f"Filtered local file content: {json_name}")

        return json_data

    def filtered_by_oe_dict(
            self, filter_oe_dict: dict[str, tuple[float | None, float | None]],
            limit: int = 100, output_file: str = 'filtered_results.json',
            update_existing: bool = False
    ) -> dict | list | None:
        """
        Retrieves TLEs by orbital elements.

        Parameters
        ----------
        filter_oe_dict : dict
            Dictionary of Orbital Elements bounds.
        limit : int
            Limit results.
        output_file : str
            Filename to save/merge results.
        update_existing : bool
            If True, appends results to 'output_file'.
            If False, overwrites 'output_file'.
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

        # Pass the merge flag directly to the helper
        return self._get_json_and_save(url, output_file, merge=update_existing)

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

            # NOTE: This function no longer saves to "temp_filtered_tles.json".
            # It simply returns the filtered values for the parent to save properly.
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