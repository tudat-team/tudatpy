import json
import getpass
import os
import requests
from collections import defaultdict
from urllib.parse import urljoin, urlparse
from tudatpy.dynamics import environment
from datetime import datetime, timedelta
import math
from tudatpy.data import get_resource_path
from tudatpy.astro.element_conversion import mean_to_true_anomaly

class SpaceTrackQuery:
    """
    Handles queries to the Space-Track.org API for retrieving TLEs and other
    space data. Manages authentication, session persistence, and local caching
    to minimise API usage.

    For working with the retrieved OMM records, see :class:`OMMUtils`.
    """

    def __init__(
            self,
            username: str | None = None,
            password: str | None = None,
            spacetrack_url: str = "https://www.space-track.org",
            tle_data_folder: str = get_resource_path() + "/tle_data",
            timeout: int = 60,
    ) -> None:
        """
        Parameters
        ----------
        username : str | None, optional
            Space-Track.org username. Prompted interactively if not provided.
        password : str | None, optional
            Space-Track.org password. Prompted interactively if not provided.
        spacetrack_url : str, optional
            Base URL. Defaults to ``'https://www.space-track.org'``.
        tle_data_folder : str, optional
            Directory for locally cached TLE files.
        """
        if not username:
            username = input("space-track username: ")
        if not password:
            password = getpass.getpass("space-track password: ")

        self.username: str = username
        self.spacetrack_url: str = spacetrack_url
        self.tle_data_folder: str = tle_data_folder

        os.makedirs(self.tle_data_folder, exist_ok=True)

        self.session: requests.Session = requests.Session()
        self.timeout = timeout
        self._login(password)
        del password  # do not keep the plaintext password on the instance

    # ------------------------------------------------------------------
    # Context manager
    # ------------------------------------------------------------------

    def __enter__(self) -> "SpaceTrackQuery":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.logout()

    # ------------------------------------------------------------------
    # Authentication
    # ------------------------------------------------------------------

    def _login(self, password: str) -> None:t
        """
        Logs in to Space-Track.org.

        Parameters
        ----------
        
        password : str
            Plaintext password. Never stored on the instance.

        Raises
        ------
        requests.exceptions.RequestException
            On network errors or bad credentials.
        """
        print("Logging into Space-Track...")
        try:
            response = self.session.post(
                urljoin(self.spacetrack_url, "/ajaxauth/login"),
                json={"identity": self.username, "password": password}, timeout=self.timeout
            )
            response.raise_for_status()
            print("Login successful.")
        except requests.exceptions.RequestException as e:
            print(f"Login failed: {e}")
            raise

    def logout(self) -> None:
        """
        Logs out and closes the session. Safe to call multiple times.
        """
        try:
            self.session.get(urljoin(self.spacetrack_url, "/ajaxauth/logout"),
                             timeout = self.timeout)
            print("Logged out of Space-Track.")
        except requests.exceptions.RequestException as e:
            print(f"Logout request failed (session may already be expired): {e}")
        finally:
            self.session.close()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_url(self, *path_parts: str) -> str:
        """
        Joins path segments onto the base URL, safe on all platforms.

        Parameters
        ----------
        *path_parts : str
            URL path segments to join.

        Returns
        -------
        str
            The full URL.
        """
        path = "/".join(part.strip("/") for part in path_parts)
        return urljoin(self.spacetrack_url.rstrip("/") + "/", path)

    def _fetch_json(self, url: str) -> list:
        """
        GET a URL and return the body as a list.

        Parameters
        ----------
        url : str
            The URL to fetch.

        Returns
        -------
        list
            The JSON response body as a list.
        """
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        data = response.json()
        return data if isinstance(data, list) else [data]

    def _get_json_and_save(
            self, url: str, json_name: str, merge: bool = False
    ) -> list | None:
        """
        Fetches JSON from the API and persists it to the local cache.

        Parameters
        ----------
        url : str
            API endpoint.
        json_name : str
            Cache filename (relative to ``tle_data_folder``).
        merge : bool
            ``True`` → merge into existing file. ``False`` → overwrite.

        Returns
        -------
        list | None
            The freshly fetched records, or ``None`` on failure.
        """
        filepath = os.path.join(self.tle_data_folder, json_name)
        try:
            print(f"Fetching data for {json_name}...")
            new_data = self._fetch_json(url)

            if merge and os.path.exists(filepath):
                print(f"Merge requested. Updating {json_name}...")
                self._save_unique_sorted(filepath, new_data)
                with open(filepath, "r") as f:
                    final = json.load(f)
                return final.get("data", final) if isinstance(final, dict) else final
            else:
                print(f"Overwrite mode. Saving {json_name}...")
                if os.path.exists(filepath):
                    os.remove(filepath)
                self._save_unique_sorted(filepath, new_data)
                return new_data

        except requests.exceptions.RequestException as e:
            print(f"API Error: {e}")
            return None

    def _save_unique_sorted(self, filepath: str, new_data: list[dict]) -> None:
        """
        Merges ``new_data`` into the cache file and updates ``last_api_hit``.

        When two records share the same ``NORAD_CAT_ID + EPOCH`` key the one
        with the later ``CREATION_DATE`` wins.

        File structure on disk::

            {"last_api_hit": "<ISO datetime>", "data": [...]}

        Parameters
        ----------
        filepath : str
            Path to the JSON cache file.
        new_data : list[dict]
            New OMM records to merge.

        Returns
        -------
        None
        """
        existing: list[dict] = []
        if os.path.exists(filepath):
            try:
                with open(filepath, "r") as f:
                    content = json.load(f)
                existing = content.get("data", content) if isinstance(content, dict) else content
            except json.JSONDecodeError:
                existing = []

        unique: dict[str, dict] = {}
        for item in existing:
            unique[f"{item['NORAD_CAT_ID']}_{item['EPOCH']}"] = item
        for item in new_data:
            key = f"{item['NORAD_CAT_ID']}_{item['EPOCH']}"
            if key not in unique or item.get("CREATION_DATE", "") > unique[key].get("CREATION_DATE", ""):
                unique[key] = item

        sorted_data = sorted(unique.values(), key=lambda x: (x["EPOCH"], x["NORAD_CAT_ID"]))
        with open(filepath, "w") as f:
            json.dump({"last_api_hit": datetime.now().isoformat(), "data": sorted_data}, f, indent=4)

    # ------------------------------------------------------------------
    # Download methods
    # ------------------------------------------------------------------

    def get_tles_for_date_range(
            self,
            norad_id: int | str,
            start_date: str,
            end_date: str,
            override_last_api_hit: bool = False,
    ) -> list[dict] | None:
        """
        Retrieves TLEs for a single satellite over a date range.

        Respects a 1.5-hour API cooldown to avoid hammering the server when
        data is absent (e.g. the satellite had not yet launched).

        Parameters
        ----------
        norad_id : int | str
            NORAD catalogue ID.
        start_date : str
            Start of range, ``"YYYY-MM-DD"``.
        end_date : str
            End of range, ``"YYYY-MM-DD"``.
        override_last_api_hit : bool
            If ``True``, bypass the cooldown and force an API call.

        Returns
        -------
        list[dict] | None
            OMM records within the requested window.
        """
        norad_id = str(norad_id)
        filepath = os.path.join(self.tle_data_folder, f"tle_{norad_id}.json")

        req_start = datetime.strptime(start_date, "%Y-%m-%d")
        req_end = datetime.strptime(end_date, "%Y-%m-%d") + timedelta(days=1, microseconds=-1)

        local_data: list[dict] = []
        last_hit: datetime | None = None
        if os.path.exists(filepath):
            with open(filepath, "r") as f:
                content = json.load(f)
            if isinstance(content, dict):
                local_data = content.get("data", [])
                last_hit = datetime.fromisoformat(content["last_api_hit"]) if "last_api_hit" in content else None
            else:
                local_data = content  # legacy list format

        needs_fetch = True
        cooldown = timedelta(hours=1.5)

        if local_data:
            local_min = datetime.fromisoformat(local_data[0]["EPOCH"])
            local_max = datetime.fromisoformat(local_data[-1]["EPOCH"])
            if req_start >= local_min and req_end <= local_max:
                print("Local cache fully covers request. Skipping API.")
                needs_fetch = False

        if needs_fetch and last_hit:
            time_since = datetime.now() - last_hit
            if not override_last_api_hit and time_since < cooldown:
                print(f"Cache is recent ({time_since} ago). Skipping API despite gaps.")
                needs_fetch = False

        if needs_fetch:
            print(f"Fetching {norad_id} from API ({start_date} -- {end_date})...")
            url = self._build_url(
                f"basicspacedata/query/class/gp_history/NORAD_CAT_ID/{norad_id}"
                f"/EPOCH/{start_date}--{end_date}/orderby/EPOCH asc/format/json"
            )
            fetched = self._fetch_json(url)
            self._save_unique_sorted(filepath, fetched or [])
            with open(filepath, "r") as f:
                local_data = json.load(f).get("data", [])

        return [
            omm for omm in local_data
            if req_start <= datetime.fromisoformat(omm["EPOCH"]) <= req_end
        ]

    def latest_on_orbit(
            self, update_existing: bool = False, filename: str | None = None
    ) -> list | None:
        """
        Retrieves the newest propagable element set for all on-orbit payloads.

        Parameters
        ----------
        update_existing : bool
            ``True`` → merge into existing file. ``False`` → overwrite.
        filename : str | None
            Optional filename override. Defaults to ``'latest_on_orbit.json'``.

        Returns
        -------
        list | None
            List of OMM records or None if the request failed.
        """
        url = self._build_url(
            "basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD"
            "/decay_date/null-val/epoch/>now-30/orderby/norad_cat_id/format/json"
        )
        return self._get_json_and_save(url, filename or "latest_on_orbit.json", merge=update_existing)

    def descending_epoch(
            self,
            N: int | None = None,
            update_existing: bool = False,
            filename: str | None = None,
    ) -> list | None:
        """
        Retrieves GP data ordered by epoch (most recent first).

        Parameters
        ----------
        N : int | None
            Optional result limit.
        update_existing : bool
            ``True`` → merge. ``False`` → overwrite.
        filename : str | None
            Optional filename override.

        Returns
        -------
        list | None
            List of OMM records or None if the request failed.
        """
        if filename:
            json_name = filename
        elif N is not None:
            json_name = f"gp_descending_limit_{N}.json"
        else:
            json_name = "gp_descending.json"

        parts = ["basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/orderby/epoch desc"]
        if N is not None:
            parts.append(f"limit/{N}")
        parts.append("format/json")
        return self._get_json_and_save(self._build_url("/".join(parts)), json_name, merge=update_existing)

    def get_tles_by_norad_ids(
            self,
            norad_ids: int | list[int],
            history: bool = False,
            orderby: str = "epoch desc",
            limit_per_object: int = 1,
            update_existing: bool = False,
            filename: str | None = None,
    ) -> list | None:
        """
        Retrieves TLEs for one or more specific NORAD IDs.

        Parameters
        ----------
        norad_ids : int | list[int]
            One or more NORAD IDs.
        history : bool
            If ``True``, queries ``gp_history`` for historical records.
        orderby : str
            API ordering string.
        limit_per_object : int
            Maximum records per object.
        update_existing : bool
            ``True`` → merge. ``False`` → overwrite.
        filename : str | None
            Force a specific cache filename.

        Returns
        -------
        list | None
            List of OMM records or None if the request failed.
        """
        if not isinstance(norad_ids, (list, tuple, set)):
            norad_ids = [norad_ids]

        id_string = ",".join(map(str, norad_ids))
        tle_class = "gp"
        use_orderby = True

        if history or limit_per_object > 1:
            tle_class = "gp_history"
            if len(norad_ids) > 1:
                print("Note: Removing 'orderby' for batch history query.")
                use_orderby = False

        if filename:
            json_name = filename
        elif len(norad_ids) == 1:
            json_name = f"{tle_class}_{norad_ids[0]}_limit_{limit_per_object}.json"
        else:
            json_name = f"{tle_class}_batch_{len(norad_ids)}_ids_limit_{limit_per_object}.json"

        parts = ["basicspacedata/query", f"class/{tle_class}", f"NORAD_CAT_ID/{id_string}"]
        if use_orderby:
            parts.append(f"orderby/{orderby}")
        parts.extend([f"limit/{limit_per_object}", "format/json"])

        json_data = self._get_json_and_save(self._build_url("/".join(parts)), json_name, merge=update_existing)

        return json_data

    def filtered_by_oe_dict(
            self,
            filter_oe_dict: dict[str, tuple[float | None, float | None]],
            limit: int = 100,
            output_file: str = "filtered_results.json",
            update_existing: bool = False,
    ) -> list | None:
        """
        Retrieves payloads filtered by orbital element bounds.

        Parameters
        ----------
        filter_oe_dict : dict
            Mapping of orbital element names to ``(min, max)`` bound tuples.
            Either bound may be ``None`` for a one-sided filter.
        limit : int
            Maximum number of results.
        output_file : str
            Cache filename.
        update_existing : bool
            ``True`` → merge. ``False`` → overwrite.

        Returns
        -------
        list | None
            List of OMM records or None if the request failed.
        """
        parts = ["basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD"]
        for oe, bounds in filter_oe_dict.items():
            min_val = bounds[0] if len(bounds) > 0 else None
            max_val = bounds[1] if len(bounds) > 1 else None
            if min_val is not None and max_val is not None:
                parts.extend([oe.upper(), f"{min_val}--{max_val}"])
            elif min_val is not None:
                parts.extend([oe.upper(), f">{min_val}"])
            elif max_val is not None:
                parts.extend([oe.upper(), f"<{max_val}"])
        parts.append(f"orderby/epoch desc/limit/{limit}/format/json")
        return self._get_json_and_save(self._build_url("/".join(parts)), output_file, merge=update_existing)

    def query_from_query_builder_url(
            self,
            query: str,
            output_file: str = "custom_query.json",
            update_existing: bool = False,
    ) -> list | None:
        """
        Executes a user-provided Space-Track query URL or query path.

        The input may be:
        - a full URL copied from the Space-Track 'Query Builder'
        - only the query path starting from 'basicspacedata/query/...'

        Parameters
        ----------
        query : str
            Full Space-Track API URL or query path.
        output_file : str
            Local cache filename.
        update_existing : bool
            True → merge with existing cache. False → overwrite.

        Returns
        -------
        list | None
            JSON records returned by the API.
        """

        query = query.strip()

        if query.startswith("http"): # Full URL provided
            parsed = urlparse(query)

            if "space-track.org" not in parsed.netloc:
                raise ValueError("Only space-track.org URLs are allowed.")

            # Extract only the path part
            path = parsed.path.lstrip("/")

        else:
            path = query.lstrip("/") # Only path provided

        # Ensure correct prefix
        if not path.startswith("basicspacedata/query"):
            raise ValueError(
                "Query must start with 'basicspacedata/query/...'"
            )

        # Ensure JSON format
        if "format/json" not in path:
            if not path.endswith("/"):
                path += "/"
            path += "format/json"

        url = self._build_url(path)

        return self._get_json_and_save(
            url,
            output_file,
            merge=update_existing
        )

class OMMUtils:
    """
    Stateless utility class for working with OMM (Orbit Mean-Elements Message)
    records in dictionary form.

    This class has no dependency on Space-Track.org or any particular data
    source — it operates on plain OMM dicts regardless of where they came from.
    All methods are static; instantiation is not required.
    """

    @staticmethod
    def save_batch_to_individual_files(
            json_data: list[dict],
            output_folder: str,
    ) -> list[str] | None:
        """
        Splits a batch list of OMM records into one JSON file per NORAD ID.

        Parameters
        ----------
        json_data : list[dict]
            Batch OMM records, each containing a ``NORAD_CAT_ID`` field.
        output_folder : str
            Directory in which to write the per-satellite files.

        Returns
        -------
        list[str] | None
            List of created filenames (basenames only), or None if no data.
        """
        if not json_data:
            print("No data to split.")
            return None

        grouped: dict[str, list] = defaultdict(list)
        for record in json_data:
            grouped[record["NORAD_CAT_ID"]].append(record)

        saved_files: list[str] = []
        for norad_id, records in grouped.items():
            filename = f"{norad_id}_from_batch.json"
            filepath = os.path.join(output_folder, filename)
            with open(filepath, "w") as f:
                json.dump(records, f, indent=4)
            saved_files.append(filename)

        print(f"Split batch data into {len(saved_files)} individual files in {output_folder}.")
        return saved_files

    @staticmethod
    def get_tles(json_data: list[dict] | dict) -> defaultdict[str, list[tuple[str, str]]]:
        """
        Extracts TLE line pairs from a list of OMM records.

        Parameters
        ----------
        json_data : list[dict] | dict
            One or more OMM records containing ``TLE_LINE1`` and ``TLE_LINE2``.

        Returns
        -------
        dict[str, tuple[str, str]]
            Mapping of NORAD ID (str) → (TLE line 1, TLE line 2).
        """
        if isinstance(json_data, dict):
            json_data = [json_data]

        final_dict = defaultdict(list)

        for entry in json_data:
            final_dict[entry["NORAD_CAT_ID"]].append(
                (entry["TLE_LINE1"], entry["TLE_LINE2"])
            )

        return final_dict



    @staticmethod
    def get_tudat_keplerian_element_set(
            omm: dict,
    ) -> tuple[float, float, float, float, float, float]:
        """
        Extracts and converts Keplerian elements from a single OMM record into
        SI units suitable for Tudat.

        Parameters
        ----------
        omm : dict
            A single OMM record (km and degrees, as returned by Space-Track or
            any compliant source).

        Returns
        -------
        tuple[float, float, float, float, float, float]
            ``(a, e, i, omega, RAAN, true_anomaly)`` where distances are in
            metres and angles are in radians.

        Raises
        ------
        TypeError
            If ``omm`` is a list rather than a single dict.
        ValueError
            If ``omm`` is empty or falsy.
        """
        if not omm:
            raise ValueError("No OMM record provided.")
        if isinstance(omm, list):
            raise TypeError(
                "omm must be a single dictionary, not a list. "
                "Pass one record at a time."
            )

        a     = float(omm["SEMIMAJOR_AXIS"])    * 1e3
        e     = float(omm["ECCENTRICITY"])
        i     = float(omm["INCLINATION"])        * math.pi / 180
        omega = float(omm["ARG_OF_PERICENTER"]) * math.pi / 180
        raan  = float(omm["RA_OF_ASC_NODE"])    * math.pi / 180
        mo    = float(omm["MEAN_ANOMALY"])       * math.pi / 180

        return a, e, i, omega, raan, mean_to_true_anomaly(e, mo)

    @staticmethod
    def tle_to_TleEphemeris_object(
            tle_line_1: str, tle_line_2: str
    ) -> environment.Tle:
        """
        Converts a TLE line pair into a Tudat ``TleEphemeris`` object.

        Parameters
        ----------
        tle_line_1 : str
            First line of the TLE.
        tle_line_2 : str
            Second line of the TLE.

        Returns
        -------
        environment.Tle
            Configured TleEphemeris object.
        """
        return environment.TleEphemeris(
            "Earth", "J2000", environment.Tle(tle_line_1, tle_line_2), False
        )

    @staticmethod
    def tle_to_Tle_object(
            tle_line_1: str, tle_line_2: str
    ) -> environment.TleEphemeris:
        """
        Converts a TLE line pair into a Tudat ``TleEphemeris`` object.

        Parameters
        ----------
        tle_line_1 : str
            First line of the TLE.
        tle_line_2 : str
            Second line of the TLE.

        Returns
        -------
        environment.Tle
            Tle object.
        """
        return environment.Tle(tle_line_1, tle_line_2)

    @staticmethod
    def clean_file(filepath: str) -> None:
        """
        Removes duplicate TLE entries from a cache file in-place.

        Supports both legacy (plain list) and current (dict with metadata)
        formats. When two records share the same EPOCH, the one with the later
        ``CREATION_DATE`` is kept.

        Parameters
        ----------
        filepath : str
            Absolute path to the cache file.

        Returns
        -------
        None
        """
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            return

        with open(filepath, "r") as f:
            content = json.load(f)

        if isinstance(content, list):
            data, metadata, is_dict = content, {}, False
        elif isinstance(content, dict) and "data" in content:
            data = content["data"]
            metadata = {k: v for k, v in content.items() if k != "data"}
            is_dict = True
        else:
            print(f"Skipping {filepath}: unknown format.")
            return

        initial_count = len(data)
        unique: dict[str, dict] = {}
        for entry in data:
            epoch = entry.get("EPOCH")
            norad_id = entry.get("NORAD_CAT_ID")
            if not epoch or not norad_id:
                continue

            composite_key = f"{norad_id}_{epoch}"

            if composite_key not in unique or entry.get("CREATION_DATE", "") > unique[composite_key].get("CREATION_DATE", ""):
                unique[composite_key] = entry

        cleaned = sorted(unique.values(), key=lambda x: (x["EPOCH"], x["NORAD_CAT_ID"]))
        output: dict | list = ({**metadata, "data": cleaned} if is_dict else cleaned)

        with open(filepath, "w") as f:
            json.dump(output, f, indent=4)

        print(f"Cleaned {filepath}. Removed {initial_count - len(cleaned)} duplicates.")
