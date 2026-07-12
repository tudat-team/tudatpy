"""Rosetta (RO) file discovery and download."""

import requests
from bs4 import BeautifulSoup
import os
import re
import time
from ._download import _REQUEST_TIMEOUT


class RoMixin:
    def get_ro_files(
        self,
        local_folder,
        start_date,
        end_date,
        radio_observation_type=None,
        skip_kernel_downloads=False,
    ):
        """
        Description:
        This function downloads and organizes various types of SPICE kernels and
        radio science files for the ROSETTA mission.
        It supports downloading kernel files related to radio science, clock data, frame data,
        SPK (spacecraft position) data, CK (orientation) data, and tropospheric/ionospheric correction files.
        The function interacts with remote FTP servers to retrieve the data and organizes them into the specified local folder.
        The data is filtered by a given date range, and the function returns the downloaded files for further use.

        Inputs:
            - local_folder (`str`): The local directory where the downloaded files will be saved.
            - start_date (`datetime`): The start date for downloading data. This will filter the data to include only those within the date range.
            - end_date (`datetime`): The end date for downloading data. This will filter the data to include only those within the date range.
            - radio_observation_type (`str`): The type of radio science files to download (e.g. commissioning, checkout, solar conjunction, Lutetia, Global Gravity etc...)

        Outputs:
            (`dict`, `dict`, `dict`): A tuple containing:

            - `kernel_files_to_load` (`dict`): A dictionary where the keys are kernel types (e.g., 'ck', 'spk', 'fk', 'sclk') and values are lists of paths to the successfully downloaded and loaded kernel files.
            - `radio_science_files_to_load` (`dict`): A dictionary where keys are categories of radio science data (e.g., 'ifms_dp2', 'dsn_dps') and values are lists of paths to the successfully downloaded radio science files.
            - `ancillary_files_to_load` (`dict`): A dictionary where keys are categories of ancillary data (e.g., 'ion', 'tropospheric') and values are lists of paths to the successfully downloaded ancillary files, such as tropospheric and ionospheric corrections.

        """

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}  # empty for now

        input_mission = "ro"

        # Fetch radio science URLs once and reuse for both tropospheric and radio science loops
        print("================================================================")
        print(f"Discovering {input_mission.upper()} Radio Science Archive URLs:")
        cached_radio_science_urls = self.get_url_ro_radio_science_files(
            start_date, end_date, radio_observation_type
        )
        print(f"Found {len(cached_radio_science_urls)} Radio Science archive URL(s).")

        # Tropospheric corrections
        print("================================================================")
        print(f"Download {input_mission.upper()} Tropospheric and Ionospheric Corrections Files")
        for url_tropo_file_new in cached_radio_science_urls:
            for folder_type in [
                "CALIB/CLOSED_LOOP/IFMS/MET/",
                "CALIB/CLOSED_LOOP/DSN/MET/",
            ]:
                url_tropo_file = url_tropo_file_new + folder_type
                response = requests.get(url_tropo_file, timeout=_REQUEST_TIMEOUT)
                if response.status_code == 200:
                    html = response.text
                    # Parse the HTML with BeautifulSoup
                    soup = BeautifulSoup(html, "html.parser")
                    # Extract file links and their names
                    wanted_tropo_files = []
                    for link in soup.find_all("a"):
                        href = link.get("href")
                        if href.endswith(".TAB") or href.endswith(".AUX"):
                            wanted_tropo_files.append(href.split("/")[-1])
                    tropo_files_to_load = self.get_kernels(
                        input_mission=input_mission,
                        url=url_tropo_file,
                        wanted_files_patterns=["*L3L1B*.TAB"],
                        custom_output=local_folder,
                    )
                else:
                    continue

                if tropo_files_to_load:
                    key = folder_type.split("/")[-2]
                    if key not in self.ancillary_files_to_load:
                        self.ancillary_files_to_load[key] = []
                    self.ancillary_files_to_load[key].extend(tropo_files_to_load)

                else:
                    print("No tropospheric or ionospheric files to download this time.")

        print("================================================================")
        print(f"Download {input_mission.upper()} Radio Science Kernels:")

        # Build all DP2 URLs
        radio_science_dp2_urls = []
        for url_radio_science_file_new in cached_radio_science_urls:
            for closed_loop_type in ["IFMS/DP2/"]:
                url_radio_science_file = (
                    url_radio_science_file_new + "DATA/LEVEL02/CLOSED_LOOP/" + closed_loop_type
                )
                radio_science_dp2_urls.append((url_radio_science_file, closed_loop_type))

        print(f"Scanning {len(radio_science_dp2_urls)} archive URL(s) for IFMS/DP2 files...")

        # Process all URLs quietly and collect results
        all_skipped = []
        all_downloaded = []
        for url_radio_science_file, closed_loop_type in radio_science_dp2_urls:
            try:
                files = self.dynamic_download_url_files_single_time(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=url_radio_science_file,
                    verbose=False,
                )
                key = f"{closed_loop_type.split('/')[0]}_{closed_loop_type.split('/')[1]}"
                if key not in self.radio_science_files_to_load:
                    self.radio_science_files_to_load[key] = []
                self.radio_science_files_to_load[key].extend(files)
                all_skipped.extend(self._last_skipped_files)
                all_downloaded.extend(self._last_downloaded_files)
            except Exception as e:
                print(f"Error downloading radio science files from {closed_loop_type}: {e}")
                continue

        # Print summary
        if all_skipped:
            formatted = self._format_existing_files(all_skipped)
            print(
                f"\nThe following Radio Science files already exist and will not be downloaded:\n\n{formatted}"
            )
        if all_downloaded:
            dest_folder = os.path.join(local_folder, "dp2")
            total = len(all_downloaded)
            print(f"\nDownloaded {total} new Radio Science file(s) to {dest_folder}/:\n")
            for i, f in enumerate(all_downloaded, 1):
                print(f"  [{i}/{total}] {os.path.basename(f)}")
            print()
        if not all_skipped and not all_downloaded:
            print("No Radio Science files found for the given date range.")

        if not skip_kernel_downloads:
            # Clock files
            print("================================================================")
            print(f"Download {input_mission.upper()} Clock Kernels:")
            url_clock_files = "https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/sclk/"

            wanted_clock_files = self.get_latest_clock_kernel_name(input_mission)

            clock_files_to_load = self.get_kernels(
                input_mission=input_mission,
                url=url_clock_files,
                wanted_files=wanted_clock_files,
                custom_output=local_folder,
            )

            if clock_files_to_load:
                self.kernel_files_to_load["sclk"] = clock_files_to_load
            else:
                print("No sclk files to download this time.")

            print("================================================================")
            print(f"Download {input_mission.upper()} Frame Kernels:")
            url_frame_files = "https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/fk/"
            wanted_frame_files = self._get_latest_versioned_files(
                url_frame_files,
                [
                    r"ROS_V(\d+)\.TF",
                    r"ROS_DSK_SURFACES_V(\d+)\.TF",
                    r"ROS_CGS_AUX_V(\d+)\.TF",
                    r"ROS_LUTETIA_RSOC_V(\d+)\.TF",
                ],
            )
            frame_files_to_load = self.get_kernels(
                input_mission=input_mission,
                url=url_frame_files,
                wanted_files=wanted_frame_files,
                custom_output=local_folder,
            )

            if frame_files_to_load:
                self.kernel_files_to_load["fk"] = frame_files_to_load
            else:
                print("No fk files to download this time.")

            # Spk files
            print("================================================================")
            print(f"Download {input_mission.upper()} SPK Kernels:")
            url_spk_files = ["https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/spk/"]

            spk_files_to_load = []
            for url_spk_file in url_spk_files:
                spk_files_to_load.extend(
                    self.dynamic_download_url_files_time_interval(
                        input_mission,
                        local_path=local_folder,
                        start_date=start_date,
                        end_date=end_date,
                        url=url_spk_file,
                    )
                )

            if spk_files_to_load:
                self.kernel_files_to_load["spk"] = spk_files_to_load
            else:
                print("No spk files to download this time.")

            # Orientation files
            print("================================================================")
            print(f"Download {input_mission.upper()} CK Kernels:")
            url_ck_files = ["https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/ck/"]

            ck_files_to_load = []
            for url_ck_file in url_ck_files:
                ck_files_to_load.extend(
                    self.dynamic_download_url_files_time_interval(
                        input_mission,
                        local_path=local_folder,
                        start_date=start_date,
                        end_date=end_date,
                        url=url_ck_file,
                    )
                )

            if ck_files_to_load:
                self.kernel_files_to_load["ck"] = ck_files_to_load
            else:
                print("No ck files to download this time.")

        print("----------------------------------------------------------------")
        print(
            "All requested, relevant and previously non-existing RO files have been now downloaded. Enjoy!"
        )

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    #########################################################################################################

    def get_url_ro_radio_science_files(
        self, start_date_ro, end_date_ro, radio_observation_type=None
    ):

        url = "https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/RSI/RO-C-RSI-1-2-3-EXT3-1881-V1.0/AAREADME.TXT"
        radio_science_base_url = (
            "https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/RSI/"
        )

        mapping_dict = self.get_ro_rsi_volume_ID_mapping(url)
        mapping_dict = self.add_ro_mission_phase_designation(mapping_dict)

        if radio_observation_type:
            mapping_dict = self.filter_mapping_dict_by_radio_observation_type(
                mapping_dict,
                radio_observation_type,
                start_date_ro,
                end_date_ro,
            )

        self.radio_science_urls = []

        rsi_volume_ID_list = self.get_ro_rsi_volume_ID(start_date_ro, end_date_ro, mapping_dict)

        if rsi_volume_ID_list:
            # Iterate all entries for each volume ID (not just the first)
            for rsi_id in rsi_volume_ID_list:
                for entry in mapping_dict[rsi_id]:
                    try:
                        # Extract target and mission phase abbreviation (Abbn)
                        target = entry.get("target")
                        abbn = entry.get("abbn")
                        rsi_volume_ID_num = entry.get("rsi_volume_id_num", "").strip()

                        # Guard against None or empty values
                        if not target or not abbn or not rsi_volume_ID_num:
                            print(
                                f"Warning: Incomplete mapping for volume {rsi_id} "
                                f"(target={target!r}, abbn={abbn!r}, "
                                f"num={rsi_volume_ID_num!r}). Skipping."
                            )
                            continue

                        target = target.strip()
                        abbn = abbn.strip()

                        # Construct the URL using the target, Abbn and rsi_volume_id.
                        # The URL pattern:
                        # "https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/RSI/RO-X-RSI-1-2-3-PPP-RRRR-V1.0/"
                        volume_ID_url = (
                            radio_science_base_url
                            + "RO-"
                            + target
                            + "-RSI-1-2-3-"
                            + abbn
                            + "-"
                            + rsi_volume_ID_num
                            + "-V1.0/"
                        )
                        # Check if URL exists with a HEAD request
                        response = requests.head(volume_ID_url, timeout=_REQUEST_TIMEOUT)
                        if response.status_code == 200:
                            print(f"URL Exists: {volume_ID_url}")
                            if volume_ID_url not in self.radio_science_urls:
                                self.radio_science_urls.append(volume_ID_url)
                        else:
                            print(f"URL does not exist: {volume_ID_url}")
                    except Exception as e:
                        print(f"Error occurred for rsi_volume_id {rsi_id}: {e}")
                        continue

        if len(self.radio_science_urls) > 0:
            return self.radio_science_urls
        else:
            raise ValueError(
                f"No url available for RO radio science files. Please check the mapping."
            )

    #########################################################################################################

    def get_ro_rsi_volume_ID(self, start_date, end_date, mapping_dict):
        """
        Given a start_date and end_date, iterate over the mapping_dict (which is keyed by rsi_volume_id)
        and return a list of rsi_volume_id values whose associated record's start_date_utc falls within the interval.

        Parameters:
            start_date (datetime): The start of the input date interval.
            end_date (datetime): The end of the input date interval.
            mapping_dict (dict): A dictionary where the keys are rsi_volume_id strings and the values
                                are lists of dictionaries. Each dictionary includes at least:
                                - "rsi_volume_id": the RSI volume ID
                                - "start_date_utc": a datetime object representing the record's start date.

        Returns:
            list: A list of rsi_volume_id strings that fall within the [start_date, end_date] interval.

        Raises:
            ValueError: If no matching rsi_volume_id is found.
        """
        rsi_volume_id_list = []

        # Iterate over each key-value pair. In this dictionary, keys are rsi_volume_id.
        for key, items in mapping_dict.items():
            # Ensure items is always a list.
            if not isinstance(items, list):
                items = [items]

            for item in items:
                try:
                    record_date = item.get("start_date_utc")
                    if record_date is None:
                        # Log or skip entries without a valid date.
                        print(f"Warning: Missing start_date_utc in item: {item}")
                        continue
                    # Check if the record_date falls within the input interval.
                    if start_date <= record_date <= end_date:
                        rsi_volume_id_list.append(key)
                except Exception as e:
                    print(f"Error processing item {item}: {e}")

        # If no valid volume IDs are found, raise an error.
        if rsi_volume_id_list:
            return list(dict.fromkeys(rsi_volume_id_list))
        else:
            raise ValueError(
                f"No RSI Volume_ID found associated with input interval: {start_date} - {end_date}."
            )

    #########################################################################################################

    def add_ro_mission_phase_designation(self, mapping_dict):
        """
        Updates mapping_dict entries by adding an 'Abbn' field (mission phase abbreviation)
        based on each entry's start_date_utc.
        The function caches the mapping from unique start dates to their mission phase.

        Parameters:
            mapping_dict (dict): Dictionary where keys are rsi_volume_id and values are lists of dictionaries, each having a 'start_date_utc' datetime object.

        Returns:
        dict: The updated mapping_dict with two additional keys for each entry:

              - "Abbn": the mission phase abbreviation.
              - "target": the target designation ("X" for non-target-specific or pre-comet phases, "C" for comet, "M" for Mars, "A" for asteroid flybys).

        """

        # Define the mission phase date ranges for the Rosetta mission, in which RSI experiments were conducted.
        mission_phases = {
            "CVP1": {"start": "2004-03-05", "end": "2004-06-06"},
            "CVP2": {"start": "2004-09-06", "end": "2004-10-16"},
            "CR2": {"start": "2005-04-05", "end": "2006-07-28"},
            "MARS": {"start": "2006-07-29", "end": "2007-05-28"},
            "EAR2": {"start": "2007-09-13", "end": "2008-01-27"},
            "CR4B": {"start": "2008-10-06", "end": "2009-09-13"},
            "AST2": {"start": "2010-05-17", "end": "2010-09-03"},
            "PRL": {"start": "2014-01-21", "end": "2014-11-19"},
            "ESC1": {"start": "2014-11-20", "end": "2015-03-10"},
            "ESC2": {"start": "2015-03-11", "end": "2015-06-30"},
            "ESC3": {"start": "2015-07-01", "end": "2015-10-21"},
            "ESC4": {"start": "2015-10-22", "end": "2015-12-31"},
            "EXT1": {"start": "2016-01-01", "end": "2016-04-05"},
            "EXT2": {"start": "2016-04-06", "end": "2016-06-30"},
            "EXT3": {"start": "2016-07-01", "end": "2016-09-30"},
        }

        for phase, dates in mission_phases.items():
            # Convert the start date and end date string using the provided method.
            dates["start_dt"] = self.format_string_to_datetime(dates["start"]).date()
            dates["end_dt"] = self.format_string_to_datetime(dates["end"]).date()

        # Cache mapping: map each unique entry start date (as date object) to its phase abbreviation.
        date_to_phase = {}
        unique_dates = set()
        for entries in mapping_dict.values():
            for entry in entries:
                dt = entry.get("start_date_utc")
                if dt is not None:
                    unique_dates.add(dt.date())

        for d in unique_dates:
            found_phase = None
            for phase, data in mission_phases.items():
                start_dt = data["start_dt"]
                end_dt = data["end_dt"]
                if start_dt <= d <= end_dt:
                    found_phase = phase
                    break
            if found_phase is None:
                print(
                    f"Warning: Date {d} does not fall within any known "
                    f"Rosetta mission phase. "
                    f"URL construction for this date will be skipped."
                )
            date_to_phase[d] = found_phase

        # Define a mapping from mission phase abbreviation to target designation.
        # (X: early/launch phases, C: comet target, M: Mars, A: asteroid)
        phase_to_target = {
            "CVP1": "X",
            "CVP2": "X",
            "CR2": "X",
            "MARS": "M",
            "EAR2": "X",
            "AST2": "A",
            "CR4B": "X",
            # After the flybys, the mission target becomes the comet (67P)
            "PRL": "C",
            "ESC1": "C",
            "ESC2": "C",
            "ESC3": "C",
            "ESC4": "C",
            "EXT1": "C",
            "EXT2": "C",
            "EXT3": "C",
        }

        # Update each entry in mapping_dict, adding both the mission phase abbreviation ("abbn")
        # and the corresponding target designation ("target").
        for entries in mapping_dict.values():
            for entry in entries:
                dt = entry.get("start_date_utc")
                if dt is not None:
                    phase = date_to_phase.get(dt.date())
                    entry["abbn"] = phase
                    entry["target"] = phase_to_target.get(phase, None)
                else:
                    entry["abbn"] = None
                    entry["target"] = None

        return mapping_dict

    #########################################################################################################

    def get_ro_rsi_volume_ID_mapping(self, url):
        """
        Description:
        Fetches data from a given URL and extracts volume ID, date range, and observation type
        for the Rosetta mission. The data is returned as a dictionary mapping date intervals
        to volume IDs and metadata.

        Inputs:
            - url (`str`): The URL from which to fetch the data (plain text format).

        Outputs:
            - `mapping_dict` (`dict`): A dictionary where keys are "rsi_volume_id",
              and values are dictionaries with:

                - `rsi_volume_id` (`str`): The rsi volume ID.
                - `volume_id` (`str`): The volume ID.
                - `start_date_file` (`str`): Start date (YYYY-MM-DD).
                - `start_date_utc` (`datetime`): Start date in UTC format.
                - `radio_observation_type` (`str`): Type of observation.
        """

        # Step 1: Fetch content from the URL (with retries)
        max_retries = 3
        for attempt in range(1, max_retries + 1):
            try:
                response = requests.get(url, timeout=_REQUEST_TIMEOUT)
                response.raise_for_status()
                break
            except requests.exceptions.RequestException as e:
                if attempt < max_retries:
                    wait_time = 2**attempt
                    print(
                        f"Attempt {attempt}/{max_retries} failed for {url}: {e}. "
                        f"Retrying in {wait_time}s..."
                    )
                    time.sleep(wait_time)
                else:
                    raise ConnectionError(
                        f"Failed to fetch RSI volume ID mapping from {url} "
                        f"after {max_retries} attempts: {e}"
                    ) from e
        aareadme_text = response.text

        # Step 2: Parse content using regex to extract the table entries
        self.mapping_dict = {}
        pattern = re.compile(
            r"^\s*(RORSI_\d{4})\s+(RORSI_\d{4})\s+(\d{4}-\d{2}-\d{2})\s+(.+?)\s*$",
            re.MULTILINE,
        )

        # Step 3: Find all matches and populate the dictionary
        for match in pattern.finditer(aareadme_text):
            full_rsi_volume_id = match.group(1)  # e.g. "RORSI_0001"
            # Remove the "RORSI_" prefix to keep only the numeric part.
            rsi_volume_id_num = full_rsi_volume_id.replace("RORSI_", "")
            volume_id = match.group(2)
            start_date_file = match.group(3) if len(match.group(3)) == 10 else None
            start_date_utc = (
                self.format_string_to_datetime(start_date_file)
                if start_date_file is not None
                else None
            )
            radio_observation_type = match.group(4).strip()

            # Add entry to dictionary
            if full_rsi_volume_id not in self.mapping_dict:
                self.mapping_dict[full_rsi_volume_id] = []  # Initialize as a list

            # Append the current entry to the list
            self.mapping_dict[full_rsi_volume_id].append(
                {
                    "rsi_volume_id_num": rsi_volume_id_num,
                    "volume_id": volume_id,
                    "start_date_file": start_date_file,
                    "start_date_utc": start_date_utc,
                    "radio_observation_type": radio_observation_type,
                }
            )
        return self.mapping_dict
