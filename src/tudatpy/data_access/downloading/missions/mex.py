"""Mars Express (MEX) file discovery and download."""

import requests
from bs4 import BeautifulSoup
import re
from ._download import _REQUEST_TIMEOUT


class MexMixin:
    def get_mex_files(self, local_folder, start_date, end_date, radio_observation_type=None):
        """
        Description:
        This function downloads and organizes various types of SPICE kernels and
        radio science files for the MEX (Mars Express) mission.
        It supports downloading kernel files related to radio science, clock data, frame data,
        SPK (spacecraft position) data, CK (orientation) data, and tropospheric/ionospheric correction files.
        The function interacts with remote FTP servers to retrieve the data and organizes them into the specified local folder.
        The data is filtered by a given date range, and the function returns the downloaded files for further use.

        Inputs:
            - local_folder (`str`): The local directory where the downloaded files will be saved.
            - start_date (`datetime`): The start date for downloading data.
              This will filter the data to include only those within the date range.
            - end_date (`datetime`): The end date for downloading data.
              This will filter the data to include only those within the date range.
            - radio_observation_type (`str`): The type of radio science files to download (e.g. phobos gravity, commissioning, etc...)

        Outputs:
            - (`dict`, `dict`, `dict`): A tuple containing:

                - `kernel_files_to_load` (`dict`): A dictionary where the keys are kernel types (e.g., 'ck', 'spk', 'fk', 'sclk') and values are lists of paths to the successfully downloaded and loaded kernel files.
                - `radio_science_files_to_load` (`dict`): A dictionary where keys are categories of radio science data (e.g., 'ifms_dp2', 'dsn_dps') and values are lists of paths to the successfully downloaded radio science files.
                - `ancillary_files_to_load` (`dict`): A dictionary where keys are categories of ancillary data (e.g., 'ion', 'tropospheric') and values are lists of paths to the successfully downloaded ancillary files, such as tropospheric and ionospheric corrections.
        """

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}  # empty for now

        input_mission = "mex"
        # Tropospheric corrections
        print("================================================================")
        print(f"Download {input_mission.upper()} Tropospheric and Ionospheric Corrections Files")
        url_tropo_files = self.get_url_mex_radio_science_files(
            start_date, end_date, radio_observation_type
        )
        for url_tropo_file_new in url_tropo_files:
            for folder_type in [
                "calib/closed_loop/ifms/met/",
                "calib/closed_loop/dsn/ion/",
                "calib/closed_loop/dsn/tro/",
            ]:
                url_flag = False
                url_tropo_file = url_tropo_file_new + folder_type
                response = requests.get(url_tropo_file)
                if response.status_code == 200:
                    url_flag = True
                    html = response.text
                    # Parse the HTML with BeautifulSoup
                    soup = BeautifulSoup(html, "html.parser")
                    # Extract file links and their names
                    wanted_tropo_files = []
                    for link in soup.find_all("a"):
                        href = link.get("href")
                        if href.endswith(".tab") or href.endswith(".aux"):
                            wanted_tropo_files.append(href.split("/")[-1])
                    tropo_files_to_load = self.get_kernels(
                        input_mission=input_mission,
                        url=url_tropo_file,
                        wanted_files_patterns=["*l3l1b*.tab"],
                        custom_output=local_folder,
                    )
                else:
                    continue

                if tropo_files_to_load:
                    key = folder_type.split("/")[-2]
                    self.ancillary_files_to_load[key] = tropo_files_to_load

                else:
                    print("No tropospheric or ionospheric files to download this time.")

        print("================================================================")
        print(f"Download {input_mission.upper()} Radio Science Kernels:")
        url_radio_science_files = self.get_url_mex_radio_science_files(
            start_date, end_date, radio_observation_type
        )
        for url_radio_science_file_new in url_radio_science_files:
            for closed_loop_type in ["ifms/dp2/", "dsn/dps/", "dsn/dpx/"]:
                try:
                    url_radio_science_file = (
                        url_radio_science_file_new + "data/level02/closed_loop/" + closed_loop_type
                    )
                    files = self.dynamic_download_url_files_single_time(
                        input_mission,
                        local_path=local_folder,
                        start_date=start_date,
                        end_date=end_date,
                        url=url_radio_science_file,
                    )
                    key = f"{closed_loop_type.split('/')[0]}_{closed_loop_type.split('/')[1]}"
                    self.radio_science_files_to_load[key] = files
                except Exception as e:
                    print(
                        f"Error downloading {closed_loop_type} radio science files from {url_radio_science_file}: {e}"
                    )
                    continue

        # Clock files
        print("================================================================")
        print(f"Download {input_mission.upper()} Clock Kernels:")
        url_clock_files = "https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/sclk/"
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
        url_frame_files = "https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/fk/"
        wanted_frame_files = [
            "MEX_V16.TF",
            "MEX_SCI_V01.TF",
            "MEX_RELAY_LOCATIONS_V03.TF",
            "MEX_DSK_SURFACES_V04.TF",
            "MEX_PFS_ROIS_V02.TF",
        ]
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
        url_spk_files = ["https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/spk/"]

        spk_files_to_load = []

        if len(url_spk_files) == 1:
            spk_files_to_load = self.dynamic_download_url_files_time_interval(
                input_mission,
                local_path=local_folder,
                start_date=start_date,
                end_date=end_date,
                url=url_spk_files[0],
            )

        else:
            for url_spk_file in url_spk_files:
                spk_files_to_load = self.dynamic_download_url_files_time_interval(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=url_spk_file,
                )

        if spk_files_to_load:
            self.kernel_files_to_load["spk"] = spk_files_to_load
        else:
            print("No spk files to download this time.")

            # Orientation files
        print("================================================================")
        print(f"Download {input_mission.upper()} CK Kernels:")
        url_ck_files = ["https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/ck/"]

        ck_files_to_load = []
        if len(url_ck_files) == 1:
            ck_files_to_load = self.dynamic_download_url_files_time_interval(
                input_mission,
                local_path=local_folder,
                start_date=start_date,
                end_date=end_date,
                url=url_ck_files[0],
            )

        else:
            for url_ck_file in url_ck_files:
                ck_files_to_load = self.dynamic_download_url_files_time_interval(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=url_ck_file,
                )

        if ck_files_to_load:
            self.kernel_files_to_load["ck"] = ck_files_to_load
        else:
            print("No spk files to download this time.")

        print("----------------------------------------------------------------")
        print(
            "All requested, relevant and previously non-existing MEX files have been now downloaded. Enjoy!"
        )

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    #########################################################################################################

    def get_mex_volume_ID(self, start_date, end_date, interval_dict):
        self.volume_id_list = []
        for (key_start, key_end), items in interval_dict.items():
            if not isinstance(items, list):
                items = [items]

            for item in items:
                # Check if either start or end of the input interval overlaps with the dictionary interval
                if key_end:
                    if start_date <= key_end and end_date >= key_start:
                        self.volume_id_list.append(item["volume_id"])
                if end_date >= key_start:
                    self.volume_id_list.append(item["volume_id"])

        if self.volume_id_list:
            return self.volume_id_list
        else:
            raise ValueError(
                f"No MEX Volume_ID found associated to input interval: {start_date} - {end_date}."
            )

    #########################################################################################################

    def get_url_mex_radio_science_files(
        self, start_date_mex, end_date_mex, radio_observation_type=None
    ):

        # url = "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-v1/mexmrs_0735/aareadme.txt"
        url = (
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext9-v1/mexmrs_4405/aareadme.txt"
        )
        radio_science_base_url = "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-v1/"
        radio_science_base_urls = [
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext1-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext2-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext3-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext4-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext5-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext6-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext7-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext8-v1/",
            "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-ext9-v1/",
        ]

        mapping_dict = self.get_mex_volume_ID_mapping(url)

        if radio_observation_type:
            mapping_dict = self.filter_mapping_dict_by_radio_observation_type(
                mapping_dict,
                radio_observation_type,
                start_date_mex,
                end_date_mex,
            )

        self.radio_science_urls = []
        volume_ID_list = self.get_mex_volume_ID(start_date_mex, end_date_mex, mapping_dict)
        if volume_ID_list:
            for volume_ID in volume_ID_list:
                for radio_science_base_url in radio_science_base_urls:
                    try:
                        volume_ID_url = radio_science_base_url + volume_ID + "/"
                        response = requests.head(
                            volume_ID_url, timeout=_REQUEST_TIMEOUT
                        )  # Use HEAD to check existence without downloading the content
                        if response.status_code == 200:
                            print(f"URL Exists: {volume_ID_url}")
                            self.radio_science_urls.append(volume_ID_url)
                        else:
                            print(f"URL does not exist: {volume_ID_url}")
                    except requests.exceptions.RequestException as e:
                        print(
                            f"Error occurred for radio science base url {radio_science_base_url} and volume ID {volume_ID}: {e}"
                        )
                        continue

        if len(self.radio_science_urls) > 0:
            return self.radio_science_urls
        else:
            raise ValueError(
                f"No url available for MEX radio science files. Please check the mapping."
            )

    #########################################################################################################

    def filter_mapping_dict_by_radio_observation_type(
        self, mapping_dict, radio_observation_type, start_date_mex, end_date_mex
    ):
        """
        Description:

        Filters a mapping dictionary to extract entries based on a specified observation type and a date range.
        The function returns a new dictionary where each key corresponds to a filtered set of entries that match
        the specified observation type and fall within the given start and end dates.

        Inputs:
            - mapping_dict (`dict`): A dictionary where keys represent categories and values are lists of entries. Each entry is expected to be a dictionary containing:
              - `start_date_utc` (`str`): The start date in UTC (format: YYYY-MM-DD).
              - `radio_observation_type` (`str`): The type of observation.

            - radio_observation_type (`str`): The type of observation to filter by (e.g., 'Phobos Gravity').
            - start_date_mex (`str`): The start date for filtering in UTC (format: YYYY-MM-DD).
            - end_date_mex (`str`): The end date for filtering in UTC (format: YYYY-MM-DD).

        Outputs:
            - `filtered_dict` (`dict`): A dictionary where keys are the same as in `mapping_dict`, and values are lists of filtered entries that match the specified observation type and fall within the date range.
        """

        filtered_dict = {
            key: [
                entry
                for entry in values
                if entry["radio_observation_type"] == radio_observation_type
                and start_date_mex <= entry["start_date_utc"] <= end_date_mex
            ]
            for key, values in mapping_dict.items()
            if any(
                entry["radio_observation_type"] == radio_observation_type
                and start_date_mex <= entry["start_date_utc"] <= end_date_mex
                for entry in values
            )
        }
        return filtered_dict

    #########################################################################################################

    def get_mex_volume_ID_mapping(self, url):
        """
        Description:
        Fetches data from a given URL and extracts volume ID, date range, and observation type
        for the Mars Express mission. The data is returned as a dictionary mapping date intervals
        to volume IDs and metadata.

        Inputs:
            - url (`str`): The URL from which to fetch the data (plain text format).

        Outputs:
            - `mapping_dict` (`dict`): A dictionary where keys are tuples of (start_date_utc, end_date_utc), and values are dictionaries with:

                - `volume_id` (`str`): The volume ID.
                - `start_date_file` (`str`): Start date (YYYY-MM-DD).
                - `end_date_file` (`str`): End date (YYYY-MM-DD).
                - `radio_observation_type` (`str`): Type of observation.

        """

        # Step 1: Fetch content from the URL
        response = requests.get(url)
        response.raise_for_status()  # Check for request errors
        aareadme_text = response.text

        # Step 2: Parse content using regex to extract the table entries
        self.mapping_dict = {}
        pattern = re.compile(
            r"^\s*(MEXMRS_\d{4})\s+(MEXMRS_\d{4}|\d{4}-\d{2}-\d{2})\s+(\d{4}-\d{2}-\d{2})\s+(.+?)\s*$",
            re.MULTILINE,
        )

        # Step 3: Find all matches and populate the dictionary
        for match in pattern.finditer(aareadme_text):
            volume_id = match.group(1)
            start_date_file = match.group(2) if len(match.group(2)) == 10 else match.group(3)
            end_date_file = match.group(3) if len(match.group(2)) == 10 else None
            start_date_utc = self.format_string_to_datetime(start_date_file)
            end_date_utc = (
                self.format_string_to_datetime(end_date_file) if end_date_file != None else None
            )
            interval_key_for_retrieval = (start_date_utc, end_date_utc)
            radio_observation_type = match.group(4).strip()

            # Add entry to dictionary
            if interval_key_for_retrieval not in self.mapping_dict:
                self.mapping_dict[interval_key_for_retrieval] = []  # Initialize as a list

            # Append the current entry to the list
            self.mapping_dict[interval_key_for_retrieval].append(
                {
                    "volume_id": volume_id,
                    "start_date_file": start_date_file,
                    "end_date_file": end_date_file,
                    "start_date_utc": start_date_utc,
                    "end_date_utc": end_date_utc,
                    "radio_observation_type": radio_observation_type,
                }
            )
        return self.mapping_dict

    ########################################################################################################################################
    ################################################# END OF MEX SECTION ###############################################################
    ########################################################################################################################################
