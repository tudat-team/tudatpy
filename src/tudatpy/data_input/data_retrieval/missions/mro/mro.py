"""Mars Reconnaissance Orbiter (MRO) file discovery and download."""

import requests
from datetime import datetime
import re


class MroMixin:
    """Mars Reconnaissance Orbiter mission file-discovery mixin."""

    def get_mro_files(self, local_folder, start_date, end_date, radio_science_file_type="odf"):
        """
        Description:
        Downloads various SPICE kernel and ancillary files for the MRO mission, including clock, frame,
        orientation, and SPK kernels, based on a specified date range. Files are saved in the provided local folder.

        Inputs:
            - local_folder (`str`): Path to the local folder where files will be saved.
            - start_date (`datetime`): The start date for the data retrieval.
            - end_date (`datetime`): The end date for the data retrieval.

        Outputs:
            - `kernel_files_to_load` (`dict`): A dictionary containing the loaded kernel files, categorized by type
              (e.g., 'sclk', 'fk', 'ck', 'spk').
            - `radio_science_files_to_load` (`dict`): An empty dictionary for now, intended for radio science files.
            - `ancillary_files_to_load` (`dict`): An empty dictionary for now, intended for ancillary files.
        """

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}

        input_mission = "mro"
        # ODF files
        print("================================================================")
        print(f"Download {input_mission.upper()} {radio_science_file_type.upper()} files:")
        url_radio_science_file_type = [
            "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/{}/".format(
                radio_science_file_type
            )
        ]
        key = radio_science_file_type
        for url_radio_science_file in url_radio_science_file_type:
            try:
                files = self.dynamic_download_url_files_single_time(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=url_radio_science_file,
                )
                self.radio_science_files_to_load[key] = files
            except Exception as e:
                print(
                    f"Error downloading mro radio science file from {url_radio_science_file}: {e}"
                )
                continue

        if not self.radio_science_files_to_load:
            print("No Radio Science files to download this time.")
        # Clock Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Clock Kernels:")
        url_clock_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/sclk/"
        )
        wanted_clock_files = self.get_latest_clock_kernel_name(input_mission)
        clock_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_clock_files,
            wanted_files=wanted_clock_files,
            custom_output=local_folder,
        )

        try:

            if clock_files_to_load:
                self.kernel_files_to_load["sclk"] = clock_files_to_load
            else:
                print("No sclk files to download this time.")
        except Exception as e:
            print(f"Error handling clock file downloads: {e}")
            raise

        # Frame Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Frame Kernels:")
        url_frame_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/fk/"
        )
        wanted_frame_files = ["mro_v16.tf"]
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

            # Planetary and Ephemeris Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} SPK Kernels:")
        url_spk_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/spk/"
        )
        wanted_spk_files = self.get_url_mro_spk_files(start_date, end_date)
        spk_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_spk_files,
            wanted_files=wanted_spk_files,
            custom_output=local_folder,
        )
        wanted_struct_files = ["mro_struct_v10.bsp"]
        struct_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_spk_files,
            wanted_files=wanted_struct_files,
            custom_output=local_folder,
        )

        spk_files_to_load.extend(struct_files_to_load)

        if spk_files_to_load:
            self.kernel_files_to_load["spk"] = spk_files_to_load
        else:
            print("No spk files to download this time.")

            # Orientation Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Orientation Kernels:")
        measured_url_ck_files = [
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/ck/"
        ]

        if len(measured_url_ck_files) == 1:
            measured_ck_files_to_load = self.dynamic_download_url_files_time_interval(
                input_mission,
                local_path=local_folder,
                start_date=start_date,
                end_date=end_date,
                url=measured_url_ck_files[0],
            )
        else:
            for measured_url_ck_file in measured_url_ck_files:
                measured_ck_files_to_load = self.dynamic_download_url_files_time_interval(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=measured_url_ck_file,
                )

        if measured_ck_files_to_load:
            self.kernel_files_to_load["ck"] = measured_ck_files_to_load
        else:
            print("No ck files to download this time.")

            # Tropospheric corrections
        print("================================================================")
        print(f"Download {input_mission.upper()} Tropospheric Corrections Files")
        url_tropo_files = (
            "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/tro/"
        )
        tropo_files_to_load = self.dynamic_download_url_files_time_interval(
            input_mission,
            local_path=local_folder,
            start_date=start_date,
            end_date=end_date,
            url=url_tropo_files,
        )

        if tropo_files_to_load:
            self.ancillary_files_to_load["tro"] = tropo_files_to_load
        else:
            print("No tropospheric files to download this time.")

            # Ionospheric corrections
        print("================================================================")
        print(f"Download {input_mission.upper()} Ionospheric Corrections Files")
        url_ion_files = (
            "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/ion/"
        )
        ion_files_to_load = self.dynamic_download_url_files_time_interval(
            input_mission,
            local_path=local_folder,
            start_date=start_date,
            end_date=end_date,
            url=url_ion_files,
        )

        try:

            if ion_files_to_load:
                self.ancillary_files_to_load["ion"] = ion_files_to_load
            else:
                print("No ionospheric files to download this time.")
        except Exception as e:
            print(f"Error handling file downloads: {e}")
            raise

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    ########################################################################################################################################

    def get_mapping_dict_spk_mro(self, url):
        """
        Description:
        Fetches and parses phase information for MRO SPK kernels from a given URL. It extracts the phase name
        along with the associated start and end dates, and returns a dictionary mapping date intervals to phase names.

        Inputs:
            - url (`str`): The URL from which to fetch the phase information.

        Outputs:
            - `phase_dict` (`dict`): A dictionary where the keys are tuples representing date intervals
              (start_date, end_date), and the values are dictionaries containing the SPK phase ID (`spk_ID`).
        """

        # Step 1: Fetch content from the URL
        response = requests.get(url)
        response.raise_for_status()  # Check for request errors
        file_text = response.text

        # Step 2: Define a dictionary to hold date intervals and phase names
        self.phase_dict = {}

        # Step 3: Regex pattern to capture phase name and dates
        pattern = re.compile(
            r"^\s*(\w+)\s+(\d{4} \w{3} \d{2})\s+(\d{4} \w{3} \d{2})\s*$",
            re.MULTILINE,
        )

        # Step 4: Extract all matches and populate the dictionary
        for match in pattern.finditer(file_text):
            phase = match.group(1)
            start_date_str = match.group(2)
            end_date_str = match.group(3)

            # Convert string dates to datetime objects for easier manipulation if needed
            start_date = datetime.strptime(start_date_str, "%Y %b %d")
            end_date = datetime.strptime(end_date_str, "%Y %b %d")

            # Use the date interval as the key and the phase as the value
            self.phase_dict[(start_date, end_date)] = {"spk_ID": phase}

        return self.phase_dict

    ########################################################################################################################################

    def get_mro_spk_ID(self, start_date, end_date, mapping_dict):
        """
        Description:
        Checks if the given date interval overlaps with any of the date intervals in the mapping dictionary,
        and returns a list of SPK IDs associated with those overlapping intervals.

        Inputs:
            - start_date (`datetime`): The start date of the input interval.
            - end_date (`datetime`): The end date of the input interval.
            - mapping_dict (`dict`): A dictionary mapping date intervals (as tuples of start and end dates)
              to SPK phase IDs (`spk_ID`).

        Outputs:
            - `spk_id_list` (`list`): A list of SPK IDs associated with the overlapping intervals.

        Raises:
            - `ValueError`: If no SPK IDs are found for the given date interval.
        """

        self.spk_id_list = []
        for (key_start, key_end), item in mapping_dict.items():
            # Check if either start or end of the input interval overlaps with the dictionary interval
            if start_date <= key_end and end_date >= key_start:
                self.spk_id_list.append(item["spk_ID"])

        if self.spk_id_list:
            return self.spk_id_list
        else:
            raise ValueError(
                f"No MRO SPK_ID found associated to input interval: {start_date} - {end_date}."
            )

    ########################################################################################################################################

    def get_url_mro_spk_files(self, start_date, end_date):
        """
        Description:
        Fetches the SPK file URLs associated with the given date range for the MRO mission. It first retrieves
        a mapping of SPK phases to date intervals, identifies the SPK IDs that overlap with the input interval,
        and then constructs the corresponding SPK file URLs.

        Inputs:
            - start_date (`datetime`): The start date of the input interval.
            - end_date (`datetime`): The end date of the input interval.

        Outputs:
            - `spk_ID_urls` (`list`): A list of SPK file URLs corresponding to the overlapping SPK IDs.

        Raises:
            - `ValueError`: If no SPK file URLs are found for the given date range.
        """

        url = "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/spk/spkinfo.txt"
        mapping_dict = self.get_mapping_dict_spk_mro(url)

        self.spk_ID_urls = []
        spk_ID_list = self.get_mro_spk_ID(start_date, end_date, mapping_dict)
        if spk_ID_list:
            for spk_ID in spk_ID_list:
                spk_ID_url = f"mro_{spk_ID}.bsp"
                self.spk_ID_urls.append(spk_ID_url)

        if len(self.spk_ID_urls) > 0:
            return self.spk_ID_urls
        else:
            raise ValueError(
                f"No url available for MEX radio science files. Please check the mapping."
            )

    ########################################################################################################################################
    ##################################################### END OF MRO SECTION ###############################################################
    ########################################################################################################################################
