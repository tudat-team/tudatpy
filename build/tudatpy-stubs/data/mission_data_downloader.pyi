import typing
from _typeshed import Incomplete
from tudatpy.interface import spice as spice

class LoadPDS:
    flag_check_existing_files: bool
    flag_load_standard_kernels: bool
    supported_patterns: Incomplete
    cassini_titan_flyby_dict: Incomplete
    type_to_extension: Incomplete
    supported_mission_odf_time_formats: Incomplete
    supported_mission_meta_kernel_url: Incomplete
    supported_mission_meta_kernel_pattern: Incomplete
    supported_mission_kernels_url: Incomplete

    def __init__(self) -> None:
        ...

    def print_titan_flyby_table(self) -> None:
        """
        Description:
        This method prints a table displaying the Titan flyby data in a readable format. It iterates over the `cassini_titan_flyby_dict`, extracts relevant information for each flyby, and constructs a table. The table is formatted with headers and colored using the `colorama` library for enhanced readability.

        Input:
            None.

        Output:
            None: This method prints the formatted table to the console.
        """

    def spice_transfer2binary(self, input_file, timeout: int=5):
        """
        Description:
            Converts transfer-file format SPICE kernels (e.g., .ckf, .spk)
            into binary SPICE kernels (e.g., .ck, .bsp) using the SPACIT utility.
            This is necessary for loading the kernels with spice.load_kernels.
            The function handles the conversion by running the SPACIT tool as a subprocess and
            manages timeouts during the conversion process.

        Input:
            - input_file (str): The path to the SPICE kernel file to be converted. It must have either a .ckf or .spk extension.
            - timeout (int, optional): The timeout duration for the SPACIT process, default is 5 seconds.

        Output:
            - str: The path to the output file (either .ck or .bsp), depending on the input file type. If the output file already exists, the function returns the existing output file path.

        Notes:
            If the conversion fails or times out, an error message is printed.
        """

    def create_pattern(self, placeholders):
        ...
    supported_time_formats: Incomplete
    date: Incomplete

    def format_string_to_datetime(self, date):
        """
        Description:
            Attempts to convert a given date string into a `datetime` object using multiple supported time formats. The function iterates over a list of predefined date formats (e.g., `%y%m%d`, `%Y%j`, `%Y_%j`) and tries each format until it successfully parses the date. This is used throughout the code to handle various date formats from different mission data files.

        Input:
            - `date` (`str`): The date string to be parsed.

        Output:
            - `datetime`: The parsed `datetime` object if a matching format is found, or `None` if no format matches.

        Notes:
            The function supports multiple date formats, such as:
            - `YYMMDD` (e.g., `220101`)
            - `YYYYjjj` (e.g., `2022323`)
            - `YYYY_MM_DD` (e.g., `2022_03_15`)
            - And others, tailored for specific mission data.
            If no matching format is found, the function silently skips to the next format without throwing an error.
        """

    def format_datetime_to_string(self, date, format_key):
        """
        Description:
            Converts a `datetime` object to a string based on a specified format key from a list of supported formats. The function maps the provided format key (e.g., "YYMMDD", "YYYY-jjj") to a corresponding `strftime` format string and then formats the given `datetime` object accordingly.

        Input:
            - `date` (`datetime`): The `datetime` object to be formatted.
            - `format_key` (`str`): A key from the `supported_time_formats` dictionary that specifies the desired output format.

        Output:
            - `str`: The formatted date string according to the specified format key.

        Notes:
            The function supports several predefined date formats, such as:
            - `"YYMMDD"`: e.g., `220101`
            - `"YYYYjjj"`: e.g., `2022323`
            - `"YYYY-MM-DD"`: e.g., `2022-03-15`
            If the provided `format_key` is not in the supported formats, an error message is printed.
        """
    files_to_load: Incomplete

    def get_kernels(self, input_mission, url, wanted_files: Incomplete | None=None, wanted_files_patterns: Incomplete | None=None, custom_output: Incomplete | None=None):
        """
        Downloads specific SPICE kernel files or files matching a pattern from a given URL to a local directory
        if they do not already exist locally.

        Input:
            - `input_mission` (`str`): The name of the mission
            - `url` (`str`): The base URL where the kernel files are hosted.
            - `wanted_files` (`list`, optional): A list of specific filenames to be downloaded from the URL.
            - `wanted_files_pattern` (`str`, optional): A pattern (e.g., '\\*.tf') to match filenames for downloading.
            - `custom_output` (`str`, optional): The local directory where the downloaded files will be stored.

        Output:
            - `list`: A list of full file paths for the downloaded (or already existing) kernel files.
        """

    def add_custom_mission_kernel_pattern(self, input_mission, custom_kernel_type, custom_pattern):
        """
        Description:
            Allows users to define and add custom regex patterns for a specific mission to the list of supported patterns. Once added, the custom pattern can be used for mission data file matching.

        Input:
            - `input_mission` (`str`): The name of the mission for which the custom pattern is defined.
            - `custom_kernel_type`: the name of kernel type (eg 'ck'., 'spk', etc...)
            - `custom_pattern` (`str`): The regular expression pattern that will be associated with the mission kernel.

        Output:
            - `dict`: The updated `supported_patterns` dictionary containing the new custom pattern.

        Notes:
            - This function updates the `supported__kernel_patterns` dictionary with the new mission and its associated pattern.
            - Custom patterns can be added dynamically at runtime.
        """

    def add_custom_mission_meta_kernel_url(self, input_mission, url):
        """
        Description:
            Allows users to define and add custom regex patterns for a specific mission to the list of supported patterns. Once added, the custom pattern can be used for mission data file matching.

        Input:
            - `input_mission` (`str`): The name of the mission for which the custom pattern is defined.
            - `url` (`str`): The url that will be associated with the mission.

        Output:
            - `dict`: The updated `supported_mission_meta_kernel_url` dictionary containing the new custom pattern.

        Notes:
            - This function updates the `supported_mission_meta_kernel_url` dictionary with the new mission and its associated pattern.
            - Custom patterns can be added dynamically at runtime.
        """

    def add_custom_mission_kernel_url(self, input_mission, url):
        """
        Description:
            Allows users to define and add custom regex patterns for a specific mission to the list of supported patterns. Once added, the custom pattern can be used for mission data file matching.

        Input:
            - `input_mission` (`str`): The name of the mission for which the custom pattern is defined.
            - `url` (`str`): The url that will be associated with the mission.

        Output:
            - `dict`: The updated `supported_mission_kernels_url` dictionary containing the new custom pattern.

        Notes:
            - This function updates the `supported_mission_kernels_url` dictionary with the new mission and its associated pattern.
            - Custom patterns can be added dynamically at runtime.
        """

    def add_custom_mission_meta_kernel_pattern(self, input_mission, custom_pattern):
        """
        Description:
            Allows users to define and add custom regex patterns for a specific mission to the list of supported patterns. Once added, the custom pattern can be used for mission data file matching.

        Input:
            - `input_mission` (`str`): The name of the mission for which the custom pattern is defined.
            - `custom_pattern` (`str`): The regular expression pattern that will be associated with the mission.

        Output:
            - `dict`: The updated `supported_mission_kernel_pattern` dictionary containing the new custom pattern.

        Notes:
            - This function updates the `supported_mission_kernel_pattern` dictionary with the new mission and its associated pattern.
            - Custom patterns can be added dynamically at runtime.
        """
    in_intervals: bool

    def is_date_in_intervals(self, date, intervals):
        """
        Description:
            Checks if a given date falls within any of the specified date intervals.

        Input:
            - `date` (`datetime`): The date to check.
            - `intervals` (`list` of tuples): A list of tuples, where each tuple contains two `datetime` objects representing the start and end of an interval.

        Output:
            - `bool`: `True` if the date falls within any of the intervals, otherwise `False`.
        """

    def clean_mission_archive(self, local_folder) -> None:
        """
        Description:
            Cleans up the mission archive by removing empty subdirectories within the specified local folder.

        Input:
            - `local_folder` (`str`): Path to the local folder where the archive is stored.

        Output:
            - `None`: This function does not return any value. It performs in-place directory cleanup.
        """

    def match_type_extension(self, data_type, filename):
        """
        Checks if the extension of a given file matches the expected extension for the specified data type.

        Args:
            data_type (str): The data type (e.g., 'ck', 'spk', 'odf', etc.).
            filename (str): The filename whose extension is to be checked.

        Returns:
            bool: True if the file extension matches any valid extension for the given data type.
        """

    def get_extension_for_data_type(self, data_type, first_only: bool=True):
        """
        Returns one or more file extensions for the given data type.

        Args:
            data_type (str): The data type (e.g., 'ck', 'spk', etc.).
            first_only (bool): If True, returns only the first matching extension.
                            If False, returns a list of all valid extensions.

        Returns:
            str or list[str] or None: Extension(s) for the given data type, or None if unsupported.
        """
    relevant_files: Incomplete

    def dynamic_download_url_files_single_time(self, input_mission, local_path, start_date, end_date, url):
        """
        Description:
            Downloads files for a single time interval based on specified mission parameters and date range.
            It checks if the files already exist locally, and if not, attempts to download them from the specified URL.
            The filenames are matched to supported patterns, and relevant files are downloaded if they are not already present.

        Inputs:
            - input_mission (`str`): The name of the mission (e.g., 'cassini', 'mro').
            - local_path (`str`): The local directory where the files should be saved.
            - start_date (`datetime`): The start date of the time interval for which files are required.
            - end_date (`datetime`): The end date of the time interval for which files are required.
            - url (`str`): The base URL where the files are hosted.

        Outputs:
            - `self.relevant_files` (`list`): A list of paths to the files that were either already present or successfully downloaded.
        """
    existing_files: Incomplete

    def check_existing_files(self, data_type, local_subfolder, start_date, end_date):
        """
        Description:
            Checks the local directory for files that match a given pattern and fall within a specified date range.
            This function filters and returns the files that already exist locally, based on their extensions and matching patterns.

        Inputs:
            - data_type (`str`): The type of data (e.g., 'ck', 'spk').
            - local_subfolder (`str`): Path to the local directory where the files are stored.
            - start_date (`datetime`): The start date of the time interval for which files are required.
            - end_date (`datetime`): The end date of the time interval for which files are required.

        Outputs:
            - `self.existing_files` (`list`): A list of paths to the files that already exist and match the given pattern.
        """

    def dynamic_download_url_files_time_interval(self, input_mission, local_path, start_date, end_date, url):
        """
        Description:
            Downloads files within a specified time interval from a given URL, checking if the files already exist locally before downloading them.
            This function is typically used to download mission data files for a specific time range, ensuring that files are only downloaded if they are not already present.

        Inputs:
            - input_mission (`str`): The name of the mission (e.g., 'cassini', 'mro').
            - local_path (`str`): The directory path where files will be stored locally.
            - start_date (`datetime`): The start date of the time interval for which files are required.
            - end_date (`datetime`): The end date of the time interval for which files are required.
            - url (`str`): The base URL for downloading the files.

        Outputs:
            - `self.relevant_files` (`list`): A list of local file paths for the files that were either found or successfully downloaded.
        """

    def download_url_files_time(self, local_path, filename_format, start_date, end_date, url, time_format, indices_date_filename):
        """
        Description:
            Downloads files within a specified time interval from a given URL. The function checks for existing files locally,
            identifies missing files for the given time interval, and downloads them if necessary.
            It handles cases where files are organized in nested folders at the target URL.

        Inputs:
            - local_path (`str`): The local directory path where files will be stored.
            - filename_format (`str`): The format of the filenames to be downloaded, with date placeholders and optional folder structure.
            - start_date (`datetime`): The start date of the time interval for which files are needed.
            - end_date (`datetime`): The end date of the time interval for which files are needed.
            - url (`str`): The base URL from which the files will be downloaded.
            - time_format (`str`): The format string used to represent dates in the filenames (e.g., '%Y%m%d').
            - indices_date_filename (`list[int]`): A list of indices indicating where date strings are embedded within the filename format.

        Outputs:
            - `self.relevant_files` (`list`): A list of local file paths for files that were either already present or successfully downloaded.

        Raises:
            - `Exception`: If the `filename_format` contains more than one folder.

        """
    start_date_utc: Incomplete
    end_date_utc: Incomplete
    date_utc: Incomplete

    def parse_filename(self, input_mission, data_type, filename):
        """
        Description:
        This function parses a filename into its components based on the mission and data type patterns.
        It uses regular expressions to match the filename against predefined patterns for different mission types.
        It then returns the parsed components as a dictionary, along with the indices of underscores within the matched groups.

        Inputs:
            - input_mission (`str`): The mission name (e.g., \'cassini\', \'mro\', etc.).
            - data_type (`str`): The type of data (e.g., "radioscience", "kernels", "all").
            - filename (`str`): The filename to be parsed.

        Outputs:
            - (`dict`, `list`): A tuple containing:
              - `dictionary` (`dict`): A dictionary where the keys are the component names
              (e.g., \'date\', \'purpose\', etc.), and the values are the parsed values from the filename.
              - `underscore_indices` (`list`): A list of indices representing
              the positions of underscores in the matched groups.
        """

    def reconstruct_filename(self, dictionary, underscore_indices):
        """
        Description:
            This function reconstructs the original filename by combining the parsed components (from the `dictionary`) and placing underscores at the specified indices (from the `underscore_indices`). It takes into account the specific format of the filename, ensuring the proper placement of underscores.

        Inputs:
            - dictionary (`dict`): A dictionary containing the parsed components of the filename, where the keys are the component names and the values are the corresponding values extracted from the filename.
            - underscore_indices (`list`): A list of indices where underscores should be placed in the reconstructed filename. These indices correspond to the position of the groups in the parsed dictionary.

        Outputs:
            - (`str`): The reconstructed filename as a string, with underscores placed at the specified positions.
        """
    kernel_files_names: Incomplete

    def extract_kernels_from_meta_kernel(self, input_mission):
        """
        Fetches a meta-kernel file from an HTTPS URL and categorizes kernel files by type
        using the type-to-extension mapping.

        Parameters:
            input_mission (str): e.g. 'mex', 'mro', 'juice', 'ro', etc...

        Returns:
            dict: A dictionary categorizing kernel files by type based on `type_to_extension`.
        """
    kernel_files_to_load: Incomplete

    def download_kernels_from_meta_kernel(self, input_mission, local_folder):
        ...
    latest_kernel: Incomplete

    def get_latest_meta_kernel(self, input_mission):
        """
        Finds the most recent meta-kernel file URL based on year and version.

        Parameters:
            base_url (str): The base URL containing meta-kernel file links.

        Returns:
            str: The URL of the most recent meta-kernel file.
        """

    def get_latest_clock_kernel_name(self, input_mission):
        ...
    all_kernel_files: Incomplete
    all_radio_science_files: Incomplete
    all_ancillary_files: Incomplete

    def get_mission_files(self, input_mission, start_date: Incomplete | None=None, end_date: Incomplete | None=None, flyby_IDs: Incomplete | None=None, custom_output: Incomplete | None=None, all_meta_kernel_files: Incomplete | None=None, load_kernels: Incomplete | None=None, radio_observation_type: Incomplete | None=None, radio_science_file_type: str='odf'):
        """
        Description:
        This function downloads and organizes mission-specific data files (kernels, radio science files, and ancillary files)
        for a specified space mission. It supports downloading data for multiple missions such as Cassini, MEX, JUICE, MRO and ROSETTA.
        The function also allows for downloading files within a specific date range or for specific flybys (in the case of Cassini).
        It automatically creates the necessary folder structure for storing the downloaded data and loads the files into the SPICE kernel system.

        Inputs:
            - input_mission (`str`): The name of the space mission for which to download the data. Valid values include:
                - 'cassini'
                - 'mex'
                - 'juice'
                - 'mro'
                - 'ro'
            - start_date (`datetime`, optional): The start date for downloading data. If not provided, data is not filtered by date.
            - end_date (`datetime`, optional): The end date for downloading data (data from the end_date will be downloaded as well, meaning data from start_date <= date <= end_date).
                                                If not provided, data is not filtered by date.
            - flyby_IDs (`list` or `str`, optional): A list of flyby IDs (e.g., ['T101', 'E303']) for Cassini missions.
                It can also include special values like 'ALL_TITAN' or 'ALL_ENCELADUS' to download all flybys for Titan or Enceladus.
            - custom_output (`str`, optional): A custom path where the downloaded files will be stored. If not provided,
                the default folder structure is used based on the mission name.
            - radio_observation_type: ('str', optional): type of (mex) radio_observation_type (e.g. phobos gravity, commissioning, occultation, etc ...)
            - radio_science_file_type ('str', optional): tnf or odf

        Outputs:
            - (`dict`, `dict`, `dict`): A tuple containing:
                - `all_kernel_files` (`dict`): A dictionary where keys are kernel types (e.g., 'ck', 'spk') and values are
                  lists of paths to the successfully loaded kernel files.
                - `all_radio_science_files` (`dict`): A dictionary where keys are radio science data types and values are
                  lists of paths to the successfully loaded radio science files.
                - `all_ancillary_files` (`dict`): A dictionary where keys are ancillary data types and values are
                  lists of paths to the successfully loaded ancillary files.
        """
    radio_science_files_to_load: Incomplete
    ancillary_files_to_load: Incomplete

    def get_mex_files(self, local_folder, start_date, end_date, radio_observation_type: Incomplete | None=None):
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
            - radio_observation_type (`str`): The type of radio science files to download (e.g. phobos gravity, commissioniong, etc...)

        Outputs:
            - (`dict`, `dict`, `dict`): A tuple containing:

                - `kernel_files_to_load` (`dict`): A dictionary where the keys are kernel types (e.g., 'ck', 'spk', 'fk', 'sclk') and values are lists of paths to the successfully downloaded and loaded kernel files.
                - `radio_science_files_to_load` (`dict`): A dictionary where keys are categories of radio science data (e.g., 'ifms_dp2', 'dsn_dps') and values are lists of paths to the successfully downloaded radio science files.
                - `ancillary_files_to_load` (`dict`): A dictionary where keys are categories of ancillary data (e.g., 'ion', 'tropospheric') and values are lists of paths to the successfully downloaded ancillary files, such as tropospheric and ionospheric corrections.
        """
    volume_id_list: Incomplete

    def get_mex_volume_ID(self, start_date, end_date, interval_dict):
        ...
    radio_science_urls: Incomplete

    def get_url_mex_radio_science_files(self, start_date_mex, end_date_mex, radio_observation_type: Incomplete | None=None):
        ...

    def filter_mapping_dict_by_radio_observation_type(self, mapping_dict, radio_observation_type, start_date_mex, end_date_mex):
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

    def filter_mapping_dict_by_radio_observation_type(self, mapping_dict, radio_observation_type, start_date_mex, end_date_mex):
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
            - `filtered_dict` (`dict`): A dictionary where keys are the same as in `mapping_dict`, and values are lists
              of filtered entries that match the specified observation type and fall within the date range.
        """
    mapping_dict: Incomplete

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

    def get_juice_files(self, local_folder, start_date, end_date):
        """
        Description:
        Downloads various SPICE kernel and ancillary files for the JUICE mission, including clock, frame,
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

    def get_mro_files(self, local_folder, start_date, end_date, radio_science_file_type: str='odf'):
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
    phase_dict: Incomplete

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
    spk_id_list: Incomplete

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
    spk_ID_urls: Incomplete

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

    def get_cassini_flyby_files(self, local_folder):
        """
        Description:
        Downloads Cassini mission data files related to a specific flyby, including kernel, ancillary, and
        radio science files. The function fetches the cumulative index table and identifies the relevant
        files for download based on the provided flyby ID. It handles the file download process, ensuring
        that files are downloaded only if they are not already present in the specified local folder.

        Inputs:
            - local_folder (`str`): The path to the local directory where files should be downloaded.

        Outputs:
            - kernel_files_to_load (`dict`): A dictionary containing kernel files (e.g., `ck`, `spk`).
            - radio_science_files_to_load (`dict`): A dictionary containing radio science files (e.g., `odf`).
            - ancillary_files_to_load (`dict`): A dictionary containing ancillary files (e.g., `ion`, `tro`, `eop`).
        """

    def cassini_flyby2experiment(self, flyby_ID):
        """
        Description:
        Returns the experiment data for a given flyby ID from the `cassini_titan_flyby_dict`. If the provided
        flyby ID exists in the dictionary, it retrieves the corresponding data; otherwise, it returns a
        message indicating that the flyby ID was not found.

        Inputs:
            - flyby_ID (`str`): The identifier for a specific flyby.

        Outputs:
            - `dict` or `str`: The experiment data associated with the provided flyby ID if found,
              otherwise a message stating that the flyby ID was not found.
        """

    def get_cassini_full_moon_flybys_list(self, moon):
        """
        Description:
        Returns a list of flyby IDs corresponding to either Titan or Enceladus. If the provided
        moon name is 'titan', it returns the flyby IDs from the `cassini_titan_flyby_dict`; if the name is
        'enceladus', it returns the flyby IDs from the `enceladus_flyby_dict`. If the input moon name
        is invalid, an error is raised.

        Inputs:
            - moon (`str`): The name of the moon ('titan' or 'enceladus') for which flyby IDs should
              be retrieved.

        Outputs:
            - `list`: A list of flyby IDs corresponding to the provided moon. If the moon name is invalid,
              a `ValueError` is raised.
        """
    pds_repo_url: Incomplete

    def get_cassini_flyby_volume_url(self, flyby_dict):
        """
        Description:
        Retrieves the PDS (Planetary Data System) repository URL for a given flyby experiment. The function
        extracts the 'experiment' value from the provided `flyby_dict` and constructs the URL based on that
        information. It then checks whether the constructed URL exists by sending a HEAD request. If the URL
        is valid (status code 200), the function returns the URL; otherwise, it returns `None` and prints
        an error message.

        Inputs:
            - flyby_dict (`dict`): A dictionary containing data for a specific flyby, including the 'experiment'
              key used to construct the URL.

        Outputs:
            - `str` or `None`: The constructed PDS repository URL if it exists; otherwise, `None` if the URL
              is not found or an error occurs.
        """
    cumindex_url: Incomplete

    def get_cassini_flyby_cumindex_url(self, flyby_dict):
        """
        Description:
        Constructs the URL for the cumulative index (CUMINDEX.TAB) file associated with a specific flyby
        experiment. The function extracts the 'experiment' and 'cumindex_repo' (or 'pds_repo') values from
        the provided `flyby_dict` to build the URL. It then checks if the URL exists by sending a HEAD request.
        If the URL is valid (status code 200), it returns the URL; otherwise, it attempts an alternative URL and
        returns it if valid. If neither URL is found, it returns `None` and prints an error message.

        Inputs:
            - flyby_dict (`dict`): A dictionary containing data for a specific flyby, including the 'experiment'
              and 'cumindex_repo' (or 'pds_repo') keys used to construct the cumulative index URL.

        Outputs:
            - `str` or `None`: The constructed cumulative index URL if it exists; otherwise, `None` if the URL
              is not found or an error occurs.
        """

    def get_cassini_volume_cumindex_table(self, cumindex_url):
        """
        Description:
        Constructs the URL for the cumulative index (CUMINDEX.TAB) file associated with a specific flyby
        experiment. The function extracts the 'experiment' and 'cumindex_repo' (or 'pds_repo') values from
        the provided `flyby_dict` to build the URL. It then checks if the URL exists by sending a HEAD request.
        If the URL is valid (status code 200), it returns the URL; otherwise, it attempts an alternative URL and
        returns it if valid. If neither URL is found, it returns `None` and prints an error message.

        Inputs:
            - flyby_dict (`dict`): A dictionary containing data for a specific flyby, including the 'experiment'
              and 'cumindex_repo' (or 'pds_repo') keys used to construct the cumulative index URL.

        Outputs:
            - `str` or `None`: The constructed cumulative index URL if it exists; otherwise, `None` if the URL
              is not found or an error occurs.
        """

    def get_grail_a_files(self, local_folder, start_date, end_date):
        """
        Description:
        Downloads various SPICE kernel and ancillary files for the GRAIL-A mission, including clock, frame,
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

    def get_grail_b_files(self, local_folder, start_date, end_date):
        """
        Description:
        Downloads various SPICE kernel and ancillary files for the GRAIL_B mission, including clock, frame,
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

    def get_ro_files(self, local_folder, start_date, end_date, radio_observation_type: Incomplete | None=None):
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
            - start_date (`datetime`): The start date for downloading data.
              This will filter the data to include only those within the date range.
            - end_date (`datetime`): The end date for downloading data.
              This will filter the data to include only those within the date range.
            - radio_observation_type (`str`): The type of radio science files to download (e.g. commissioning, checkout, solar conjuction, Lutetia, Global Gravity etc...)

        Outputs:
            - (`dict`, `dict`, `dict`): A tuple containing:
                - `kernel_files_to_load` (`dict`): A dictionary where the keys are kernel types
                (e.g., 'ck', 'spk', 'fk', 'sclk') and values are lists of paths to the successfully downloaded and loaded kernel files.
                - `radio_science_files_to_load` (`dict`): A dictionary where keys are categories of
                radio science data (e.g., 'ifms_dp2', 'dsn_dps') and values are lists of paths
                to the successfully downloaded radio science files.
                - `ancillary_files_to_load` (`dict`): A dictionary where keys are categories
                of ancillary data (e.g., 'ion', 'tropospheric') and values are lists of paths
                to the successfully downloaded ancillary files, such as tropospheric and ionospheric corrections.
        """

    def get_url_ro_radio_science_files(self, start_date_ro, end_date_ro, radio_observation_type: Incomplete | None=None):
        ...

    def get_ro_rsi_volume_ID(self, start_date, end_date, mapping_dict):
        """
        Given a start_date and end_date, iterate over the mapping_dict (which is keyed by rsi_volume_id)
        and return a list of rsi_volume_id values whose associated record\'s start_date_utc falls within the interval.

        Parameters:
            start_date (datetime): The start of the input date interval.
            end_date (datetime): The end of the input date interval.
            mapping_dict (dict): A dictionary where the keys are rsi_volume_id strings and the values
                                are lists of dictionaries. Each dictionary includes at least:
                                - "rsi_volume_id": the RSI volume ID
                                - "start_date_utc": a datetime object representing the record\'s start date.

        Returns:
            list: A list of rsi_volume_id strings that fall within the [start_date, end_date] interval.

        Raises:
            ValueError: If no matching rsi_volume_id is found.
        """

    def add_ro_mission_phase_designation(self, mapping_dict):
        """
        Updates mapping_dict entries by adding an \'Abbn\' field (mission phase abbreviation)
        based on each entry\'s start_date_utc.
        The function caches the mapping from unique start dates to their mission phase.

        Parameters:
            mapping_dict (dict): Dictionary where keys are rsi_volume_id and values are lists
                                of dictionaries, each having a \'start_date_utc\' datetime object.

        Returns:
        dict: The updated mapping_dict with two additional keys for each entry:
              - "Abbn": the mission phase abbreviation.
              - "target": the target designation ("X" for non-target-specific or pre-comet phases,
                          "C" for comet, "M" for Mars, "A" for asteroid flybys).
    """

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