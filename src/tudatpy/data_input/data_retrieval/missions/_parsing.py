"""Filename-pattern construction and date/string parsing helpers shared by
every mission."""

from datetime import datetime
import re
from dateutil.relativedelta import relativedelta


class _ParsingMixin:
    def create_pattern(self, placeholders):
        """Create a regular expression from filename placeholders.

        Parameters
        ----------
        placeholders : list[str]
            Ordered filename placeholders.

        Returns
        -------
        str
            Regular-expression pattern.
        """
        final_pattern = ""
        for placeholder in placeholders:
            escaped = re.escape(placeholder)
            # Replace specific dynamic sections with regex patterns
            # Replace sequences of digits with a numeric pattern
            pattern = re.sub(r"\d+", r"[0-9]+", escaped)
            # Replace fixed-length sequences of lowercase letters (e.g., "sc", "hg") with a general pattern

            if (
                placeholder == "start_date_file"
                or placeholder == "end_date_file"
                or placeholder == "version"
            ):
                pattern = r"[0-9]+"

            else:
                pattern = re.sub(r"[a-z]{1,}", r"[a-z]{1,}", pattern, count=1)

            if placeholder == "_":
                placeholder_mapping = "_"
            elif placeholder == "version":
                placeholder_mapping = f"(?P<{placeholder}>v{pattern})"
            elif placeholder == "extension":
                placeholder_mapping = f"(?P<{placeholder}>.{pattern})"
            elif placeholder == "start_date_file":
                placeholder_mapping = f"(?P<{placeholder}>{pattern})"
            else:
                placeholder_mapping = f"(?P<{placeholder}>{pattern})"
            final_pattern += placeholder_mapping

        # Add start and end anchors
        return f"^{final_pattern}$"

    #########################################################################################################

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

        self.supported_time_formats = {
            "YYMMDD": (r"%y%m%d"),
            "YYYYjjj": (r"%Y%j"),
            "YYYY_jjj": (r"%Y_%j"),
            "YYjjjhhmm": (r"%y%j%H%M"),
            "YYYY_JJJ_hhmm": (r"%Y_%j_%H%M"),
            "YYYY_MM_DD": (r"%Y_%m_%d"),
            "YYYY-MM-DD": (r"%Y-%m-%d"),
            "YYYY-jjjThhmmss": (r"%Y-%jT%H:%M:%S"),
            "YYYYMMDD_": (
                r"%Y%m%d_",
            ),  # added this for mex_sa_date_trick (formats %Y_%j_%H%M and %Y%m%d might be confused, so I add an artificial underscore.)
        }

        for key, str_format in self.supported_time_formats.items():
            try:
                self.date = datetime.strptime(date, str_format)
                if self.date is not None:
                    return self.date
            except (ValueError, TypeError):
                continue

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

        self.supported_time_formats = {
            "YYMMDD": (r"%y%m%d"),
            "YYYYjjj": (r"%Y%j"),
            "YYYY_jjj": (r"%Y_%j"),
            "YYjjjhhmm": (r"%y%j%H%M"),
            "YYYY_JJJ_hhmm": (r"%Y_%j_%H%M"),
            "YYYY_MM_DD": (r"%Y_%m_%d"),
            "YYYY-MM-DD": (r"%Y-%m-%d"),
            "YYYYMMDD_": (
                r"%Y%m%d_"
            ),  # added this for mex_sa_date_trick (formats %Y_%j_%H%M and %Y%m%d might be confused, so I add an artificial underscore.)
        }

        if format_key in self.supported_time_formats:
            # Get the format string and use strftime to format the datetime
            str_format = self.supported_time_formats[format_key]
            self.date = date.strftime(str_format)
            return self.date
        else:
            print("Format key not supported")
        #########################################################################################################

    def match_type_extension(self, data_type, filename):
        """
        Checks if the extension of a given file matches the expected extension for the specified data type.

        Args:
            data_type (str): The data type (e.g., 'ck', 'spk', 'odf', etc.).
            filename (str): The filename whose extension is to be checked.

        Returns:
            bool: True if the file extension matches any valid extension for the given data type.
        """

        # Get the file extension and compare it
        extension = filename.split(".")[-1].lower()
        valid_exts = self.type_to_extension.get(data_type.lower())

        if not valid_exts:
            return False

        if isinstance(valid_exts, str):
            valid_exts = [valid_exts]

        return extension in valid_exts

    #########################################################################################################

    def get_extension_for_data_type(self, data_type, first_only=True):
        """
        Returns one or more file extensions for the given data type.

        Args:
            data_type (str): The data type (e.g., 'ck', 'spk', etc.).
            first_only (bool): If True, returns only the first matching extension.
                            If False, returns a list of all valid extensions.

        Returns:
            str or list[str] or None: Extension(s) for the given data type, or None if unsupported.
        """
        ext = self.type_to_extension.get(data_type.lower())
        if not ext:
            return None

        if isinstance(ext, list):
            return ext[0] if first_only else ext
        return ext

    #########################################################################################################

    def parse_filename(self, input_mission, data_type, filename):
        """
        Description:
        This function parses a filename into its components based on the mission and data type patterns.
        It uses regular expressions to match the filename against predefined patterns for different mission types.
        It then returns the parsed components as a dictionary, along with the indices of underscores within the matched groups.

        Inputs:
            - input_mission (`str`): The mission name (e.g., 'cassini', 'mro', etc.).
            - data_type (`str`): The type of data (e.g., "radioscience", "kernels", "all").
            - filename (`str`): The filename to be parsed.

        Outputs:
            - (`dict`, `list`): A tuple containing:
              - `dictionary` (`dict`): A dictionary where the keys are the component names
              (e.g., 'date', 'purpose', etc.), and the values are the parsed values from the filename.
              - `underscore_indices` (`list`): A list of indices representing
              the positions of underscores in the matched groups.
        """

        input_mission = input_mission.lower()
        data_type_lower = data_type.lower()
        all_supported_types = list(self.type_to_extension.keys())

        if input_mission in self.supported_patterns:
            level_one_object = self.supported_patterns[input_mission]
        else:
            raise ValueError("Selected Mission Not Supported (yet!) Aborting ...")

        if data_type_lower == "all":
            data_type_types = all_supported_types

        elif data_type_lower in all_supported_types:
            data_type_types = [data_type_lower]

        else:
            print(f"Specified data type: {data_type} not supported.")

        # Initialize the dictionary and underscore index list
        underscore_indices = []

        if type(level_one_object) is not str:  # deals with multiple keys
            for data_type_type in level_one_object.keys():
                for data_type_type in data_type_types:
                    if self.match_type_extension(data_type_type, filename):
                        pattern = level_one_object[data_type_type]
                        match = re.match(pattern, filename)

                        if match:
                            # Populate the dictionary with matched groups
                            dictionary = match.groupdict()

                            # Add underscore index tracking based on group positions
                            last_pos = 0
                            group_index = 0
                            for key in dictionary.keys():
                                current_pos = match.start(
                                    key
                                )  # Get the start position of the matched group
                                if dictionary[key] is not None:
                                    if (
                                        last_pos != current_pos and last_pos != 0
                                    ):  # Check if there's a gap since last valid group
                                        underscore_indices.append(
                                            group_index
                                        )  # Track the index of the underscore
                                    last_pos = match.end(
                                        key
                                    )  # Move to the end of the current matched group
                                group_index += 1  # Increment the group index

                            # If present, convert start and end dates in utc (both will only be present in spice file names).
                            if "start_date_file" in dictionary and "end_date_file" in dictionary:
                                if dictionary["start_date_file"] is not None:

                                    if (
                                        dictionary["start_date_file"][0] == "P"
                                    ):  # this deals with MEX ORMF files

                                        self.start_date_utc = self.format_string_to_datetime(
                                            dictionary["start_date_file"][1:7],
                                        )
                                        self.end_date_utc = self.start_date_utc + relativedelta(
                                            months=1
                                        )
                                        dictionary["start_date_utc"] = self.start_date_utc
                                        dictionary["end_date_utc"] = self.end_date_utc

                                    elif (
                                        len(dictionary["start_date_file"]) == 4
                                        and dictionary["purpose"] == "SA"
                                    ):
                                        mex_sa_date_trick = dictionary["start_date_file"] + "0101_"
                                        self.start_date_utc = self.format_string_to_datetime(
                                            mex_sa_date_trick
                                        )
                                        self.end_date_utc = self.start_date_utc + relativedelta(
                                            years=1
                                        )
                                        dictionary["start_date_utc"] = self.start_date_utc
                                        dictionary["end_date_utc"] = self.end_date_utc

                                    else:
                                        self.start_date_utc = self.format_string_to_datetime(
                                            dictionary["start_date_file"],
                                        )
                                        self.end_date_utc = (
                                            self.format_string_to_datetime(
                                                dictionary["end_date_file"],
                                            )
                                            if dictionary["end_date_file"] != "000000"
                                            else self.start_date_utc + relativedelta(months=1)
                                        )  # this deals with MEX ORMM files
                                        dictionary["start_date_utc"] = self.start_date_utc
                                        dictionary["end_date_utc"] = self.end_date_utc

                            # If present, convert date_file in utc (only one date is present in the Radio Science file names)
                            elif "date_file" in dictionary:
                                self.date_utc = self.format_string_to_datetime(
                                    dictionary["date_file"]
                                )
                                dictionary["date_utc"] = self.date_utc

                            return dictionary, underscore_indices

                        else:
                            # print(f'Filename: {filename} does not match any supported patterns.')
                            continue
                    else:
                        continue

    #########################################################################################################

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

        # Extract the keys (group names) from the dictionary in order
        group_values = list(dictionary.values())

        # Initialize an empty list to build the filename
        reconstructed_filename = []

        # Iterate through the group values and append to the list
        for i, group in enumerate(group_values):
            if isinstance(group, str):
                reconstructed_filename.append(group)  # Add the current group

            # Check if the current index is in the underscore_indices list
            if (
                i + 1 in underscore_indices
            ):  # Underscore comes *after* the current group (between i and i+1)
                reconstructed_filename.append("_")

        # Join the list into a single string (the reconstructed filename)
        return "".join(reconstructed_filename)

    #########################################################################################################

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
        self.in_intervals = False
        for interval in intervals:
            if (date.date() >= interval[0].date()) and (date.date() <= interval[1].date()):
                self.in_intervals = True
        return self.in_intervals
