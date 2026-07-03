"""Generic download engine shared by every mission (HTTP fetch, versioned-file
lookup, SPACIT kernel conversion, date-range crawling)."""

import requests
from bs4 import BeautifulSoup
import os
from urllib.request import urlretrieve
from datetime import timedelta
import glob
import re
import subprocess


_REQUEST_TIMEOUT = 30  # seconds


class _DownloadEngineMixin:
    def spice_transfer2binary(self, input_file, timeout=5):
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

        if input_file.lower().endswith(".ckf"):
            output_file = input_file.split(".")[0] + ".ck"

            if not os.path.exists(output_file):
                proc = subprocess.Popen(
                    ["spacit"],
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

                # Write the necessary inputs to the process' stdin stream
                proc.stdin.write("T\n")  # Command to convert transfer files to binary
                proc.stdin.write(f"{input_file}\n")  # Input file
                proc.stdin.write(f"{output_file}\n")  # Output file

                try:
                    # Wait for the process to finish, with a timeout
                    stdout, stderr = proc.communicate(timeout=timeout)

                    # Check if the process was successful
                    if proc.returncode != 0:
                        print(f"Error converting {input_file} to binary. Error: {stderr}")
                    else:
                        print(
                            f"Successfully converted {input_file} to {output_file}. Output: {stdout}"
                        )

                except subprocess.TimeoutExpired:
                    proc.kill()  # Kill the process if it times out

            return output_file

        elif input_file.lower().endswith(".spk"):
            output_file = input_file.split(".")[0] + ".bsp"

            if not os.path.exists(output_file):
                proc = subprocess.Popen(
                    ["spacit"],
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

                # Write the necessary inputs to the process' stdin stream
                proc.stdin.write("T\n")  # Command to convert transfer files to binary
                proc.stdin.write(f"{input_file}\n")  # Input file
                proc.stdin.write(f"{output_file}\n")  # Output file

                try:
                    # Wait for the process to finish, with a timeout
                    stdout, stderr = proc.communicate(timeout=timeout)

                    # Check if the process was successful
                    if proc.returncode != 0:
                        print(f"Error converting {input_file} to binary. Error: {stderr}")
                    else:
                        print(
                            f"Successfully converted {input_file} to {output_file}. Output: {stdout}"
                        )

                except subprocess.TimeoutExpired:
                    proc.kill()  # Kill the process if it times out

            return output_file

        else:
            output_file = input_file
            return output_file

    def _download_file(self, url, local_path, timeout=60):
        """Download a file via requests with timeout and atomic write.

        Downloads to a temporary file first, then atomically renames to the
        final path on success.  On failure the temp file is cleaned up.
        """
        tmp_path = local_path + ".tmp"
        try:
            response = requests.get(url, timeout=timeout, stream=True)
            response.raise_for_status()
            with open(tmp_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            os.replace(tmp_path, local_path)
        except Exception:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            raise

    #########################################################################################################

    def _get_latest_versioned_files(self, url, patterns):
        """Fetch an HTML listing page and find the latest version per pattern.

        Parameters:
            url (str): URL of the directory listing page.
            patterns (list[str]): Regex patterns each containing one capture
                group for the version number, e.g. ``r"ROS_V(\\d+)\\.TF"``.

        Returns:
            list[str]: One filename per pattern — the one with the highest
            version number.  Patterns with no match are silently skipped.
        """
        try:
            response = requests.get(url, timeout=_REQUEST_TIMEOUT)
            response.raise_for_status()
        except Exception as e:
            print(
                f"Warning: Failed to fetch versioned file listing from {url}: {e}. "
                f"Returning empty list."
            )
            return []
        html = response.text

        latest_files = []
        for pattern in patterns:
            best_version = -1
            best_filename = None
            for m in re.finditer(pattern, html, re.IGNORECASE):
                version_num = int(m.group(1))
                if version_num > best_version:
                    best_version = version_num
                    best_filename = m.group(0)
            if best_filename is not None:
                latest_files.append(best_filename)
            else:
                print(f"Warning: No file matched pattern {pattern!r} at {url}")
        return latest_files

    #########################################################################################################

    def get_kernels(
        self,
        input_mission,
        url,
        wanted_files=None,
        wanted_files_patterns=None,
        custom_output=None,
    ):
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

        input_mission = input_mission.lower()
        # Determine data type and file extension
        data_type = url.split("/")[-2].lower()
        ext = self.get_extension_for_data_type(data_type)

        if not ext and wanted_files:
            ext = wanted_files[0].split(".")[-1]
            if ext == "asc":
                data_type = "maneuver"
            elif ext == "level_1b":
                data_type = "antenna_switch"
            elif ext == "level_0":
                data_type = "maneuver"

        if custom_output:
            if input_mission in {
                "grail-a",
                "grail-b",
            }:  # both grail-a and -b found in grail archive
                base_folder = f"{custom_output}/{input_mission}"
            else:
                base_folder = f"{custom_output}"

        else:
            if input_mission in {
                "grail-a",
                "grail-b",
            }:  # both grail-a and -b found in grail archive
                base_folder = f"grail_archive/{input_mission}"
            else:
                base_folder = f"{input_mission}_archive/"

        os.makedirs(os.path.join(base_folder, data_type), exist_ok=True)

        # Retrieve files matching pattern if specified
        matched_files_list = []
        if wanted_files_patterns:
            # Parse the URL to find matching files
            reqs = requests.get(url, timeout=_REQUEST_TIMEOUT)
            soup = BeautifulSoup(reqs.text, "html.parser")

            for wanted_files_pattern in wanted_files_patterns:
                # Extract all links that match the pattern
                regex_pattern = re.escape(wanted_files_pattern).replace(r"\*", ".*")
                matched_files = [
                    os.path.basename(link.get("href"))
                    for link in soup.find_all("a", href=True)
                    if re.match(regex_pattern, link.get("href"))
                ]
                matched_files_list.extend(matched_files)

        # Combine explicitly specified files and pattern-matched files
        all_files_to_download = set(wanted_files or []) | set(matched_files_list)

        self.files_to_load = [
            os.path.join(base_folder, data_type, file) for file in all_files_to_download
        ]

        # Download each file if not already present
        skipped_files = []
        for file, local_file in zip(all_files_to_download, self.files_to_load):
            if not os.path.exists(local_file):
                full_url = str(os.path.join(url, file))
                print(f"Downloading File: {full_url} to: {local_file}")
                self._download_file(full_url, local_file)
            else:
                skipped_files.append(os.path.basename(local_file))

        if skipped_files:
            formatted = self._format_existing_files(
                skipped_files  # already basenames, _format_existing_files handles both
            )
            print(f"The following files already exist and will not be downloaded:\n\n{formatted}\n")

        return self.files_to_load

    #########################################################################################################

    def clean_mission_archive(self, local_folder):
        """
        Description:
            Cleans up the mission archive by removing empty subdirectories within the specified local folder.

        Input:
            - `local_folder` (`str`): Path to the local folder where the archive is stored.

        Output:
            - `None`: This function does not return any value. It performs in-place directory cleanup.
        """

        print(f"Cleaning Mission Archive...")
        if os.path.exists(local_folder):
            for directory in os.listdir(local_folder):
                directory_path = os.path.join(local_folder, directory)
                if os.path.isdir(directory_path) and not os.listdir(directory_path):
                    os.rmdir(directory_path)
        print(f"Done.")

    #########################################################################################################

    def dynamic_download_url_files_single_time(
        self,
        input_mission,
        local_path,
        start_date,
        end_date,
        url,
        verbose=True,
    ):
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
            - verbose (`bool`): If True, print progress messages. If False, operate quietly.

        Outputs:
            - `self.relevant_files` (`list`): A list of paths to the files that were either already present or successfully downloaded.
        """

        input_mission = input_mission.lower()
        self.relevant_files = []
        self._last_skipped_files = []
        self._last_downloaded_files = []
        if verbose:
            print("Checking URL:", url)
        data_type = url.split("/")[-2]
        local_subfolder = os.path.join(local_path, data_type.lower())

        try:
            supported_pattern = self.supported_patterns[input_mission][data_type.lower()]

        except KeyError:
            raise ValueError("Pattern not found among supported patterns.")

        existing_files = self.check_existing_files(data_type, local_subfolder)
        if not existing_files:
            existing_files = set()

        all_dates = [
            start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)
        ]
        # Initialize a dictionary to hold files from the HTML response
        files_url_dict = {}

        try:
            reqs = requests.get(url, timeout=_REQUEST_TIMEOUT)
            reqs.raise_for_status()
        except Exception as e:
            raise ValueError(f"Error fetching data from {url}: {e}")

        # Parse links from the HTML response
        for link in BeautifulSoup(reqs.text, "html.parser").find_all("a"):
            if not "[To Parent Directory]" in link:
                full_link = link.get("href")  # Extract the href attribute
                if isinstance(full_link, str):
                    # Check if the link is a full URL or a relative path
                    if ("/") in full_link:
                        # It's a full URL
                        filename = os.path.basename(full_link)  # Get the filename from the URL
                    else:
                        # It's a relative path
                        filename = full_link

                    # Check if the filename matches the specified format
                    if self.match_type_extension(data_type, filename):
                        try:
                            RS_dict, RS_underscores = self.parse_filename(
                                input_mission, data_type, filename
                            )
                            filename_to_download = self.reconstruct_filename(
                                RS_dict, RS_underscores
                            )
                            # Determine the date string format
                            date_key = (
                                RS_dict["date_file"][:5]
                                if input_mission in ("mex", "ro")
                                else RS_dict["date_file"][:8]
                            )
                            files_url_dict[date_key] = filename_to_download

                            # Downloading only the latest version of a file
                            version = RS_dict.get(
                                "version"
                            )  # Define default version if not provided
                            if not version:
                                files_url_dict[(start_date, end_date)] = filename_to_download

                            else:
                                ext = RS_dict.get("extension")
                                # Extract the base filename without the version
                                base_name_no_version_no_ext = filename_to_download.replace(
                                    version, ""
                                ).replace(ext, "")
                                current_version = int(
                                    version[1:]
                                )  # Extract numeric version (e.g., v02 -> 2)
                                # Only store the highest version
                                if not any(
                                    base_name_no_version_no_ext in value
                                    for value in files_url_dict.values()
                                ):
                                    files_url_dict[date_key] = filename_to_download
                                else:
                                    stored_filename = files_url_dict[date_key]
                                    stored_version_str = stored_filename.replace(
                                        base_name_no_version_no_ext, ""
                                    ).replace(ext, "")
                                    stored_version = int(stored_version_str[1:])
                                    if current_version >= stored_version:
                                        latest_filename_to_download = (
                                            f"{base_name_no_version_no_ext}{version}{ext}"
                                        )
                                        files_url_dict[date_key] = latest_filename_to_download

                        except Exception as e:
                            if input_mission != "grail-a" and input_mission != "grail-b":
                                print(f"Could not parse file: {filename} - Error: {e}")
                                continue
                    else:
                        continue
            else:
                continue

        # Download missing files and collect existing ones that match the date range
        for date in all_dates:
            if input_mission in ("mex", "ro"):
                date_string = f"{date.year % 100:02d}{date.timetuple().tm_yday:03d}"

            else:
                if input_mission in self.supported_mission_odf_time_formats.keys():
                    format_key = self.supported_mission_odf_time_formats[input_mission]
                    date_string = self.format_datetime_to_string(date, format_key)
                else:
                    print(
                        f"No ODF time format associated to input mission: {input_mission}. Please provide it in self.supported_mission_odf_time_formats."
                    )

            if date_string in files_url_dict:
                download_file = files_url_dict[date_string]
                full_download_url = os.path.join(url, download_file)

                full_local_path = os.path.join(local_subfolder, download_file)
                os.makedirs(local_subfolder, exist_ok=True)
                if full_local_path in existing_files:
                    # File already exists locally and matches date range — include it
                    if full_local_path not in self.relevant_files:
                        self.relevant_files.append(full_local_path)
                        self._last_skipped_files.append(full_local_path)
                        if verbose:
                            print(
                                f"File already exists, skipping: {os.path.basename(download_file)}"
                            )
                else:
                    try:
                        if verbose:
                            print(f"Downloading: {full_download_url} to: {full_local_path}")
                        self._download_file(
                            full_download_url,
                            full_local_path,
                        )
                        self.relevant_files.append(full_local_path)
                        self._last_downloaded_files.append(full_local_path)
                    except Exception as e:
                        print(f"!! Failed to download {full_download_url}: {e} !!")

        if len(self.relevant_files) == 0:
            if verbose:
                print("Nothing to download.")

        if verbose:
            print("...Done.")

        return self.relevant_files

    #########################################################################################################

    def check_existing_files(self, data_type, local_subfolder):
        """
        Description:
            Checks the local directory for files that match a given pattern.
            This function filters and returns the files that already exist locally, based on their extensions and matching patterns.

        Inputs:
            - data_type (`str`): The type of data (e.g., 'ck', 'spk').
            - local_subfolder (`str`): Path to the local directory where the files are stored.

        Outputs:
            - `self.existing_files` (`set`): A set of paths to the files that already exist and match the given pattern.
        """
        # Get all existing files that match the filename format
        ext = self.get_extension_for_data_type(data_type)
        self.existing_files = {
            f
            for f in glob.glob(f"{local_subfolder}/*")
            if re.search(rf"\.{ext}$", f, re.IGNORECASE)
        }
        if self.existing_files:
            self.flag_check_existing_files = True
        return self.existing_files

    #########################################################################################################

    @staticmethod

    def _format_existing_files(existing_files):
        """Format a set of file paths as a sorted bulleted list of basenames."""
        sorted_names = sorted(os.path.basename(f) for f in existing_files)
        return "\n".join(f"  - {name}" for name in sorted_names)

    #########################################################################################################

    def dynamic_download_url_files_time_interval(
        self, input_mission, local_path, start_date, end_date, url
    ):
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
        input_mission = input_mission.lower()
        self.relevant_files = []
        print("Checking URL:", url)
        data_type = url.split("/")[-2]

        if not os.path.exists(local_path):
            print(f"Creating Local Folder: {local_path}")
            os.makedirs(local_path, exist_ok=True)

        local_subfolder = os.path.join(local_path, data_type)

        if not os.path.exists(local_subfolder):
            print(f"Creating Local Subfolder: {local_subfolder}")
            os.makedirs(local_subfolder, exist_ok=True)

        # Prepare the date range for searching existing files
        all_dates = [
            start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)
        ]

        try:
            supported_pattern = self.supported_patterns[input_mission][data_type]
        except KeyError:
            raise ValueError(f"Pattern not found among supported patterns.")

        existing_files = self.check_existing_files(data_type, local_subfolder)
        if not existing_files:
            existing_files = set()

        if existing_files:
            print(
                f"--------------------------------------- EXISTING FILES CHECK ---------------------------------------------\n"
            )
            print(
                f"The following files already exist in the folder and will not be downloaded:\n\n{self._format_existing_files(existing_files)}\n"
            )

        # Initialize a dictionary to hold files from the HTML response
        files_url_dict = {}

        # Get the content of the URL
        reqs = requests.get(url, timeout=_REQUEST_TIMEOUT)

        # Parse links from the HTML response
        for link in BeautifulSoup(reqs.text, "html.parser").find_all("a"):
            full_link = link.get("href")  # Extract the href attribute
            if isinstance(full_link, str):
                # Check if the link is a full URL or a relative path
                if ("/") in full_link:
                    # It's a full URL
                    filename = os.path.basename(full_link)  # Get the filename from the URL
                else:
                    # It's a relative path
                    filename = full_link
                    # Check if the filename matches the specified format
                if self.match_type_extension(data_type, filename):
                    try:
                        # Parse the information from the filename
                        dictionary, underscores = self.parse_filename(
                            input_mission, data_type, filename
                        )
                        # Reconstruct the filename and get time intervals
                        filename_to_download = self.reconstruct_filename(dictionary, underscores)
                        start_time = dictionary["start_date_utc"]
                        end_time = dictionary["end_date_utc"]

                        # Downloading only the latest version of a file
                        version = dictionary.get(
                            "version"
                        )  # Define default version if not provided

                        if not version:
                            files_url_dict[(start_time, end_time)] = filename_to_download

                        else:
                            ext = dictionary.get("extension")

                            if start_time is None or end_time is None:
                                print(
                                    f"Unwanted filename found: {filename_to_download}. Skipping... [Do not worry! ;)]"
                                )
                                continue

                            # Extract the base filename without the version
                            base_name_no_version_no_ext = filename_to_download.replace(
                                version, ""
                            ).replace(ext, "")
                            current_version = int(
                                version[1:]
                            )  # Extract numeric version (e.g., v02 -> 2)
                            # Only store the highest version
                            if not any(
                                base_name_no_version_no_ext in value
                                for value in files_url_dict.values()
                            ):
                                files_url_dict[(start_time, end_time)] = filename_to_download
                            else:
                                stored_filename = files_url_dict[(start_time, end_time)]
                                stored_version_str = stored_filename.replace(
                                    base_name_no_version_no_ext, ""
                                ).replace(ext, "")
                                stored_version = int(stored_version_str[1:])
                                if current_version >= stored_version:
                                    latest_filename_to_download = (
                                        f"{base_name_no_version_no_ext}{version}{ext}"
                                    )
                                    files_url_dict[(start_time, end_time)] = (
                                        latest_filename_to_download
                                    )

                    except Exception:
                        continue  # Skip to the next link

        # Download files for all intervals from the HTML response
        # Use direct interval overlap test: [new_start, new_end] overlaps [start_date, end_date]
        # iff new_start <= end_date AND new_end >= start_date
        for new_interval, filename_to_download in files_url_dict.items():
            new_start, new_end = new_interval
            # Check if the file's interval overlaps with the requested date range
            if not (new_start.date() <= end_date.date() and new_end.date() >= start_date.date()):
                continue

            full_local_path = os.path.join(local_subfolder, filename_to_download)

            if full_local_path in existing_files:
                # File already exists locally and overlaps date range — include it
                if full_local_path not in self.relevant_files:
                    self.relevant_files.append(full_local_path)
            else:
                print(
                    f"Downloading: {os.path.join(url,filename_to_download)} to: {full_local_path}"
                )  # Print which file is being downloaded
                self._download_file(
                    os.path.join(url, filename_to_download), full_local_path
                )  # Download the file
                self.relevant_files.append(full_local_path)

        if len(self.relevant_files) == 0:
            print("Nothing to download.")

        print("...Done.")  # Indicate completion

        return self.relevant_files

    def download_url_files_time(
        self,
        local_path,
        filename_format,
        start_date,
        end_date,
        url,
        time_format,
        indices_date_filename,
    ):
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

        os.makedirs(local_path, exist_ok=True)

        # Retrieve all dates contained within the time interval defined by the start and end dates provided as inputs
        all_dates = [
            start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)
        ]

        # Split the file name at the symbol "/", to separate between the name of the file and the folder in which it might be contained
        filename_split = filename_format.split("/")
        folder = ""
        if len(filename_split) > 2:
            raise Exception(
                "In the current implementation, the filename format cannot contain more than one folder."
            )
        elif len(filename_split) == 2:
            folder = filename_split[0]

        # In the reduced file name (after removing the folder part), replace the wildcards - indicated by "\w" - by "*", which will later allow the
        # BeautifulSoup package to look for all pattern-matching names at the targeted url (without any a priori information on the date and/or
        # wildcard present in the file name)
        reduced_filename = filename_split[-1]
        reduced_filename = reduced_filename.replace(r"\w", "*")

        # Retrieve all filenames present at the "local_path" location that match the specified filename format

        existing_files = glob.glob(local_path + reduced_filename)
        if len(existing_files) > 0:
            print(
                f"The following files already exist in the folder:\n\n {existing_files}\n\n and will not be downloaded."
            )

        # Identify dates of interest that are not covered by existing files
        self.relevant_files = []
        dates_without_file = []
        for date in all_dates:
            # Create string object corresponding to the date under investigation
            date_str = date.strftime(time_format)

            # Reconstruct the expected filename for this particular date
            current_filename = ""
            index = 0
            for ind in indices_date_filename:
                current_filename += filename_format[index:ind] + date_str
                index = ind + 1
            current_filename += filename_format[index:]

            # Parse existing files and check whether the current file name can be found
            current_files = [
                x
                for x in existing_files
                if re.match(local_path + current_filename.split("/")[-1], x)
            ]
            # If so, add the identified file to the list of relevant files to be loaded
            if len(current_files) > 0:
                for file in current_files:
                    self.relevant_files.append(current_files[0])
            # If not, mark the current date as non-covered by any file yet (i.e., date for which a file is missing)
            else:
                dates_without_file.append(date)

        # Retrieve the missing files from the specified url
        if len(dates_without_file) > 0:

            # List containing the names of the files to be downloaded
            files_url = []

            # Parse all files contained at the targeted url
            reqs = requests.get(url)
            for link in BeautifulSoup(reqs.text, "html.parser").find_all("a"):

                # Retrieve full url link for each of these files
                full_link = link.get("href")

                # Check whether the file of interest is nested within an additional folder
                if len(folder) == 0:
                    current_filename = full_link.split("/")[-1]
                else:
                    current_filename = full_link.split("/")[-2]

                # Store the name of the file to be downloaded
                files_url.append(current_filename)

        # Parse all dates for which a file was originally missing and download missing files from the specified url.
        for date in dates_without_file:

            # List of the files to be downloaded
            files_to_download = []

            # Reconstruct expected filename for the date under consideration
            date_str = date.strftime(time_format)
            current_filename = ""
            index = 0
            for ind in indices_date_filename:
                current_filename += filename_format[index:ind] + date_str
                index = ind + 1
            current_filename += filename_format[index:]

            # Check whether a matching file was found at the targeted url for this specific date, and split the filename at "/"
            # to account for the possibility that the targeted file is stored in a nested folder
            file_to_download = [x for x in files_url if re.match(current_filename.split("/")[0], x)]

            # If the file is directly stored at the specified url (no nested folder), then the filename can be stored directly
            if len(folder) == 0:
                files_to_download = file_to_download

            # Otherwise, explore additional folder layer
            if len(folder) > 0 and len(file_to_download) > 0:
                reqs2 = requests.get(url + file_to_download[0])

                # Parse all files within the current folder
                for nested_link in BeautifulSoup(reqs2.text, "html.parser").find_all("a"):
                    nested_file = nested_link.get("href")

                    # Retrieve all matching file names within the current folder
                    relevant_link = [
                        x
                        for x in [nested_file]
                        if re.match(current_filename.split("/")[-1], x.split("/")[-1])
                    ]

                    # If a match is found, store the filename that should be downloaded (now including the extra folder layer)
                    if len(relevant_link) == 1:
                        files_to_download.append(
                            file_to_download[0] + "/" + relevant_link[0].split("/")[-1]
                        )

            # Download all relevant files from the targeted url
            for file in files_to_download:
                print("Downloading ", url + file)
                urlretrieve(url + file, local_path + file.split("/")[-1])
                self.relevant_files.append(local_path + file.split("/")[-1])

        # Return the list of all relevant files that should be loaded to cover the time interval of interest
        return self.relevant_files

