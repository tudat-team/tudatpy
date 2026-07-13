"""Cassini Titan-flyby file discovery and download."""

import requests
import os
from urllib.request import urlretrieve
from tabulate import tabulate
from colorama import Fore


class CassiniMixin:
    def print_titan_flyby_table(self):
        """
        Description:
        This method prints a table displaying the Titan flyby data in a readable format. It iterates over the `cassini_titan_flyby_dict`, extracts relevant information for each flyby, and constructs a table. The table is formatted with headers and colored using the `colorama` library for enhanced readability.

        Input:
            None.

        Output:
            None: This method prints the formatted table to the console.
        """

        data = []
        for key, value in self.cassini_titan_flyby_dict.items():
            row = [
                key,
                value["experiment"],
                value["pds_repo"],
                value["ancillary_repo"],
                value["cumindex_repo"],
                value["date"],
                value["doy"],
                value["year"],
            ]
            data.append(row)

        headers = [
            "Flyby ID",
            "Experiment",
            "PDS Repository",
            "Ancillary Repository",
            "Cumulative Index Repository",
            "Date",
            "DOY",
            "Year",
        ]

        # Adding color to the header
        colored_headers = [f"{Fore.MAGENTA}{header}{Fore.RESET}" for header in headers]

        print(tabulate(data, colored_headers, tablefmt="fancy_grid", stralign="center"))

    #########################################################################################################

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

        input_mission = "cassini"

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}

        flyby_ID = local_folder.split("/")[-1]

        if not os.path.exists(local_folder):
            os.mkdir(local_folder)

        filenames_to_download = []
        flyby_dict = self.cassini_flyby2experiment(flyby_ID)

        # Ensure 'cumindex_repo' is set to 'pds_repo' if 'cumindex_repo' is None
        cumindex_repo = (
            flyby_dict["cumindex_repo"].lower()
            if flyby_dict["cumindex_repo"] is not None
            else flyby_dict["pds_repo"].lower()
        )

        # Get URLs and cumulative index table
        volume_url = self.get_cassini_flyby_volume_url(flyby_dict)
        cumindex_url = self.get_cassini_flyby_cumindex_url(flyby_dict)
        cumindex_dict = self.get_cassini_volume_cumindex_table(cumindex_url)

        # Create dictionary of filenames by pds_repo
        filenames_dict = {pds_repo: data["file_name"] for pds_repo, data in cumindex_dict.items()}

        for pds_repo, filenames_list in filenames_dict.items():
            wanted_filenames = [
                volume_url + pds_repo.lower() + "/" + filename for filename in filenames_list
            ]
            filenames_to_download.extend(wanted_filenames)

        print("================================================================")
        print(
            f"Download {input_mission.upper()} Kernels (ck, spk) Ancillary Files (eop, ion, tro) and Radio Science (odf) files from PDS Atmosphere Node:"
        )
        # Download each file if not already present
        for filename in filenames_to_download:
            # Extract the actual filename from the URL
            actual_filename = filename.split("/")[-1]

            if actual_filename.lower().endswith(".ckf"):
                file_ext = "ck"
                file_folder = file_ext
                local_file_path = os.path.join(local_folder, file_folder, actual_filename)

            else:
                file_ext = actual_filename.lower().split(".")[-1]
                file_folder = file_ext
                local_file_path = os.path.join(local_folder, file_folder, actual_filename)

            # Check if the file already exists in the local folder
            if os.path.exists(local_file_path):
                print(
                    f"File: {actual_filename} already exists in {local_file_path} and will not be downloaded."
                )
                # Add file paths for 'ck' type files

                if file_ext.lower() in ["ion", "tro"]:
                    self.ancillary_files_to_load.setdefault(file_ext, []).append(local_file_path)

                elif file_ext.lower() in ["odf"]:
                    self.radio_science_files_to_load.setdefault(file_ext, []).append(
                        local_file_path
                    )
                else:
                    self.kernel_files_to_load.setdefault(file_ext, []).append(local_file_path)

            else:
                os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
                try:
                    # Download the file if it doesn't exist
                    urlretrieve(filename, local_file_path)
                    print(f"Downloading: '{filename}' to: {local_file_path}")

                    if file_ext.lower() in ["ion", "tro", "eop"]:
                        self.ancillary_files_to_load.setdefault(file_ext, []).append(
                            local_file_path
                        )

                    elif file_ext.lower() in ["odf"]:
                        self.radio_science_files_to_load.setdefault(file_ext, []).append(
                            local_file_path
                        )
                    else:
                        self.kernel_files_to_load.setdefault(file_ext, []).append(local_file_path)

                except Exception as e:
                    try:
                        urlretrieve(filename.lower(), local_file_path)
                        print(f"Downloading: '{filename.lower()}' to: {local_file_path}")

                        if file_ext.lower() in ["ion", "tro", "eop"]:
                            self.ancillary_files_to_load.setdefault(file_ext, []).append(
                                local_file_path
                            )

                        elif file_ext.lower() in ["odf"]:
                            self.radio_science_files_to_load.setdefault(file_ext, []).append(
                                local_file_path
                            )
                        else:
                            self.kernel_files_to_load.setdefault(file_ext, []).append(
                                local_file_path
                            )

                    except Exception as e_lower:
                        print(
                            f"Error downloading {filename}: {e}. "
                            f"Retry with {filename.lower()} also failed: {e_lower}"
                        )

        # Frame Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Frame Kernels from NAIF:")
        url_frame_files = "https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/fk/"
        wanted_frame_files = ["cas_v43.tf"]
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

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    ########################################################################################################################################

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

        # Check if the flyby_ID exists in the dictionary
        if flyby_ID in self.cassini_titan_flyby_dict:
            # Return the corresponding data for the given flyby_ID
            return self.cassini_titan_flyby_dict[flyby_ID]
        else:
            # If the flyby_ID is not found, return a message indicating so
            return f"Flyby ID {flyby_ID} not found."

    ########################################################################################################################################

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
        if moon.lower() == "titan":
            # Return the keys of the cassini_titan_flyby_dict
            return list(self.cassini_titan_flyby_dict.keys())
        elif moon.lower() == "enceladus":
            # Return the keys of the enceladus_flyby_dict
            return list(self.enceladus_flyby_dict.keys())
        else:
            raise ValueError("Invalid moon name. Please provide a valid Saturn Moon name.")

    ########################################################################################################################################

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

        # Extract the 'experiment' value from the flyby_dict
        experiment = flyby_dict["experiment"]

        # Define the URL template
        pds_repo_url = f"https://atmos.nmsu.edu/pdsd/archive/data/co-ssa-rss-1-{experiment}-v10/"

        # Check if the URL exists by sending a HEAD request
        try:
            response = requests.head(pds_repo_url, allow_redirects=True)

            if response.status_code == 200:
                # If the status code is 200, the URL exists
                self.pds_repo_url = pds_repo_url
                return self.pds_repo_url
            else:
                # If the status code is not 200, URL doesn't exist
                print(
                    f"URL not found: {pds_repo_url}. Please check the template URL in: get_cassini_flyby_volume_url"
                )
                return None

        except requests.exceptions.RequestException as e:
            # Catch any exceptions (e.g., network errors, invalid URL)
            print(f"Error checking URL: {e}")
            return None

    ########################################################################################################################################

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

        # Extract the 'experiment' value from the flyby_dict
        experiment = flyby_dict["experiment"]
        cumindex_repo = (
            flyby_dict["cumindex_repo"]
            if flyby_dict["cumindex_repo"] is not None
            else flyby_dict["pds_repo"]
        )

        # Define the URL template

        cumindex_url = f"https://atmos.nmsu.edu/pdsd/archive/data/co-ssa-rss-1-{experiment}-v10/{cumindex_repo}/INDEX/CUMINDEX.TAB"
        # Check if the URL exists by sending a HEAD request
        try:
            response = requests.head(cumindex_url, allow_redirects=True)  # Sends a HEAD request

            if response.status_code == 200:
                # If the status code is 200, the URL exists
                self.cumindex_url = cumindex_url
                return self.cumindex_url
            else:
                # If the status code is not 200, URL doesn't exist
                cumindex_url = f"https://atmos.nmsu.edu/pdsd/archive/data/co-ssa-rss-1-{experiment}-v10/{cumindex_repo}/index/cumindex.tab"
                response = requests.head(cumindex_url, allow_redirects=True)  # Sends a HEAD request
                if response.status_code == 200:
                    # If the status code is 200, the URL exists
                    self.cumindex_url = cumindex_url
                    return self.cumindex_url
                else:
                    print(
                        f"URL not found: {cumindex_url}. Please check the template URL in: get_cassini_flyby_cumindex_url"
                    )
                    return None  # Or return an appropriate message indicating the issue

        except requests.exceptions.RequestException as e:
            # Catch any exceptions (e.g., network errors, invalid URL)
            print(f"Error checking URL: {e}")
            return None  # Or return an error message

    ########################################################################################################################################

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

        # Fetch the .tab file content from the URL
        response = requests.get(cumindex_url)
        response.raise_for_status()  # Ensure the request was successful

        # Split the response content into lines
        lines = response.text.splitlines()

        # Prepare the dictionary to store the result
        cumindex_table = {}

        # Iterate over each line in the file
        for line in lines:
            # Split the line by commas (assuming the .tab file is comma-separated)
            cols = line.split(",")

            # Skip lines that don't have the expected number of columns
            if len(cols) != 7:
                continue

            # Extract the columns
            pds_repo = cols[0].replace('"', "").replace("'", "").strip()
            file_label = cols[1].replace('"', "").replace("'", "").strip()
            file_label_path = file_label.split("/")[0] + "/" + file_label.split("/")[1]

            file_name = file_label_path + "/" + cols[2].replace('"', "").replace("'", "").strip()
            external_file_name = cols[3].replace('"', "").replace("'", "").strip()
            start_date_utc = cols[4].replace('"', "").replace("'", "").strip()[:-2]
            end_date_utc = cols[5].replace('"', "").replace("'", "").strip()[:-2]
            creation_date_utc = cols[6].replace('"', "").replace("'", "").strip()

            keywords = [
                ("tigm", "odf"),
                ("tigf", "odf"),
                ("ancillary", "eop"),
                ("ancillary", "tro"),
                ("ancillary", "ion"),
                ("ancillary", "spk"),
                ("ancillary", "ckf"),
            ]

            # Check if both words in any pair appear in the file_label
            if any(all(word in file_label.lower() for word in pair) for pair in keywords):
                try:
                    start_date_utc = self.format_string_to_datetime(start_date_utc)
                    end_date_utc = self.format_string_to_datetime(end_date_utc)
                    creation_date_utc = self.format_string_to_datetime(creation_date_utc)
                except ValueError:
                    print("Skipping time conversion due to invalid date format.")
                    continue  # Skip rows with invalid date format

                # Use setdefault to ensure the key exists and initialize lists if not
                if pds_repo not in cumindex_table:
                    cumindex_table[pds_repo] = {
                        "file_label": [],
                        "file_name": [],
                        "external_file_name": [],
                        "start_date_utc": [],
                        "end_date_utc": [],
                        "creation_date_utc": [],
                    }

                # Append data to the lists for the current pds_repo
                cumindex_table[pds_repo]["file_label"].append(file_label)
                cumindex_table[pds_repo]["file_name"].append(file_name)
                cumindex_table[pds_repo]["external_file_name"].append(external_file_name)
                cumindex_table[pds_repo]["start_date_utc"].append(start_date_utc)
                cumindex_table[pds_repo]["end_date_utc"].append(end_date_utc)
                cumindex_table[pds_repo]["creation_date_utc"].append(creation_date_utc)

        # Return the cumindex_table
        return cumindex_table

    ########################################################################################################################################
    ########################################## END OF CASSINI SECTION #####################################################################
    ########################################################################################################################################
