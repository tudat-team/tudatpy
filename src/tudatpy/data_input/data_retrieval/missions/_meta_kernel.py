"""Meta-kernel discovery/download and the add_custom_mission_* extension
points for registering new missions at runtime."""

import requests
from bs4 import BeautifulSoup
import os
from urllib.request import urlretrieve
import re
from collections import defaultdict


class _MetaKernelMixin:
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
        input_mission = input_mission.lower()
        custom_dict = {}
        custom_dict[input_mission] = {custom_kernel_type: custom_pattern}
        self.supported_patterns.update(custom_dict)

        return self.supported_patterns

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

        input_mission = input_mission.lower()
        custom_dict = {}
        custom_dict[input_mission] = url
        self.supported_mission_meta_kernel_url.update(custom_dict)

        return self.supported_mission_meta_kernel_url

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
        input_mission = input_mission.lower()
        custom_dict = {}
        custom_dict[input_mission] = url
        self.supported_mission_kernels_url.update(custom_dict)

        return self.supported_mission_kernels_url

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

        input_mission = input_mission.lower()
        custom_dict = {}
        custom_dict[input_mission] = re.compile(f"{custom_pattern}")
        self.supported_mission_meta_kernel_pattern.update(custom_dict)

        return self.supported_mission_meta_kernel_pattern

    #########################################################################################################

    def extract_kernels_from_meta_kernel(self, input_mission):
        """
        Fetches a meta-kernel file from an HTTPS URL and categorizes kernel files by type
        using the type-to-extension mapping.

        Parameters:
            input_mission (str): e.g. 'mex', 'mro', 'juice', 'ro', etc...

        Returns:
            dict: A dictionary categorizing kernel files by type based on `type_to_extension`.
        """

        input_mission = input_mission.lower()
        self.kernel_files_names = defaultdict(list)

        meta_kernel_url = self.get_latest_meta_kernel(input_mission)

        try:
            # Fetch the meta-kernel file content
            response = requests.get(meta_kernel_url)
            response.raise_for_status()  # Raise an error for bad status codes
            lines = response.text.splitlines()

            for line in lines:
                line = line.strip()

                # Ignore comment lines and empty lines
                if line.startswith("\\") or not line:
                    continue

                # Extract the kernel path between quotes
                if "'" in line or '"' in line:
                    quote_char = "'" if "'" in line else '"'
                    start_idx = line.find(quote_char)
                    end_idx = line.rfind(quote_char)

                    if start_idx != -1 and end_idx != -1:
                        relative_kernel_path = line[start_idx + 1 : end_idx]
                        full_kernel_path = relative_kernel_path.replace(
                            "$KERNELS/",
                            self.supported_mission_kernels_url[input_mission],
                        )

                        # Extract the file extension
                        file_extension = full_kernel_path.split(".")[-1].lower()

                        # Match extension to kernel type
                        matched_key = next(
                            (
                                key
                                for key, exts in self.type_to_extension.items()
                                if file_extension in (exts if isinstance(exts, list) else [exts])
                            ),
                            None,  # No fallback, unmatched extensions will be ignored
                        )

                        # Add to results if matched
                        if matched_key:
                            self.kernel_files_names[matched_key].append(full_kernel_path)

            # Add the meta-kernel URL to the dictionary
            file_extension = meta_kernel_url.split(".")[-1].lower()

            # Match extension to kernel type
            matched_key = next(
                (
                    key
                    for key, exts in self.type_to_extension.items()
                    if file_extension in (exts if isinstance(exts, list) else [exts])
                ),
                None,  # No fallback, unmatched extensions will be ignored
            )
            # Add to results if matched
            if matched_key:
                self.kernel_files_names[matched_key].append(meta_kernel_url)

            return self.kernel_files_names

        except requests.exceptions.RequestException as e:
            print(f"An error occurred while fetching the meta-kernel: {e}")
            return {}
        except Exception as e:
            print(f"An error occurred in meta kernel extraction: {e}")
            return {}

    #########################################################################################################

    def download_kernels_from_meta_kernel(self, input_mission, local_folder):
        """Download kernels listed by the latest mission meta-kernel.

        Parameters
        ----------
        input_mission : str
            Mission identifier.
        local_folder : str
            Local output directory.

        Returns
        -------
        dict[str, list[str]]
            Kernel URLs grouped by kernel type.
        """

        input_mission = input_mission.lower()
        self.kernel_files_to_load = self.extract_kernels_from_meta_kernel(input_mission)

        skipped_files = []
        for kernel_type, kernel_urls in self.kernel_files_to_load.items():
            for kernel_url in kernel_urls:
                if (
                    not kernel_type == "mk"
                ):  # some meta-kernel paths are different from kernel paths
                    url_kernel_path = kernel_url[
                        len(self.supported_mission_kernels_url[input_mission]) :
                    ]
                    local_file_path = os.path.join(local_folder, url_kernel_path)
                else:
                    url_kernel_path = os.path.join(
                        kernel_type,
                        kernel_url[len(self.supported_mission_meta_kernel_url[input_mission]) :],
                    )
                    local_file_path = os.path.join(local_folder, url_kernel_path)

                local_kernel_folder = os.path.dirname(local_file_path)
                os.makedirs(local_kernel_folder, exist_ok=True)

                # Meta-kernels should always be re-downloaded to assure the latest version
                if kernel_type == "mk" or not os.path.exists(local_file_path):
                    action = "Re-downloading" if kernel_type == "mk" else "Downloading"
                    print(f"{action}: '{kernel_url}' to: {local_file_path}")
                    urlretrieve(kernel_url, local_file_path)
                    # Patch PATH_VALUES in meta-kernel files so they point to the local directory
                    if kernel_type == "mk":
                        self._patch_meta_kernel_path_values(local_file_path, local_folder)
                else:
                    skipped_files.append(url_kernel_path)

        if skipped_files:
            formatted = self._format_existing_files(skipped_files)
            print(
                f"The following meta-kernel files already exist in {local_folder} and will not be downloaded:\n\n{formatted}\n"
            )

        return self.kernel_files_to_load

    #########################################################################################################

    def _patch_meta_kernel_path_values(self, meta_kernel_path, local_folder):
        """
        Patches the PATH_VALUES in a downloaded meta-kernel (.TM) file to point
        to the local directory where kernels have been downloaded.

        Parameters:
            meta_kernel_path (str): Path to the meta-kernel file.
            local_folder (str): The local folder where kernels are stored.
        """
        try:
            with open(meta_kernel_path, "r") as f:
                content = f.read()

            # The mk/ directory is a subdirectory of local_folder, so PATH_VALUES
            # should point to local_folder (parent of mk/)
            resolved_path = os.path.abspath(local_folder)

            # Replace PATH_VALUES block: match the pattern PATH_VALUES = ( '...' )
            patched_content = re.sub(
                r"(PATH_VALUES\s*=\s*\(\s*')[^']*('\s*\))",
                rf"\g<1>{resolved_path}\2",
                content,
            )

            if patched_content != content:
                with open(meta_kernel_path, "w") as f:
                    f.write(patched_content)
                print(f"Patched PATH_VALUES in {meta_kernel_path} to: {resolved_path}")
            else:
                print(f"No PATH_VALUES found to patch in {meta_kernel_path}")

        except Exception as e:
            print(f"Warning: Could not patch meta-kernel PATH_VALUES: {e}")

    #########################################################################################################

    def get_latest_meta_kernel(self, input_mission):
        """
        Finds the most recent meta-kernel file URL based on year and version.

        Parameters:
            base_url (str): The base URL containing meta-kernel file links.

        Returns:
            str: The URL of the most recent meta-kernel file.
        """
        input_mission = input_mission.lower()
        base_url = self.supported_mission_meta_kernel_url[input_mission]
        meta_kernel_pattern = self.supported_mission_meta_kernel_pattern[
            input_mission
        ]  # the pattern for mex, ro, and juice is an exact name

        try:
            # Fetch the HTML page
            response = requests.get(base_url)
            response.raise_for_status()  # Raise an error for bad status codes
            soup = BeautifulSoup(response.text, "html.parser")

            # Extract all matching filenames and their year/version
            meta_kernels = []
            for link in soup.find_all("a", href=True):
                match = meta_kernel_pattern.match(link["href"])
                if match:
                    # Find the most recent meta-kernel based on year and version
                    try:
                        year = int(match.group(2))
                        version = int(match.group(3))
                        meta_kernels.append(
                            (year, version, base_url + link["href"])
                        )  # iterating through the past list of meta kernels
                    except (IndexError, ValueError):
                        self.latest_kernel = (
                            base_url + link["href"]
                        )  # there is just one exact name for juice and mex (no-brainer)
                        return self.latest_kernel

            if meta_kernels:
                self.latest_kernel = max(meta_kernels, key=lambda x: (x[0], x[1]))
                return self.latest_kernel[2]  # Return the URL of the most recent file
            else:
                print("No meta-kernels found matching the pattern.")
                return None

        except requests.exceptions.RequestException as e:
            print(f"An error occurred while fetching the meta-kernel list: {e}")
            return None
        except Exception as e:
            print(f"An error occurred in getting latest meta-kernel: {e}")
            return None

    #########################################################################################################

    def get_latest_clock_kernel_name(self, input_mission):
        """Return clock-kernel file names from the latest meta-kernel.

        Parameters
        ----------
        input_mission : str
            Mission identifier.

        Returns
        -------
        list[str]
            Clock-kernel file names.
        """

        input_mission = input_mission.lower()
        kernels_from_meta_kernel = self.extract_kernels_from_meta_kernel(input_mission)

        clock_files = kernels_from_meta_kernel.get("sclk", [])

        # The shared GRAIL meta-kernel contains the clock kernels for both
        # spacecraft.  Select the kernel belonging to the requested spacecraft
        # instead of treating the expected pair as an ambiguity.
        grail_clock_prefixes = {"grail-a": "gra_", "grail-b": "grb_"}
        if input_mission in grail_clock_prefixes:
            clock_prefix = grail_clock_prefixes[input_mission]
            clock_files = [
                clock_file
                for clock_file in clock_files
                if clock_file.rsplit("/", 1)[-1].lower().startswith(clock_prefix)
            ]

        return [clock_file.rsplit("/", 1)[-1] for clock_file in clock_files]
