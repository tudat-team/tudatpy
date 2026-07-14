"""Composed LoadPDS entry point: wires the generic download/parsing engine
together with each mission's file-discovery mixin, plus the top-level
get_mission_files dispatcher."""

import copy
import os
from datetime import timedelta
import re
from tudatpy.data_input.environment_data import spice

from ._patterns import SUPPORTED_PATTERNS, CASSINI_TITAN_FLYBY_DICT
from ._download import _DownloadEngineMixin
from ._parsing import _ParsingMixin
from ._meta_kernel import _MetaKernelMixin
from .mex import MexMixin
from .juice import JuiceMixin
from .mro import MroMixin
from .cassini import CassiniMixin
from .grail import GrailMixin
from .ro import RoMixin


class LoadPDS(
    _DownloadEngineMixin,
    _ParsingMixin,
    _MetaKernelMixin,
    MexMixin,
    JuiceMixin,
    MroMixin,
    CassiniMixin,
    GrailMixin,
    RoMixin,
):
    """Mission data downloader for PDS and SPICE resources."""

    def __init__(self):
        self.flag_check_existing_files = False
        self.flag_load_standard_kernels = False

        self.supported_patterns = copy.deepcopy(SUPPORTED_PATTERNS)
        self.cassini_titan_flyby_dict = copy.deepcopy(CASSINI_TITAN_FLYBY_DICT)

        # Mapping of data types to their expected file extensions
        self.type_to_extension = {
            "ck": ["bc"],
            "dsk": ["bds"],
            "spk": ["bsp"],
            "fk": ["tf"],
            "mk": ["tm"],
            "ik": ["ti"],
            "lsk": ["tls"],
            "pck": ["bpc", "tpc"],  # supports both PCK types
            "sclk": ["tsc"],
            "odf": ["odf"],
            "tnf": ["tnf"],
            "ifms": ["tab"],
            "dp2": ["tab"],
            "dps": ["tab"],
            "dpx": ["tab"],
            "ion": ["ion"],
            "tro": ["tro"],
            "eop": ["eop"],
            "maneuver": ["asc"],
            "antenna_switch": ["asc"],
        }
        # Mapping of data types to their expected file extensions
        self.supported_mission_odf_time_formats = {
            "mro": "YYYY_jjj",
            "grail-a": "YYYY_jjj",
            "grail-b": "YYYY_jjj",
        }

        self.supported_mission_meta_kernel_url = {
            "mro": "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/extras/mk/",
            "mex": "https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/mk/",
            "juice": "https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/mk/",
            "cassini": "https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/extras/mk/",
            "grail-a": "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/extras/mk/",
            "grail-b": "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/extras/mk/",
            "ro": "https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/mk/",
        }

        self.supported_mission_meta_kernel_pattern = {
            "mro": re.compile(r"^(mro)_(\d{4})_v(\d{2})\.tm$"),
            "mex": re.compile(r"MEX_OPS.TM"),
            "juice": re.compile(r"juice_ops.tm"),
            "cassini": re.compile(r"^cas_(\d{4})_v(\d{2})\.tm$"),
            "grail-a": re.compile(r"^grail_v(\d{2})\.tm$"),
            "grail-b": re.compile(r"^grail_v(\d{2})\.tm$"),
            "ro": re.compile(r"ROS_OPS.TM"),
        }

        self.supported_mission_kernels_url = {
            "mro": "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/",
            "mex": "https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/",
            "juice": "https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/",
            "cassini": "https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/data/",
            "grail-a": "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/",
            "grail-b": "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/",
            "ro": "https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/",
        }

    #########################################################################################################

    def get_mission_files(
        self,
        input_mission,
        start_date=None,
        end_date=None,
        flyby_IDs=None,
        custom_output=None,
        all_meta_kernel_files=None,
        load_kernels=None,
        radio_observation_type=None,
        radio_science_file_type="odf",
    ):
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
        input_mission = input_mission.lower()
        self.all_kernel_files = {}
        self.all_radio_science_files = {}
        self.all_ancillary_files = {}

        if start_date and end_date:
            all_dates = [
                start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)
            ]
        print(
            f"======================================= Downloading {input_mission.upper()} Data ==============================================\n"
        )

        if custom_output:
            if (
                input_mission == "grail-a" or input_mission == "grail-b"
            ):  # both grail-a and -b found in grail archive
                base_folder = f"{custom_output}/{input_mission}"
            else:
                base_folder = f"{custom_output}"

        else:
            if (
                input_mission == "grail-a" or input_mission == "grail-b"
            ):  # both grail-a and -b found in grail archive
                base_folder = f"grail_archive/{input_mission}"
            else:
                base_folder = f"{input_mission}_archive/"

        local_folder_list = (
            []
        )  # this will be a multiple-element list (length != 1)ONLY if "Cassini" and ONLY if len(flyby_IDs) > 1
        if input_mission == "cassini":
            if not flyby_IDs:
                self.print_titan_flyby_table()
                flyby_IDs = input(
                    f"You are about to download data for the Cassini Spacecraft. What flyby would you like to download? (Check Provided Table for Reference.)\n"
                )
                # Split the input string by commas and remove any leading/trailing spaces from each item
                flyby_IDs = [flyby.strip() for flyby in flyby_IDs.split(",")]
                # Output the list of flyby IDs
                print(f"You selected the following flyby_IDs: {flyby_IDs}.")

            if isinstance(flyby_IDs, list):
                # Create a flag to track if Titan or Enceladus has been processed
                processed_moons = []

                # Process Titan flybys if 'ALL_TITAN' is in the list
                if "ALL_TITAN" in flyby_IDs:
                    print("removing")
                    flyby_IDs.remove("ALL_TITAN")  # Remove 'ALL_TITAN'
                    full_moon_flybys_list = self.get_cassini_full_moon_flybys_list("TITAN")
                    flyby_IDs.extend(full_moon_flybys_list)
                    # Remove duplicates
                    flyby_IDs = list(dict.fromkeys(flyby_IDs))
                    # Add Titan folder creation to the list
                    processed_moons.append("TITAN")

                # Process Enceladus flybys if 'ALL_ENCELADUS' is in the list
                if "ALL_ENCELADUS" in flyby_IDs:
                    flyby_IDs.remove("ALL_ENCELADUS")  # Remove 'ALL_ENCELADUS'
                    full_moon_flybys_list = self.get_cassini_full_moon_flybys_list("ENCELADUS")
                    flyby_IDs.extend(full_moon_flybys_list)
                    # Remove duplicates
                    flyby_IDs = list(dict.fromkeys(flyby_IDs))
                    # Add Enceladus folder creation to the list
                    processed_moons.append("ENCELADUS")

                # At this point, `flyby_IDs` contains the full list of flyby IDs without 'ALL_TITAN' or 'ALL_ENCELADUS'
                # Now create the local folders

                if len(processed_moons) != 0:
                    for moon in processed_moons:
                        for flyby_ID in flyby_IDs:
                            local_folder = os.path.join(base_folder, moon, flyby_ID)
                            local_folder_list.append(local_folder)  # Append to local_folder_list
                else:
                    for flyby_ID in flyby_IDs:
                        if flyby_ID.startswith("T"):
                            moon = "TITAN"
                        elif flyby_ID.startswith("E"):
                            moon = "ENCELADUS"
                        local_folder = os.path.join(base_folder, moon, flyby_ID)
                        local_folder_list.append(local_folder)  # Append to local_folder_list

            else:
                for moon in ["TITAN", "ENCELADUS"]:
                    if f"ALL_{moon}" == flyby_IDs:
                        flyby_IDs.remove(f"ALL_{moon}")
                        full_moon_flybys_list = self.get_cassini_full_moon_flybys_list(moon)
                        flyby_IDs.extend(full_moon_flybys_list)
                        flyby_IDs = list(set(flyby_IDs))  # This removes duplicates

                        for flyby_ID in flyby_IDs:
                            local_folder = os.path.join(base_folder, moon, flyby_ID)
                            local_folder_list.append(local_folder)  # append to local_folder_list

                if flyby_IDs.startswith("T"):
                    moon = "TITAN"
                elif flyby_IDs.startswith("E"):
                    moon = "ENCELADUS"

                local_folder = os.path.join(base_folder, moon, flyby_IDs)
                local_folder_list.append(local_folder)  # in this case, it is just a single folder
        else:
            local_folder_list.append(base_folder)  # in this case, it is just a single folder

        print(
            f"=========================================== Folder(s) Creation =================================================="
        )

        if os.path.exists(base_folder) == False:
            print(f"Creating Local Folder: {base_folder} and its subfolders (kernels and radio)...")
            os.makedirs(base_folder)
            for local_folder in local_folder_list:
                if os.path.exists(local_folder) == False:
                    os.makedirs(local_folder)
                else:
                    print(f"Folder: {local_folder} already exists and will not be overwritten.")
        else:
            print(f"Folder: {base_folder} already exists and will not be overwritten.")

        print(
            f"===========================================================================================\n"
        )

        kernel_files_to_load = None
        for local_folder in local_folder_list:
            self.flag_check_existing_files = False
            if all_meta_kernel_files:
                print(f"Downloading all Kernels from {input_mission} Latest Meta-Kernel File...")
                kernel_files_to_load = self.download_kernels_from_meta_kernel(
                    input_mission, local_folder
                )

            if input_mission == "ro":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = self.get_ro_files(
                        local_folder,
                        start_date,
                        end_date,
                        radio_observation_type,
                        skip_kernel_downloads=True,
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_ro_files(
                        local_folder,
                        start_date,
                        end_date,
                        radio_observation_type,
                    )

            if input_mission == "mex":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = self.get_mex_files(
                        local_folder,
                        start_date,
                        end_date,
                        radio_observation_type,
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_mex_files(
                        local_folder,
                        start_date,
                        end_date,
                        radio_observation_type,
                    )

            elif input_mission == "juice":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = self.get_juice_files(
                        local_folder, start_date, end_date
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_juice_files(local_folder, start_date, end_date)

            elif input_mission == "mro":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = self.get_mro_files(
                        local_folder,
                        start_date,
                        end_date,
                        radio_science_file_type=radio_science_file_type,
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_mro_files(
                        local_folder,
                        start_date,
                        end_date,
                        radio_science_file_type=radio_science_file_type,
                    )
            elif input_mission == "cassini":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = (
                        self.get_cassini_flyby_files(local_folder)
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_cassini_flyby_files(local_folder)

            elif input_mission == "grail-a":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = (
                        self.get_grail_a_files(local_folder, start_date, end_date)
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_grail_a_files(local_folder, start_date, end_date)
            elif input_mission == "grail-b":
                if kernel_files_to_load:
                    _, radio_science_files_to_load, ancillary_files_to_load = (
                        self.get_grail_b_files(local_folder, start_date, end_date)
                    )

                else:
                    (
                        kernel_files_to_load,
                        radio_science_files_to_load,
                        ancillary_files_to_load,
                    ) = self.get_grail_b_files(local_folder, start_date, end_date)

        if kernel_files_to_load:
            meta_kernel_present = "mk" in kernel_files_to_load

            # Always populate self.all_kernel_files with all converted kernels
            for kernel_type, kernel_files in kernel_files_to_load.items():
                self.all_kernel_files[kernel_type] = []
                for kernel_file in kernel_files:
                    converted_kernel_file = self.spice_transfer2binary(kernel_file)
                    self.all_kernel_files[kernel_type].append(converted_kernel_file)

            # If meta-kernel is present and load_kernels is True: load only meta-kernel
            if meta_kernel_present and load_kernels:
                try:
                    meta_kernel_file = self.all_kernel_files["mk"][0]
                    spice.load_kernel(meta_kernel_file)
                except Exception as e:
                    print(f"Failed to load meta-kernel file: {meta_kernel_file}, Error: {e}")

            # If meta-kernel is not present: load individual kernels if requested
            elif not meta_kernel_present and load_kernels:
                for kernel_type, kernel_files in self.all_kernel_files.items():
                    for converted_kernel_file in kernel_files:
                        try:
                            spice.load_kernel(converted_kernel_file)
                        except Exception as e:
                            print(
                                f"!! Failed to load kernel: {converted_kernel_file}, Error: {e} !!"
                            )

        else:
            print("No Kernel Files to Load.")

        if ancillary_files_to_load:
            for (
                ancillary_type,
                ancillary_files,
            ) in ancillary_files_to_load.items():
                for ancillary_file in ancillary_files:  # Iterate over each file in the list
                    try:
                        if ancillary_type not in self.all_ancillary_files.keys():
                            self.all_ancillary_files[ancillary_type] = [ancillary_file]
                            if load_kernels:
                                spice.load_kernel(ancillary_file)
                        else:
                            self.all_ancillary_files[ancillary_type].append(ancillary_file)
                            if load_kernels:
                                spice.load_kernel(ancillary_file)

                    except Exception as e:
                        print(f"!! Failed to load kernel: {ancillary_file}, Error: {e} !!")
        else:
            print("No Ancillary Files to Load.")

        if radio_science_files_to_load:
            for (
                radio_science_type,
                radio_science_files,
            ) in radio_science_files_to_load.items():
                for radio_science_file in radio_science_files:  # Iterate over each file in the list
                    if radio_science_type not in self.all_radio_science_files.keys():
                        self.all_radio_science_files[radio_science_type] = [radio_science_file]
                    else:
                        self.all_radio_science_files[radio_science_type].append(radio_science_file)
        else:
            print("No Radio Science Files to Load.")

        n_kernels = spice.get_total_count_of_kernels_loaded()
        if not self.flag_load_standard_kernels:
            if load_kernels:
                print("================================================================")
                print(f"Number of Loaded Existing + Downloaded Kernels: {n_kernels}")
                std_kernels = spice.load_standard_kernels()
                self.flag_load_standard_kernels = True
                n_standard_kernels = spice.get_total_count_of_kernels_loaded() - n_kernels
                print(f"Number of Loaded Standard Kernels: {n_standard_kernels}")
                print("================================================================")
        else:
            print(f"Number of Loaded Existing + Downloaded + Standard Kernels: {n_kernels}")
            print("================================================================")

        self.clean_mission_archive(local_folder)

        return (
            self.all_kernel_files,
            self.all_radio_science_files,
            self.all_ancillary_files,
        )
