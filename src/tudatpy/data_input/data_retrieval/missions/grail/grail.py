"""GRAIL-A / GRAIL-B file discovery and download."""

import os


class GrailMixin:
    """GRAIL mission file-discovery mixin."""

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

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}

        input_mission = "grail-a"

        # ODF files
        print("================================================================")
        print(f"Download {input_mission.upper()} ODF files:")
        url_radio_science_files = [
            "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/grail_0201/odf/"
        ]
        for url_radio_science_file in url_radio_science_files:
            try:
                files = self.dynamic_download_url_files_single_time(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=url_radio_science_file,
                )
                key = "odf"
                self.radio_science_files_to_load[key] = files
            except Exception as e:
                print(f"Error downloading grail-a odf file from {url_radio_science_file}: {e}")
                continue

        if not self.radio_science_files_to_load:
            print("No Radio Science files to download this time.")

        # Clock Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Clock Kernels:")
        url_clock_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/sclk/"
        )
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

        # Frame Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Frame Kernels:")
        url_frame_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/fk/"
        )
        wanted_frame_files = ["grail_v07.tf"]
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

            # maneuver Files
        print("================================================================")
        print(f"Download {input_mission.upper()} Maneuver Files:")
        url_man_files = "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-2-edr-v1/grail_0001/level_0/2012_04_06/"
        wanted_man_files = ["mas00_2012_04_06_a_04.asc"]
        man_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_man_files,
            wanted_files=wanted_man_files,
            custom_output=local_folder,
        )

        if man_files_to_load:
            self.ancillary_files_to_load["maneuver"] = man_files_to_load
        else:
            print("No maneuver files to download this time.")

            # Antenna Switch Files
        print("================================================================")
        print(f"Download {input_mission.upper()} Antenna Switch Files:")

        url_antenna_switch_files = (
            "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-3-cdr-v1/grail_0101/level_1b/"
        )
        antenna_switch_files_to_load = self.download_url_files_time(
            local_path=os.path.join(local_folder, "antenna_switch/"),
            filename_format="*/vgs1b_*_a_04.asc",
            start_date=start_date,
            end_date=end_date,
            url=url_antenna_switch_files,
            time_format="%Y_%m_%d",
            indices_date_filename=[0, 8],
        )

        if antenna_switch_files_to_load:
            self.ancillary_files_to_load["antenna_switch"] = antenna_switch_files_to_load
        else:
            print("No Antenna Switch files to download this time.")

            # Planetary and Ephemeris Kernels
        spk_files_to_load = []
        print("================================================================")
        print(f"Download {input_mission.upper()} SPK Kernels:")
        url_spk_files = [
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/spk/"
        ]
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
            self.kernel_files_to_load["ck"] = spk_files_to_load
        else:
            print("No spk files to download this time.")

            # Orientation Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Orientation Kernels:")
        measured_url_ck_files = [
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/ck/"
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
            "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/grail_0201/ancillary/tro/"
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
            "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/grail_0201/ancillary/ion/"
        )
        ion_files_to_load = self.dynamic_download_url_files_time_interval(
            input_mission,
            local_path=local_folder,
            start_date=start_date,
            end_date=end_date,
            url=url_ion_files,
        )

        if ion_files_to_load:
            self.ancillary_files_to_load["ion"] = ion_files_to_load
        else:
            print("No ionospheric files to download this time.")

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    ########################################################################################################################################
    ########################################## END OF GRAIL-A SECTION ####################################################################
    ########################################################################################################################################

    ########################################################################################################################################
    ########################################## START OF GRAIL-B SECTION ####################################################################
    ########################################################################################################################################

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

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}

        input_mission = "grail-b"

        # ODF files
        print("================================================================")
        print(f"Download {input_mission.upper()} ODF files:")
        url_radio_science_files = [
            "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/grail_0201/odf/"
        ]
        for url_radio_science_file in url_radio_science_files:
            try:
                files = self.dynamic_download_url_files_single_time(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=url_radio_science_file,
                )
                key = "odf"
                self.radio_science_files_to_load[key] = files
            except Exception as e:
                print(f"Error downloading grail-b odf file from {url_radio_science_file}: {e}")
                continue

        if not self.radio_science_files_to_load:
            print("No Radio Science files to download this time.")

        # Clock Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Clock Kernels:")
        url_clock_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/sclk/"
        )
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

        # Frame Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Frame Kernels:")
        url_frame_files = (
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/fk/"
        )
        wanted_frame_files = ["grail_v07.tf"]
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

            # maneuver Files
        print("================================================================")
        print(f"Download {input_mission.upper()} Maneuver Files:")
        url_man_files = (
            "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-2-edr-v1/grail_0001/level_0/"
        )
        wanted_man_files = ["mas00_2012_04_06_a_04.asc"]
        man_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_man_files,
            wanted_files_patterns=wanted_man_files,
            custom_output=local_folder,
        )

        if man_files_to_load:
            self.ancillary_files_to_load["maneuver"] = man_files_to_load
        else:
            print("No maneuver files to download this time.")

        url_antenna_switch_files = (
            "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-3-cdr-v1/grail_0101/level_1b/"
        )
        wanted_antenna_switch_files_patterns = ["*/vgs1b_*_a_04.asc"]
        antenna_switch_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_antenna_switch_files,
            wanted_files_patterns=wanted_antenna_switch_files_patterns,
            custom_output=local_folder,
        )

        if antenna_switch_files_to_load:
            self.ancillary_files_to_load["antenna_switch"] = antenna_switch_files_to_load
        else:
            print("No Antenna Switch files to download this time.")

            # Planetary and Ephemeris Kernels
        spk_files_to_load = []
        print("================================================================")
        print(f"Download {input_mission.upper()} SPK Kernels:")
        url_spk_files = [
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/spk/"
        ]
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
            self.kernel_files_to_load["ck"] = spk_files_to_load
        else:
            print("No spk files to download this time.")

            # Orientation Kernels
        print("================================================================")
        print(f"Download {input_mission.upper()} Orientation Kernels:")
        measured_url_ck_files = [
            "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/grlsp_1000/data/ck/"
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
            "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/grail_0201/ancillary/tro/"
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
            "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/grail_0201/ancillary/ion/"
        )
        ion_files_to_load = self.dynamic_download_url_files_time_interval(
            input_mission,
            local_path=local_folder,
            start_date=start_date,
            end_date=end_date,
            url=url_ion_files,
        )

        if ion_files_to_load:
            self.ancillary_files_to_load["ion"] = ion_files_to_load
        else:
            print("No ionospheric files to download this time.")

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    ########################################################################################################################################
    ########################################## END OF GRAIL-B SECTION ####################################################################
    ########################################################################################################################################
