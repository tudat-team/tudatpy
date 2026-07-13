"""JUpiter ICy-moons Explorer (JUICE) file discovery and download."""


class JuiceMixin:
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

        self.radio_science_files_to_load = {}  # empty for now, since we need fdets
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}  # empty for now

        input_mission = "juice"
        # Clock files
        print("================================================================")
        print(f"Download {input_mission.upper()} Clock Files:")
        url_clock_files = "https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/sclk/"
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
        print(f"Download {input_mission.upper()} Frame Files:")
        url_frame_files = "https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/"
        wanted_frame_files = [
            "juice_v41.tf",
            "juice_events_crema_5_1_150lb_23_1_v02.tf",
            "juice_ops_v11.tf",
            "juice_stations_topo_v01.tf",
            "juice_dsk_surfaces_v11.tf",
            "juice_sci_v17.tf",
            "juice_roi_v02.tf",
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

        print("================================================================")
        print(f"Download {input_mission.upper()} Orientation Kernels:")
        ck_files_to_load = []
        wanted_ck_files = [
            "juice_mga_crema_5_1_150lb_23_1_baseline_v04.bc",
            "juice_sa_crema_5_1_150lb_23_1_baseline_v04.bc",
            "juice_sc_crema_5_1_150lb_23_1_baseline_v03.bc",
        ]

        url_planned_ck_files = "https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/"
        planned_ck_files_to_load = self.get_kernels(
            input_mission=input_mission,
            url=url_planned_ck_files,
            wanted_files=wanted_ck_files,
            custom_output=local_folder,
        )

        if planned_ck_files_to_load:
            for file in planned_ck_files_to_load:
                ck_files_to_load.append(file)
        else:
            print("No Crema Planned CK files to download this time.")

        measured_url_ck_files = ["https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/"]
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
            for file in measured_ck_files_to_load:
                ck_files_to_load.append(file)
        else:
            print("No Measured ck files to download this time.")

        if len(ck_files_to_load) > 0:
            self.kernel_files_to_load["ck"] = ck_files_to_load
        else:
            print("No Overall CK files to download this time.")

            # SPK files
        spk_files_to_load = []
        print("================================================================")
        print(f"Download {input_mission.upper()} SPK Kernels:")
        measured_url_spk_files = ["https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/"]

        if len(measured_url_spk_files) == 1:
            measured_spk_files_to_load = self.dynamic_download_url_files_time_interval(
                input_mission,
                local_path=local_folder,
                start_date=start_date,
                end_date=end_date,
                url=measured_url_spk_files[0],
            )
        else:
            for measured_url_spk_file in measured_url_spk_files:
                measured_spk_files_to_load = self.dynamic_download_url_files_time_interval(
                    input_mission,
                    local_path=local_folder,
                    start_date=start_date,
                    end_date=end_date,
                    url=measured_url_spk_file,
                )

        if measured_spk_files_to_load:
            for file in measured_spk_files_to_load:
                spk_files_to_load.append(file)
        else:
            print("No spk files to download this time.")

        print("---------------------------------------------")

        if len(spk_files_to_load) > 0:
            self.kernel_files_to_load["spk"] = spk_files_to_load
        else:
            print("No Overall SPK files to download this time.")

            # Tropospheric corrections
        # print(f'=======================================================================================\n')
        # print(f'Download {input_mission.upper()} tropospheric corrections files\n')
        # url_tropo_files = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/tro/"
        # tropo_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
        #    local_path=local_folder, start_date=start_date,
        #    end_date=end_date, url=url_tropo_files)

        # if tropo_files_to_load:
        #    self.ancillary_files_to_load['tro'] = tropo_files_to_load
        # else:
        #    print('No tropospheric files to download this time.')

        # Ionospheric corrections
        # print(f'=======================================================================================\n')
        # print(f'Download {input_mission.upper()} ionospheric corrections files\n')
        # url_ion_files = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/ion/"
        # ion_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
        #    local_path=local_folder, start_date=start_date, end_date=end_date,
        #    url=url_ion_files)

        # if ion_files_to_load:
        #    self.ancillary_files_to_load['ion'] = ion_files_to_load
        # else:
        #    print('No ionospheric files to download this time.')

        return (
            self.kernel_files_to_load,
            self.radio_science_files_to_load,
            self.ancillary_files_to_load,
        )

    ########################################################################################################################################
    ################################################### END OF JUICE SECTION ###############################################################
    ########################################################################################################################################
