"""Legacy TNF processor with the original collection return and environment setup."""

import trk234
from pandas import concat as pd_concat
from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data.tnf import TnfTrackingDataProcessor, OpenRampHandling
from tudatpy.dynamics.environment import (
    PiecewiseLinearFrequencyInterpolator,
    SystemOfBodies,
    FrequencyGapHandling,
)
from tudatpy.estimation.observations import (
    ObservationCollection,
    create_observation_collection_from_tracking_data,
)


class Trk234Processor(TnfTrackingDataProcessor):
    def process(self) -> ObservationCollection:
        """
        Process all TNF files provided at initialization. For each file, decode the SFDU data,
        and for each requested radiometric data type, extract data via the converter's extract method.
        Ramp data is also extracted separately if a ramp converter is available.
        Then, merge the per-file outputs and process them to produce the final outputs.

        Returns
        -------
        ObservationCollection
            An ObservationCollection containing all radiometric observation sets, if none were extracted, an empty collection is returned.
        """
        # Accumulate outputs for radiometric converters.
        extracted_data = {key: [] for key in self.converters.keys()}

        for file_path in self.tnf_file_paths:
            reader = trk234.Reader(file_path)
            reader.decode(sec_chdo=False, trk_chdo=False)
            sfdu_list = reader.sfdu_list

            # Extract data.
            for dtype, converter in self.converters.items():
                extracted = converter.extract(sfdu_list)
                if not extracted.empty:
                    extracted_data[dtype].append(extracted)

        # Process observables data: merge extracted DataFrames and process them.
        observation_sets = []
        for dtype, converter in self.converters.items():
            if extracted_data[dtype]:
                merged_df = pd_concat(extracted_data[dtype], ignore_index=True)
                if not merged_df.empty:
                    observation_sets.extend(converter.process(merged_df, self.spacecraft_name))

        return create_observation_collection_from_tracking_data(observation_sets, SystemOfBodies())

    def set_tnf_information_in_bodies(
        self,
        bodies: SystemOfBodies,
        gap_handling: FrequencyGapHandling = FrequencyGapHandling.extrapolate_at_gaps,
        open_ramp_handling: OpenRampHandling = OpenRampHandling.print_warning_once,
    ) -> None:
        """
        Update stations in bodies by setting the frequency interpolators from the ramp data.
        Set the transponder turnaround ratio for the spacecraft.

        NOTE: It's not optimal to set the transponder turnaround ratio here, but it's done for now

        Parameters
        ----------
        bodies : SystemOfBodies
            The simulation bodies container.
        gap_handling : FrequencyGapHandling
            The gap handling strategy for the frequency interpolators. Defaults to ``FrequencyGapHandling.extrapolate_at_gaps``.
        open_ramp_handling : OpenRampHandling
            Strategy for closing open-ended ramp intervals. Defaults to ``OpenRampHandling.print_warning_once``.
        """

        ramp_data_list = []
        for file_path in self.tnf_file_paths:
            reader = trk234.Reader(file_path)
            reader.decode(sec_chdo=False, trk_chdo=False)
            sfdu_list = reader.sfdu_list
            ramp_extracted = self.ramp_converter.extract(sfdu_list)
            ramp_data_list.append(ramp_extracted)

        all_ramps = pd_concat(ramp_data_list, ignore_index=True)
        all_ramps.sort_values("epoch", inplace=True)
        all_ramps.reset_index(drop=True, inplace=True)

        ramp_df = self.ramp_converter.process(all_ramps)
        ramp_df = self.ramp_converter.handle_open_ramps(ramp_df, open_ramp_handling)

        ramp_df["start_time_seconds"] = ramp_df["start_time"].apply(
            lambda x: time_representation.DateTime.from_python_datetime(x).to_epoch()
        )
        ramp_df["end_time_seconds"] = ramp_df["end_time"].apply(
            lambda x: time_representation.DateTime.from_python_datetime(x).to_epoch()
        )
        earth = bodies.get("Earth")
        for station in ramp_df["station"].unique():
            station_df = ramp_df[ramp_df["station"] == station]
            frequency_interpolator = PiecewiseLinearFrequencyInterpolator(
                station_df["start_time_seconds"].tolist(),
                station_df["end_time_seconds"].tolist(),
                station_df["rate"].tolist(),
                station_df["freq"].tolist(),
                gap_handling=gap_handling,
            )
            ground_station = earth.get_ground_station(station)

            if not ground_station.has_frequency_calculator():
                ground_station.set_transmitting_frequency_calculator(frequency_interpolator)
            else:
                ground_station.transmitting_frequency_calculator.add_frequency_interpolator(
                    frequency_interpolator
                )

        if self.spacecraft_name:
            spacecraft = bodies.get(self.spacecraft_name)
            spacecraft.system_models.set_default_transponder_turnaround_ratio_function()
