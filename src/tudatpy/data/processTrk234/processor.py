import trk234
from . import converters as cnv
from pandas import concat as pd_concat
from tudatpy.estimation.observations import (
    ObservationCollection,
    ObservationDataset,
    create_observation_collection_from_dataset,
)
from tudatpy.astro import time_representation
from tudatpy.dynamics.environment import (
    PiecewiseLinearFrequencyInterpolator,
    SystemOfBodies,
    FrequencyGapHandling,
)
from .converters.ramp import OpenRampHandling


class Trk234Processor:
    """
    Processor for TNF files using pytrk234.

    For a given set of requested observables types (e.g. ['doppler', 'range']),
    this processor iterates file-by-file, uses each converter’s extract method to obtain per-file data,
    merges the outputs, and then calls each converter’s process method to produce final observation
    dataset. If simulation bodies are provided, ramp data are processed and used to set the
    stations' frequency interpolator with the set_tnf_information_in_bodies() method.

    Examples
    --------

    .. code-block:: python

        from tudatpy.data import Trk234Processor

        # Define TNF file paths
        tnf_files = ["mro_kernels/mromagr2012_002_1426xmmmv1.tnf"]

        # Create processor for both Doppler and range data
        tnf_processor = Trk234Processor(
            tnf_files,
            ["doppler", "range"],
            spacecraft_name="MRO"
        )

        # Process observations
        observation_dataset = tnf_processor.process()

        # Set frequency information in the bodies assuming you have a bodies object tudatpy.dynamics.environment.SystemOfBodies
        tnf_processor.set_tnf_information_in_bodies(bodies)

    """

    def __init__(
        self,
        tnf_file_paths: list[str],
        requested_types: list[str],
        spacecraft_name: str | None = None,
    ) -> None:
        """
        Parameters
        ----------
        tnf_file_paths : list[str]
            List of TNF file paths to be processed.
        requested_types : list[str]
            List of requested radiometric data types, e.g., ['doppler', 'range'].
            Note: "ramp" should NOT be included here.
        spacecraft_name : str, optional
            The spacecraft name for building link definitions.
        """
        self.tnf_file_paths = tnf_file_paths
        self.spacecraft_name = spacecraft_name

        # Initialize observables converters.
        self.converters: dict[str, cnv.Converter] = {}
        if "doppler" in requested_types:
            self.converters["doppler"] = cnv.DerivedDopplerConverter()
        if "range" in requested_types:
            self.converters["range"] = cnv.DerivedSraRangeConverter()

        # Initialize ramp converter if needed.
        self.ramp_converter = cnv.RampConverter()

    def process(self) -> ObservationDataset:
        """
        Process all TNF files provided at initialization. For each file, decode the SFDU data,
        and for each requested radiometric data type, extract data via the converter's extract method.
        Ramp data is also extracted separately if a ramp converter is available.
        Then, merge the per-file outputs and process them to produce the final outputs.

        Returns
        -------
        ObservationDataset
            Dataset containing all radiometric observation sets. If no observations were extracted,
            an empty dataset is returned.
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
        observation_dataset = ObservationDataset()
        for dtype, converter in self.converters.items():
            if extracted_data[dtype]:
                merged_df = pd_concat(extracted_data[dtype], ignore_index=True)
                if not merged_df.empty:
                    converter_dataset = converter.process(merged_df, self.spacecraft_name)
                    for set_id in range(converter_dataset.number_of_observation_sets):
                        observation_dataset.add_observation_set_from_dataset(
                            converter_dataset, set_id
                        )

        return observation_dataset

    def process_observation_collection(self) -> ObservationCollection:
        """
        Process all TNF files and return the result as a legacy ObservationCollection.

        Returns
        -------
        ObservationCollection
            Backwards-compatible collection facade created from the processed dataset.
        """
        return create_observation_collection_from_dataset(self.process())

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
