import trk234
from . import converters as cnv
from pandas import concat as pd_concat
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.astro import time_representation
from tudatpy.dynamics.environment import (
    PiecewiseLinearFrequencyInterpolator,
)


class Trk234Processor:
    """
    Processor for TNF files using pytrk234.

    For a given set of requested observables types (e.g. ['doppler', 'range']),
    this processor iterates file-by-file, uses each converter’s extract method to obtain per-file data,
    merges the outputs, and then calls each converter’s process method to produce final observation
    collection. If simulation bodies are provided, ramp data are processed and used to set the
    stations' frequency interpolator with the set_tnf_information_in_bodies() method.

    Examples
    --------
    >>> from tudatpy.data import Trk234Processor
    >>>
    >>> # Define TNF file paths
    >>> tnf_files = ["mro_kernels/mromagr2012_002_1426xmmmv1.tnf"]
    >>>
    >>> # Create processor for both Doppler and range data
    >>> tnf_processor = Trk234Processor(
    ...     tnf_files,
    ...     ["doppler", "range"],
    ...     spacecraft_name="MRO"
    ... )
    >>>
    >>> # Process observations
    >>> observations = tnf_processor.process()
    >>>
    >>> # Set frequency information in the bodies assuming you have a bodies object tudatpy.dynamics.environment.SystemOfBodies
    >>> tnf_processor.set_tnf_information_in_bodies(bodies)
    """

    def __init__(
        self, tnf_file_paths, requested_types, spacecraft_name=None, bodies=None
    ):
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
        bodies : object, optional
            The simulation bodies container (used for setting ramp data).
        """
        self.tnf_file_paths = tnf_file_paths
        self.spacecraft_name = spacecraft_name
        self.bodies = bodies

        # Initialize observables converters.
        self.converters = {}
        if "doppler" in requested_types:
            self.converters["doppler"] = cnv.DerivedDopplerConverter()
        if "range" in requested_types:
            self.converters["range"] = cnv.DerivedSraRangeConverter()

        # Initialize ramp converter if needed.
        self.ramp_converter = cnv.RampConverter()

    def process(self):
        """
        Process all TNF files provided at initialization. For each file, decode the SFDU data,
        and for each requested radiometric data type, extract data via the converter's extract method.
        Ramp data is also extracted separately if a ramp converter is available.
        Then, merge the per-file outputs and process them to produce the final outputs.

        Returns
        -------
        ObservationCollection or None
            An ObservationCollection containing all radiometric observation sets,
            or None if none were processed.
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
                    observation_sets.extend(
                        converter.process(merged_df, self.spacecraft_name)
                    )

        return ObservationCollection(observation_sets)

    def set_tnf_information_in_bodies(self, bodies):
        """
        Update stations in bodies by setting the frequency interpolators from the ramp data.
        Set the transponder turnaround ratio for the spacecraft.

        NOTE: It's not optimal to set the transponder turnaround ratio here, but it's done for now

        Parameters
        ----------
        ramp_df : pd.DataFrame
            DataFrame containing ramp data.
        spacecraft_name : str
            The spacecraft name used for setting the transponder turnaround ratio.
        bodies : object
            The simulation bodies container.
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

        ramp_df["start_time_seconds"] = ramp_df["start_time"].apply(
            lambda x: time_representation.DateTime.from_python_datetime(x).to_epoch()
        )
        ramp_df["end_time_seconds"] = ramp_df["end_time"].apply(
            lambda x: time_representation.DateTime.from_python_datetime(x).to_epoch()
        )
        earth = bodies.get("Earth")
        for station in ramp_df["station"].unique():
            station_df = ramp_df[ramp_df["station"] == station]
            frequencyInterpolator = PiecewiseLinearFrequencyInterpolator(
                station_df["start_time_seconds"].tolist(),
                station_df["end_time_seconds"].tolist(),
                station_df["rate"].tolist(),
                station_df["freq"].tolist(),
            )
            groundStation = earth.get_ground_station(station)
            groundStation.set_transmitting_frequency_calculator(frequencyInterpolator)

        if self.spacecraft_name:
            spacecraft = bodies.get(self.spacecraft_name)
            spacecraft.system_models.set_default_transponder_turnaround_ratio_function()
