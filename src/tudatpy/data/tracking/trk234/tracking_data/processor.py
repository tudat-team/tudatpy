import trk234
from . import converters as cnv
from pandas import concat as pd_concat
from tudatpy.kernel.data.tracking import TrackingData
from tudatpy.kernel.data.environment import (
    TrackingSupplementaryData,
    RampedFrequencySupplementaryData,
)
from tudatpy.astro import time_representation
from tudatpy.dynamics.environment import SystemOfBodies
from .converters.ramp import OpenRampHandling


class Trk234TrackingDataProcessor:
    """
    Processor for TNF files using pytrk234.

    For a given set of requested observables types (e.g. ['doppler', 'range']),
    this processor iterates file-by-file, uses each converter's extract method to obtain per-file
    data, merges the outputs, and then calls each converter's process method to produce
    :class:`~tudatpy.data.TrackingData` objects. Ramp data are always extracted and merged into
    :class:`~tudatpy.data.TrackingSupplementaryData` objects (one per ground station), holding the
    stations' frequency ramps.

    Examples
    --------

    .. code-block:: python

        from tudatpy.data import Trk234TrackingDataProcessor
        from tudatpy.estimation.observations import (
            create_observation_collection,
            set_tracking_supplementary_data_in_bodies,
        )

        # Define TNF file paths
        tnf_files = ["mro_kernels/mromagr2012_002_1426xmmmv1.tnf"]

        # Create processor for both Doppler and range data
        tnf_processor = Trk234TrackingDataProcessor(
            tnf_files,
            ["doppler", "range"],
            spacecraft_name="MRO"
        )

        # Process files into tracking data and supplementary (ramp) data
        tracking_data, supplementary_data = tnf_processor.process()

        # Convert the tracking data into an ObservationCollection, and set the frequency ramps in
        # the bodies, assuming you have a bodies object tudatpy.dynamics.environment.SystemOfBodies
        observations = create_observation_collection(tracking_data, bodies)
        set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
        tnf_processor.set_transponder_turnaround_ratio(bodies)

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

        # Initialize ramp converter.
        self.ramp_converter = cnv.RampConverter()

    def process(
        self,
        open_ramp_handling: OpenRampHandling = OpenRampHandling.print_warning_once,
    ) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:
        """
        Process all TNF files provided at initialization. For each file, decode the SFDU data,
        and for each requested radiometric data type, extract data via the converter's extract method.
        Ramp data is also extracted, to build the stations' frequency ramp supplementary data.
        The per-file outputs are then merged and processed to produce the final outputs.

        Parameters
        ----------
        open_ramp_handling : OpenRampHandling
            Strategy for closing open-ended ramp intervals. Defaults to ``OpenRampHandling.print_warning_once``.

        Returns
        -------
        tuple[list[TrackingData], list[TrackingSupplementaryData]]
            A tuple ``(tracking_data, supplementary_data)``. ``tracking_data`` contains one entry per
            requested radiometric observable group; ``supplementary_data`` contains one entry per
            ground station with associated frequency ramps. Either list is empty if no matching data
            were extracted.
        """
        # Accumulate outputs for radiometric converters.
        extracted_data = {key: [] for key in self.converters.keys()}
        ramp_data_list = []

        for file_path in self.tnf_file_paths:
            reader = trk234.Reader(file_path)
            reader.decode(sec_chdo=False, trk_chdo=False)
            sfdu_list = reader.sfdu_list

            # Extract observables data.
            for dtype, converter in self.converters.items():
                extracted = converter.extract(sfdu_list)
                if not extracted.empty:
                    extracted_data[dtype].append(extracted)

            # Extract ramp data.
            ramp_data_list.append(self.ramp_converter.extract(sfdu_list))

        # Process observables data: merge extracted DataFrames and process them.
        tracking_data_list: list[TrackingData] = []
        for dtype, converter in self.converters.items():
            if extracted_data[dtype]:
                merged_df = pd_concat(extracted_data[dtype], ignore_index=True)
                if not merged_df.empty:
                    tracking_data_list.extend(converter.process(merged_df, self.spacecraft_name))

        # Process ramp data into tracking supplementary data.
        supplementary_data_list = self._build_supplementary_data(ramp_data_list, open_ramp_handling)

        return tracking_data_list, supplementary_data_list

    def _build_supplementary_data(
        self,
        ramp_data_list: list,
        open_ramp_handling: OpenRampHandling,
    ) -> list[TrackingSupplementaryData]:
        """
        Merge the per-file ramp data and build one ``TrackingSupplementaryData`` object per station,
        holding a ``RampedFrequencySupplementaryData`` container with that station's merged ramps.
        """
        all_ramps = pd_concat(ramp_data_list, ignore_index=True)
        if all_ramps.empty:
            return []

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

        supplementary_data_list = []
        for station in ramp_df["station"].unique():
            station_df = ramp_df[ramp_df["station"] == station]
            frequency_data = RampedFrequencySupplementaryData()
            for _, row in station_df.iterrows():
                frequency_data.add_frequency_ramp(
                    row["start_time_seconds"],
                    row["end_time_seconds"],
                    row["freq"],
                    row["rate"],
                )

            supplementary_data = TrackingSupplementaryData("Earth", station)
            supplementary_data.set_frequency_supplementary_data([frequency_data])
            supplementary_data_list.append(supplementary_data)

        return supplementary_data_list

    def set_transponder_turnaround_ratio(self, bodies: SystemOfBodies) -> None:
        """
        Set the default transponder turnaround ratio function for the spacecraft, if a spacecraft
        name was provided at initialization.

        NOTE: It's not optimal to set the transponder turnaround ratio here, but it's done for now.

        Parameters
        ----------
        bodies : SystemOfBodies
            The simulation bodies container.
        """
        if self.spacecraft_name:
            spacecraft = bodies.get(self.spacecraft_name)
            spacecraft.system_models.set_default_transponder_turnaround_ratio_function()
