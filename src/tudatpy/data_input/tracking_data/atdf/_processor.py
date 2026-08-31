import os
from pathlib import Path

import pandas as pd
from atdf2ascii import main as atdf2ascii_main

from tudatpy.data_input.tracking_data import (
    TrackingData,
    TrackingSupplementaryData,
)

from ._converters import (
    AtdfNwayDopplerConverter,
    AtdfNwayRangeConverter,
    AtdfRampConverter,
)

_DEFAULT_PROC_COUNT = max(1, (os.cpu_count() or 2) // 2)


class AtdfTrackingDataProcessor:
    """
    Processor for ATDF (TRK-2-25) files using ``atdf2ascii``.

    For the observable groups enabled through the ``doppler_*``/``range_*``
    flags, the processor first runs ``atdf2ascii`` to decode the ATDF files
    into intermediate ``.msr``/``.ramp`` ASCII tables, reads those tables back
    into ``pandas`` DataFrames, and then dispatches them to the range,
    Doppler, and ramp converters. The converters produce
    :class:`~tudatpy.data_input.tracking_data.TrackingData` objects for the
    enabled observable groups, and one
    :class:`~tudatpy.data_input.tracking_data.TrackingSupplementaryData` object
    per ground station, holding that station's frequency ramps.

    Parameters
    ----------
    atdf_file_path : list[pathlib.Path]
        List of ATDF file paths to be processed.
    spacecraft_name : str
        Spacecraft body name used in generated link definitions.
    doppler_one_way, doppler_two_way, doppler_three_way, range_one_way, range_two_way : bool
        Observable groups to be decoded by ``atdf2ascii``. ``doppler_one_way``
        and ``range_one_way`` are reserved for future support and currently
        raise ``NotImplementedError`` if set to ``True``, since no converter
        exists yet for 1-way Doppler/range data.
    """

    atdf_time_tag_format = "%d-%b-%Y %H:%M:%S.%f"

    def __init__(
        self,
        atdf_file_path: list[Path],
        spacecraft_name: str,
        proc_count: int = _DEFAULT_PROC_COUNT,
        doppler_one_way: bool = False,
        doppler_two_way: bool = True,
        doppler_three_way: bool = True,
        range_one_way: bool = False,
        range_two_way: bool = True,
    ):
        self.atdf_file_path = atdf_file_path

        self.output_dir = Path("output")
        self.spacecraft_name = spacecraft_name
        self.proc_count = proc_count

        if doppler_one_way:
            raise NotImplementedError(
                "doppler_one_way=True decodes 1-way Doppler data, but no converter "
                "exists yet for 'Data Type' == '1-Way-Doppler'."
            )
        if range_one_way:
            raise NotImplementedError(
                "range_one_way=True decodes 1-way range data, but no converter "
                "exists yet for 'Data Type' == '1-Way-Range'."
            )

        self.doppler_one_way = doppler_one_way
        self.doppler_two_way = doppler_two_way
        self.doppler_three_way = doppler_three_way
        self.range_one_way = range_one_way
        self.range_two_way = range_two_way

    def convert_atdf_to_ascii(
        self,
        output_dir: Path | None,
        count_time: list[float] | None = None,
    ):
        if output_dir is None:
            output_dir = self.output_dir

        output_dir.mkdir(parents=True, exist_ok=True)

        if count_time is not None:
            count_time = list(count_time)

        for atdf_file in self.atdf_file_path:
            atdf2ascii_main(
                input_file=atdf_file.absolute().as_posix(),
                output_dir=output_dir.as_posix(),
                proc_count=self.proc_count,
                count_time=count_time,
                doppler_one_way=self.doppler_one_way,
                doppler_two_way=self.doppler_two_way,
                doppler_three_way=self.doppler_three_way,
                range_one_way=self.range_one_way,
                range_two_way=self.range_two_way,
            )

    def _process_msr_dataframe(self, df_raw: pd.DataFrame):

        df_processed = df_raw.copy()

        df_processed["time_tag (DT)"] = pd.to_datetime(
            df_processed["time_tag (UTC)"].str.strip(),
            format=self.atdf_time_tag_format,
        )

        return df_processed

    def _process_ramp_dataframe(self, df_raw: pd.DataFrame):

        df_processed = df_raw.copy()

        time_cols = ["Start-Time", "End-Time"]

        for col in time_cols:
            df_processed[f"{col} (DT)"] = df_processed[col].apply(
                lambda col: pd.to_datetime(col.strip(), format=self.atdf_time_tag_format)
            )

        return df_processed

    @staticmethod
    def _read_atdf_ascii_table(atdf_files: list[Path], output_dir: Path, suffix: str):

        df_list = []

        for atdf_file in atdf_files:
            table_file = output_dir / f"{atdf_file.stem}.{suffix}"
            with open(table_file) as f:
                header_line = ""
                for line in f:
                    if line.startswith("#"):
                        header_line = line
                    else:
                        break
            cols = [c.strip() for c in header_line.lstrip("#").split(",")]

            if "Observed (Hz)" in cols:
                cols[cols.index("Observed (Hz)")] = "Observed"
            elif "Observed (RU)" in cols:
                cols[cols.index("Observed (RU)")] = "Observed"

            df_raw = pd.read_csv(
                table_file,
                comment="#",  # skips all '#' lines (preamble + column header)
                header=None,
                names=cols,
                skipinitialspace=True,
            )
            df_list.append(df_raw)

        return pd.concat(df_list, ignore_index=True)

    def read_atdf_ascii_msr(self, output_dir: Path):

        self.df_raw_msr = self._read_atdf_ascii_table(self.atdf_file_path, output_dir, "msr")
        self.df_processed_msr = self._process_msr_dataframe(self.df_raw_msr)

    def read_atdf_ascii_ramp(self, output_dir: Path):

        self.df_raw_rmp = self._read_atdf_ascii_table(self.atdf_file_path, output_dir, "ramp")
        self.df_processed_rmp = self._process_ramp_dataframe(self.df_raw_rmp)

    def nway_doppler_enabled(self) -> bool:
        return self.doppler_two_way or self.doppler_three_way

    def nway_range_enabled(self) -> bool:
        return self.range_two_way

    def any_observable_enabled(self) -> bool:
        return self.nway_doppler_enabled() or self.nway_range_enabled()

    def process(
        self,
        output_dir: Path,
        earth_name: str = "Earth",
    ) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:

        if self.any_observable_enabled():
            self.read_atdf_ascii_msr(output_dir)

        if self.nway_range_enabled():
            rng_converter = AtdfNwayRangeConverter()
            n_way_range_obs = rng_converter.extract(self.df_processed_msr, self.spacecraft_name)
            rng_tracking_data = rng_converter.process(n_way_range_obs, earth_name)
        else:
            rng_tracking_data = []

        if self.nway_doppler_enabled():
            doppler_converter = AtdfNwayDopplerConverter()
            nway_doppler_obs = doppler_converter.extract(
                self.df_processed_msr, self.spacecraft_name
            )
            doppler_tracking_data = doppler_converter.process(nway_doppler_obs, earth_name)
        else:
            doppler_tracking_data = []

        self.read_atdf_ascii_ramp(output_dir)

        ramp_converter = AtdfRampConverter()
        ramp_df = ramp_converter.extract(self.df_processed_rmp)
        supplementary_data_list = ramp_converter.process(ramp_df, earth_name)

        return doppler_tracking_data + rng_tracking_data, supplementary_data_list
