from ..implementations.oem import OEMMetadata
from ..implementations.tdm import TDMMetadata, validate_tdm_metadata_for_data
from ..implementations.omm import OMMMetadata, OMMData, validate_omm_data_for_theory
from typing import Any
from ..implementations.common import CCSDSHeader
from collections import defaultdict
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime

TIMESCALE_CONVERTER = time_representation.default_time_scale_converter()


class Segment:
    """A single block consisting of one Metadata section and one Data section."""

    def __init__(self, metadata: OEMMetadata | TDMMetadata | OMMMetadata, data: Any):
        if isinstance(metadata, TDMMetadata):
            # TDM has per-observable metadata requirements (Table 3-5) --
            # OEM's data section has a single fixed schema, so there's
            # nothing to cross-check for it.
            validate_tdm_metadata_for_data(metadata, data)
        elif isinstance(metadata, OMMMetadata) and isinstance(data, OMMData):
            # OMM's TLE-related data fields are conditionally required
            # based on the metadata's MEAN_ELEMENT_THEORY (Table 4-3).
            validate_omm_data_for_theory(metadata, data)
        self.metadata = metadata
        self.data = data
        self.covariance = None
        # OMM-only (Table 4-3's User-Defined Parameters block): open-ended
        # USER_DEFINED_x keywords, which don't fit a fixed-field pydantic
        # model -- see OMMData's docstring. None unless the source OMM
        # actually had a USER_DEFINED_* keyword.
        self.user_defined: dict[str, str] | None = None


class NDMMessage:
    """The high-level object representing the entire parsed file."""

    def __init__(self, header: CCSDSHeader, segments: list[Segment]):
        self.header = header
        self.segments = segments

    def _add_time_columns_to_df(self, table: pd.DataFrame) -> pd.DataFrame:
        """
        Shared internal helper to add absolute numeric epoch tracking columns
        (epoch_seconds_UTC and epoch_seconds_TDB) to any parsed CCSDS dataframe.
        """
        augmented_table = table.copy()

        # Build UTC epoch seconds via Tudat DateTime representation element-wise.
        # + t.nanosecond / 1e9 recovers sub-microsecond precision when `t` is
        # a pandas Timestamp (nanosecond resolution); plain datetime.datetime
        # objects have no such attribute, but every caller of this method
        # currently passes a DatetimeIndex built via pd.to_datetime(...), so
        # `t` is always a Timestamp here.
        augmented_table["epoch_seconds_UTC"] = [
            DateTime(
                t.year,
                t.month,
                t.day,
                t.hour,
                t.minute,
                t.second + t.microsecond / 1e6 + t.nanosecond / 1e9,
            ).epoch()
            for t in table.index
        ]

        # Convert UTC scale to true dynamic TDB seconds scale since J2000
        augmented_table["epoch_seconds_TDB"] = [
            TIMESCALE_CONVERTER.convert_time(
                input_scale=time_representation.utc_scale,
                output_scale=time_representation.tdb_scale,
                input_value=t_utc,
            )
            for t_utc in augmented_table["epoch_seconds_UTC"]
        ]

        return augmented_table


class OEMMessage(NDMMessage):
    """Specific container for Orbit Ephemeris Messages."""

    def has_covariance(self, segment_idx: int = 0) -> bool:
        """
        Parameters
        ----------
        segment_idx: int
            check whether a given segment has covariance

        Returns
        -------
        bool
            True or False
        """
        return hasattr(self.segments[segment_idx], "covariance") and bool(
            self.segments[segment_idx].covariance
        )

    def time_bounds(self, segment_idx: int | None = None) -> list[str]:
        """
        Returns the earliest and latest timestamps in the message or a specific segment.

        Parameters
        ----------
        segment_idx : int, optional
            The index of the segment to check. If None, checks all segments.

        Returns
        -------
        list[str]
            A list containing [min_time, max_time] as ISO strings.
        """
        if segment_idx is None:
            all_times = []
            for seg in self.segments:
                all_times.extend(seg.data.keys())
            segment_times = all_times
        else:
            segment_times = self.segments[segment_idx].data.keys()

        return [min(segment_times), max(segment_times)]

    def to_ccsds_state_history(self, segment_idx: int | None = None) -> dict[datetime, np.ndarray]:
        """
        Extracts ephemeris data in strict CCSDS compliance (Kilometers and UTC strings).

        Parameters
        ----------
        segment_idx : int, optional
            The index of the segment to extract. If None (default), all segments
            are seamlessly concatenated chronologically.

        Returns
        -------
        dict[datetime, np.ndarray]
            A dictionary mapping native Python datetime epochs to a 6-element NumPy array
            [X, Y, Z, VX, VY, VZ] in km and km/s.
        """
        indices = [segment_idx] if segment_idx is not None else range(len(self.segments))

        # Flatten segments into a temporary dictionary to clear up edge overrides
        flattened_raw = {}
        for idx in indices:
            seg = self.segments[idx]
            if seg.data:
                flattened_raw.update(seg.data)

        # Ensure chronological order of the returned mapping dictionary.
        sorted_epoch_strs = sorted(flattened_raw.keys())
        parsed_epochs = pd.to_datetime(sorted_epoch_strs)

        ccsds_history = {}
        for epoch_dt, epoch_str in zip(parsed_epochs, sorted_epoch_strs):
            state_data = flattened_raw[epoch_str]
            ccsds_history[epoch_dt] = np.array(
                [
                    state_data["x"],
                    state_data["y"],
                    state_data["z"],
                    state_data["vx"],
                    state_data["vy"],
                    state_data["vz"],
                ]
            )

        return ccsds_history

    def to_ccsds_covariance_history(
        self, segment_idx: int | None = None
    ) -> dict[datetime, np.ndarray]:
        """
        Returns the covariance history in strict CCSDS compliance (Kilometers squared and UTC strings).

        Parameters
        ----------
        segment_idx : int, optional
            The index of the segment to extract. If None (default), all segments
            containing covariance data are seamlessly concatenated chronologically.

        Returns
        -------
        dict[datetime, np.ndarray]
            A dictionary mapping native Python datetime epochs to 6x6 NumPy covariance arrays
            in km^2, (km/s)^2, etc.
        """
        indices = [segment_idx] if segment_idx is not None else range(len(self.segments))

        # Flatten all relevant segment covariance blocks into a temporary mapping dictionary
        flattened_raw_cov = {}
        for idx in indices:
            if hasattr(self.segments[idx], "covariance") and self.segments[idx].covariance:
                flattened_raw_cov.update(self.segments[idx].covariance)

        if not flattened_raw_cov:
            return {}

        # Ensure chronological matrix mapping order.
        sorted_epoch_strs = sorted(flattened_raw_cov.keys())

        # See to_ccsds_state_history above: kept as pandas Timestamps
        # (nanosecond precision) rather than lossily converted via
        # .to_pydatetime() (microsecond ceiling), and parsed as one batch
        # call instead of once per row for the same reason.
        parsed_epochs = pd.to_datetime(sorted_epoch_strs)

        ccsds_cov_history = {}
        for epoch_dt, epoch_str in zip(parsed_epochs, sorted_epoch_strs):
            ccsds_cov_history[epoch_dt] = np.array(flattened_raw_cov[epoch_str], dtype=float)

        return ccsds_cov_history


class TDMMessage(NDMMessage):
    """
    Specific container for Tracking Data Messages (TDM).

    A TDM consists of a header and one or more segments. Each segment contains
    a metadata section (TDMMetadata) and a data section containing tracking
    observations (e.g., Range, Doppler, Angles) indexed by time.
    """

    def number_of_segments(self) -> int:
        """
        Returns the total count of metadata/data pairs in the message.

        Returns
        -------
        int
            The number of segments.
        """
        return len(self.segments)

    def time_bounds(self, segment_idx: int | None = None, actual_data: bool = True) -> list[str]:
        """
        Returns the earliest and latest timestamps in the message or a specific segment.

        Parameters
        ----------
        segment_idx : int, optional
            The index of the segment to check. If None, checks all segments.
        actual_data : bool, optional
            If True (default), it looks at the actual tracking observation epochs.
            If False, it looks at the Metadata START_TIME and STOP_TIME fields.

        Returns
        -------
        list[str]
            A list containing [min_time, max_time] as ISO strings.
            Returns ["N/A", "N/A"] if no timestamps are found.
        """
        if not self.segments:
            return ["N/A", "N/A"]

        # 1. Gather relevant timestamps
        all_epochs = []
        indices = [segment_idx] if segment_idx is not None else range(len(self.segments))

        for i in indices:
            seg = self.segments[i]
            if actual_data:
                # Use actual observation epochs
                if seg.data:
                    all_epochs.extend(seg.data.keys())
            else:
                # Use official Metadata bounds
                if hasattr(seg.metadata, "START_TIME") and seg.metadata.START_TIME:
                    all_epochs.append(seg.metadata.START_TIME)
                if hasattr(seg.metadata, "STOP_TIME") and seg.metadata.STOP_TIME:
                    all_epochs.append(seg.metadata.STOP_TIME)

        # 2. Safety check for empty data
        if not all_epochs:
            return ["N/A", "N/A"]

        # 3. Return min and max (More efficient than full sort)
        return [min(all_epochs), max(all_epochs)]


#    @classmethod
#    def from_dataframe(
#        cls, df: pd.DataFrame, norad_id: str, obs_code: str, angle_type: str = "RADEC"
#    ) -> "TDMMessage":
#        """
#        Factory method to create a TDMMessage directly from a pandas DataFrame
#        (dataframe = output of calculate_pass_data)
#
#        This method standardizes time columns, renames physical measurement columns to
#        CCSDS TDM shorthand, and packages the data into a validated TDMMessage.
#
#        Time is assumed to be in UTC already.
#
#        Parameters
#        ----------
#        df : pd.DataFrame
#            Input data containing tracking observations. Optionally include
#            'Sigma_Right_Ascension'/'Sigma_Declination' (RADEC) or
#            'Sigma_Azimuth'/'Sigma_Elevation' (AZEL) columns to also write
#            per-observation 1-sigma angle uncertainty, in the same units as
#            the corresponding angle column, as the non-standard
#            `SIGMA_ANGLE_1`/`SIGMA_ANGLE_2` TDM data keywords.
#            Omitted entirely if these columns aren't present.
#        norad_id : str
#            The ID of the target spacecraft (or observer/station) (i.e. PARTICIPANT_1).
#        obs_code : str
#            The ID of the observer/station (or spacecraft) (i.e. PARTICIPANT_2).
#        angle_type : str, optional
#            The coordinate system for angles ('RADEC' or 'AZEL'). Defaults to "RADEC".
#
#        Returns
#        -------
#        TDMMessage
#            A populated TDMMessage instance.
#        """
#        # Create a copy so we don't mutate the user's original DataFrame
#        working_df = df.copy()
#
#        # 1. Standardize Time to be the DatetimeIndex
#        if "Time" in working_df.columns:
#            working_df = working_df.set_index(pd.to_datetime(working_df["Time"]))
#        elif "UTC_TIME" in working_df.columns:
#            working_df = working_df.set_index(pd.to_datetime(working_df["UTC_TIME"]))
#
#        # 2. Standardize Physics Columns to TDM shorthand
#        working_df = working_df.rename(
#            columns={
#                "Right_Ascension": "RA",
#                "Declination": "DEC",
#                "Azimuth": "AZ",
#                "Elevation": "EL",
#                "Sigma_Right_Ascension": "SIGMA_RA",
#                "Sigma_Declination": "SIGMA_DEC",
#                "Sigma_Azimuth": "SIGMA_AZ",
#                "Sigma_Elevation": "SIGMA_EL",
#            }
#        )
#
#        grouped_data = defaultdict(dict)
#
#        # 3. Group by Index (UTC Time)
#        for timestamp, row in working_df.iterrows():
#            # Format the pandas Timestamp index back to a CCSDS-compliant ISO string
#            if isinstance(timestamp, pd.Timestamp):
#                time_str = timestamp.isoformat()
#                # Ensure the UTC 'Z' designation is present if the timestamp is naive
#                if "+" not in time_str and not time_str.endswith("Z"):
#                    time_str += "Z"
#            else:
#                time_str = str(timestamp)
#
#            if angle_type == "RADEC":
#                grouped_data[time_str]["ANGLE_1"] = row["RA"]
#                grouped_data[time_str]["ANGLE_2"] = row["DEC"]
#                if "SIGMA_RA" in working_df.columns:
#                    grouped_data[time_str]["SIGMA_ANGLE_1"] = row["SIGMA_RA"]
#                if "SIGMA_DEC" in working_df.columns:
#                    grouped_data[time_str]["SIGMA_ANGLE_2"] = row["SIGMA_DEC"]
#            elif angle_type == "AZEL":
#                grouped_data[time_str]["ANGLE_1"] = row["AZ"]
#                grouped_data[time_str]["ANGLE_2"] = row["EL"]
#                if "SIGMA_AZ" in working_df.columns:
#                    grouped_data[time_str]["SIGMA_ANGLE_1"] = row["SIGMA_AZ"]
#                if "SIGMA_EL" in working_df.columns:
#                    grouped_data[time_str]["SIGMA_ANGLE_2"] = row["SIGMA_EL"]
#            else:
#                raise ValueError(f"Unsupported angle type: {angle_type}.")
#
#        if angle_type == "RADEC":
#            reference_frame = "ICRF"
#        else:
#            reference_frame = None
#
#        for time_tag, measurements in grouped_data.items():
#            for measurement_name, value in measurements.items():
#                if measurement_name == "ANGLE_1":
#                    # Extract the raw float, wrap it according to astronomical definitions of RA and Dec
#                    normalized_az_or_ra_angle = value % 360.0
#                    grouped_data[time_tag][measurement_name] = normalized_az_or_ra_angle
#
#        # 2. Build the Pydantic Metadata object
#        metadata = TDMMetadata(
#            TIME_SYSTEM="UTC",
#            START_TIME=min(grouped_data.keys()),
#            STOP_TIME=max(grouped_data.keys()),
#            PARTICIPANT_1=str(norad_id),
#            PARTICIPANT_2=str(obs_code),
#            PATH="1,2",  # Choose this path as default
#            MODE="SEQUENTIAL",
#            ANGLE_TYPE=angle_type,
#            REFERENCE_FRAME=reference_frame,
#        )
#
#        # 3. Create the Segment with the now-normalized data
#        segment = Segment(metadata=metadata, data=dict(grouped_data))
#        # 6. Create a default Header
#        header = CCSDSHeader(
#            CCSDS_VERS="2.0",
#            CREATION_DATE=datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
#        )
#
#        return cls(header=header, segments=[segment])

#    def to_dataframe(
#        self,
#        segment_idx: int | None = None,
#    ) -> pd.DataFrame:
#        """
#        Convert TDM segments (assumed to be in CCSDS-compliant units, and in UTC scale)
#        to a pandas DataFrame.
#
#        This method handles both Day-of-Year (DOY) and ISO timestamps, ensuring the
#        resulting DataFrame is indexed by UTC datetime objects.
#
#        Parameters
#        ----------
#        segment_idx : int, optional
#            The index of the segment to convert. If None (default), all segments
#            are concatenated.
#
#        Returns
#        -------
#        pd.DataFrame
#            A DataFrame indexed by 'UTC_TIME' containing the tracking observations.
#        """
#
#        angle_column_map = {
#            "AZEL": {"ANGLE_1": "AZ", "ANGLE_2": "EL"},
#            "RADEC": {"ANGLE_1": "RA", "ANGLE_2": "DEC"},
#            "X_Y": {"ANGLE_1": "X_ANGLE", "ANGLE_2": "Y_ANGLE"},
#        }
#
#        # -------------------------
#        # Unified time parsing
#        # -------------------------
#        def _parse_datetime(series: pd.Series) -> pd.Series:
#            """Parse DOY and ISO timestamps into pandas datetime."""
#
#            def _parse_single(t):
#                if not isinstance(t, str):
#                    return pd.NaT
#                try:
#                    # Try DOY
#                    return pd.to_datetime(t, format="%Y-%jT%H:%M:%S.%f")
#                except Exception:
#                    # ISO Timestamps
#                    return pd.to_datetime(t, errors="coerce")
#
#            return series.map(_parse_single)
#
#        # -------------------------
#        # Segment selection
#        # -------------------------
#        indices = [segment_idx] if segment_idx is not None else range(len(self.segments))
#        all_dfs = []
#
#        for idx in indices:
#            seg = self.segments[idx]
#            if not seg.data:
#                continue
#
#            angle_type = getattr(seg.metadata, "ANGLE_TYPE", None)
#            col_rename = angle_column_map.get(angle_type, {})
#
#            # --- Correct Participant Mapping Based on PATH ---
#            part1 = getattr(seg.metadata, "PARTICIPANT_1", "UNKNOWN")
#            part2 = getattr(seg.metadata, "PARTICIPANT_2", "UNKNOWN")
#
#            # Strip spaces safely and get the standard PATH format
#            path = getattr(seg.metadata, "PATH", "1,2").replace(" ", "")
#
#            # If path starts with "1" (e.g. "1,2"), Participant 1 is RSO and 2 is Station
#            if path.startswith("1"):
#                rso_name = part1
#                station_name = part2
#            else:
#                # E.g., "2,1", so Participant 2 is RSO and 1 is Station
#                rso_name = part2
#                station_name = part1
#
#            rows = []
#            for time_tag, observations in seg.data.items():
#                row = {"UTC_TIME": time_tag}
#                for keyword, value in observations.items():
#                    col_name = col_rename.get(keyword, keyword)
#                    row[col_name] = value
#
#                # Consistently assign the resolved participant names
#                row["PARTICIPANT_1"] = rso_name
#                row["PARTICIPANT_2"] = station_name
#
#                rows.append(row)
#
#            df = pd.DataFrame(rows)
#
#            if df.empty:
#                continue
#
#            # -------------------------
#            # Time handling
#            # -------------------------
#            df["UTC_TIME"] = _parse_datetime(df["UTC_TIME"])
#            df = df.set_index("UTC_TIME")
#
#            # Convert numeric columns
#            for col in df.columns:
#                if col not in ["PARTICIPANT_1", "PARTICIPANT_2"]:
#                    df[col] = pd.to_numeric(df[col], errors="coerce")
#
#            if segment_idx is None:
#                df["segment_index"] = idx
#
#            all_dfs.append(df)
#
#        if not all_dfs:
#            return pd.DataFrame()
#
#        combined = pd.concat(all_dfs)
#        return self._add_time_columns_to_df(combined).sort_index()


class OMMMessage(NDMMessage):
    """
    Specific container for Orbit Mean-Elements Messages (OMM).

    Unlike OEM/TDM, an OMM is always exactly one segment: CCSDS 502.0-B-3
    4.1.5 states "The OMM shall be a plain text file consisting of orbit
    data for a single object" and, unlike OEM (5.2.1.2), Table 4-3 defines
    no mechanism for repeated metadata/data blocks within one file. The
    `metadata`/`data`/`covariance` properties below are just a convenience
    over `segments[0]` for that reason -- there's never a second segment
    to disambiguate.

    An OMM's data section is itself flat KVN (7.4.1) rather than the
    tabular ephemeris rows OEM uses, so -- unlike `OEMMessage`, which
    stores each segment's `data` as a raw `dict[epoch, state]` -- an
    `OMMMessage` segment's `data` is a typed `OMMData` pydantic model (see
    its docstring for why covariance and user-defined parameters are
    handled separately instead).
    """

    @property
    def metadata(self) -> OMMMetadata:
        """The (only) segment's metadata."""
        return self.segments[0].metadata

    @property
    def data(self) -> OMMData:
        """The (only) segment's mean-elements/spacecraft/TLE data."""
        return self.segments[0].data

    @property
    def covariance(self) -> dict[str, np.ndarray] | None:
        """
        The (only) segment's covariance, if present: a single-entry
        `{EPOCH: 6x6 matrix}` dict (kept dict-shaped, like OEM's
        multi-epoch covariance, even though an OMM only ever has one
        epoch -- see `NDMParser._finalize_omm_segment`).
        """
        return self.segments[0].covariance

    @property
    def user_defined(self) -> dict[str, str] | None:
        """The (only) segment's USER_DEFINED_x parameters, if any."""
        return self.segments[0].user_defined
