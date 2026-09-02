from ..implementations.oem import OEMMetadata, OEMStateVector
from ..implementations.tdm import TDMMetadata
from ..implementations.omm import (
    OMMMetadata,
    OMMData,
    OMM_HEADER_KEYWORDS,
    OMM_METADATA_KEYWORDS,
)
from ..implementations.common import CCSDSHeader
from .message import (
    NDMMessage,
    TDMMessage,
    OEMMessage,
    OMMMessage,
)
from .message import Segment
from ..ccsds_utils import (
    reconstruct_covariance,
    covariance_matrix_from_named_keywords,
)
from pathlib import Path

# CCSDS 502.0-B-3 Tables 4-1 / 4-2: the OMM header and metadata keywords.
# Unlike OEM/TDM, the OMM has no META_START/META_STOP or DATA_START/
# DATA_STOP delimiters at all (7.4.1, 7.8.8) -- an OMM file is one flat run
# of "KEYWORD = VALUE" lines, and header/metadata/data section boundaries
# are inferred purely from which CCSDS table a keyword belongs to (see
# `NDMParser._process_omm_line`). Everything else -- mean elements,
# spacecraft parameters, TLE-related parameters, the 21 named covariance
# keywords, USER_DEFINED_x -- is "data" by elimination. OMM_HEADER_KEYWORDS/
# OMM_METADATA_KEYWORDS live in `implementations.omm` (derived from
# CCSDSHeader/OMMMetadata's own fields) rather than here, so they can't
# drift out of sync with those models.


class NDMParser:
    """
    A parser for CCSDS Navigation Data Messages (NDM): Orbit Ephemeris
    Messages (OEM), Tracking Data Messages (TDM), and Orbit Mean-Elements
    Messages (OMM).

    OEM/TDM parsing is state-based, keyed on the *_START/*_STOP delimiter
    lines defined in their respective CCSDS blue books. OMM parsing is
    keyword-based instead, since the OMM has no such delimiters -- see
    `_process_omm_line`/`_finalize_omm_segment`.
    """

    def __init__(self):
        self.message_type = None  # This will store the Class (OEMMessage, TDMMessage, etc...)
        self.reset()

    def reset(self):
        """Reset internal state for a fresh parse."""
        self.header_raw = {}
        self.segments = []
        self._state = "HEADER"
        self._current_meta = {}
        self._current_data = {}
        self._current_cov_raw = {}
        self._reconstructed_cov = None
        self._current_cov_epoch = None

    def parse(self, file_path: str | Path) -> NDMMessage | TDMMessage | OEMMessage | OMMMessage:
        """
        Parse a CCSDS message file into its corresponding message object.

        Parameters
        ----------
        file_path : str | Path
            Path of the file to parse.

        Returns
        -------
        NDMMessage | TDMMessage | OEMMessage | OMMMessage
            The populated message object containing header and segments.
        """
        self.reset()

        # 1. IDENTIFY TYPE FIRST: Necessary for correct state transitions
        self.message_type = self._get_message_type(file_path)

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("COMMENT"):
                    # A genuine 'COMMENT = <value>' is only
                    # captured while inside a TDM/OEM META block, where
                    # "COMMENT" unambiguously belongs to that segment's
                    # metadata, so that state tracking makes this safe.
                    # Everything else starting with COMMENT is still dropped:
                    # header-level comments (CCSDSHeader.COMMENT is
                    # list-typed, so a bare string would fail validation)
                    # and freeform 'COMMENT <text>' lines (no '=', the
                    # canonical CCSDS comment syntax). OMM comments are
                    # dropped too, for a different reason: there is
                    # no positional information to tell a header comment
                    # from a metadata one by keyword alone.
                    if self._state == "META" and "=" in line:
                        self._process_line(line)
                    continue
                self._process_line(line)

        # 2. FINALIZE:
        #    - OEM has no DATA_STOP, so close the last segment at EOF.
        #    - OMM has no delimiters at all, so it's always finalized at EOF.
        if self.message_type is OEMMessage and self._current_data:
            self._finalize_current_segment()
        elif self.message_type is OMMMessage and self._current_meta:
            self._finalize_omm_segment()

        # 3. BUILD FINAL MESSAGE
        header = CCSDSHeader(**self._normalize_version_keyword(self.header_raw))
        return self.message_type(header=header, segments=self.segments)

    @staticmethod
    def _normalize_version_keyword(header_raw: dict) -> dict:
        """
        Renames the parsed header's type-specific version keyword
        (CCSDS_OEM_VERS, CCSDS_OMM_VERS, CCSDS_TDM_VERS -- 7.9.1) to the
        generic `CCSDS_VERS` key `CCSDSHeader` actually declares as a
        field.

        Without this, pydantic silently ignores the unrecognized
        `CCSDS_xxx_VERS` key (extra kwargs are ignored, not forbidden, by
        default) and `CCSDS_VERS` always falls back to its "2.0" default
        instead of the version actually read from the file.
        """
        normalized = dict(header_raw)
        for key in list(normalized):
            if key.startswith("CCSDS_") and key.endswith("_VERS") and key != "CCSDS_VERS":
                normalized["CCSDS_VERS"] = normalized.pop(key)
        return normalized

    def _get_message_type(self, file_path: str):
        """
        Identify message class based on file extension.

        Parameters
        ----------
        file_path : str | Path
            The path to the file being parsed.

        Returns
        -------
        Type[NDMMessage]
            The message class (TDMMessage, OEMMessage, or OMMMessage).

        Raises
        ------
        ValueError
            If the file extension does not contain 'TDM', 'OEM', or 'OMM'.
        """
        ext = str(file_path).split(".")[-1].upper()
        if "TDM" in ext:
            return TDMMessage
        elif "OEM" in ext:
            return OEMMessage
        elif "OMM" in ext:
            return OMMMessage
        raise ValueError(f"Unsupported file type: {file_path}")

    def _finalize_current_segment(self):
        """
        Validate metadata and bundle it with data/covariance into a Segment.

        This method instantiates the appropriate metadata model (OEM or TDM),
        creates a Segment object, attaches any reconstructed covariance data,
        and appends the segment to the internal list.
        """
        meta_model = OEMMetadata if self.message_type is OEMMessage else TDMMetadata
        meta_obj = meta_model(**self._current_meta)

        new_segment = Segment(meta_obj, self._current_data)

        # Attach uncertainty data if found
        if self._reconstructed_cov:
            new_segment.covariance = self._reconstructed_cov
            self._reconstructed_cov = None  # Clear for potential next segment

        self.segments.append(new_segment)

    def _finalize_omm_segment(self):
        """
        Builds the OMMMetadata/OMMData pair (plus covariance and any
        USER_DEFINED_x parameters) from the flat keyword dicts accumulated
        by `_process_omm_line`, and appends the resulting (single) Segment.

        `self._current_data` at this point holds every non-header,
        non-metadata keyword seen: the mean-elements/spacecraft/TLE
        fields OMMData models, the 21 named covariance keywords (if any),
        and any USER_DEFINED_x keywords -- all three are split back out
        here since OMMData only models the first group.
        """
        meta_obj = OMMMetadata(**self._current_meta)

        user_defined = {
            key: val for key, val in self._current_data.items() if key.startswith("USER_DEFINED_")
        }
        omm_data_fields = {
            key: val for key, val in self._current_data.items() if key in OMMData.model_fields
        }
        data_obj = OMMData(**omm_data_fields)

        new_segment = Segment(meta_obj, data_obj)

        covariance = covariance_matrix_from_named_keywords(self._current_data)
        if covariance is not None:
            # Dict-shaped (single entry) to mirror OEM's multi-epoch
            # covariance, even though an OMM only ever has the one EPOCH.
            new_segment.covariance = {data_obj.EPOCH: covariance}
        if user_defined:
            new_segment.user_defined = user_defined

        self.segments.append(new_segment)

    def _process_line(self, line: str):
        """
        Process a single line of the CCSDS message, updating state or data.

        This method handles state transitions (e.g., META_START, COVARIANCE_STOP)
        and parses key-value pairs or positional data based on the current state.

        Parameters
        ----------
        line : str
            The cleaned line from the CCSDS file.
        """
        # OMM has no *_START/*_STOP delimiters to drive a state machine
        # with -- classify by keyword instead. See _process_omm_line.
        if self.message_type is OMMMessage:
            self._process_omm_line(line)
            return

        # --- State Transitions ---
        if line == "META_START":
            # For multi-segment files, finalize the previous one before starting new
            if self.message_type is OEMMessage and self._current_data:
                self._finalize_current_segment()
            self._state = "META"
            self._current_meta = {}
            self._current_data = {}
            return

        elif line == "META_STOP":
            # OEM transitions to DATA; TDM returns to HEADER to wait for DATA_START
            self._state = "DATA" if self.message_type is OEMMessage else "HEADER"
            return

        elif line == "DATA_START" and self.message_type is TDMMessage:
            self._state = "DATA"
            return

        elif line == "DATA_STOP" and self.message_type is TDMMessage:
            self._finalize_current_segment()
            self._state = "HEADER"
            return

        elif line == "COVARIANCE_START":
            self._state = "COVARIANCE"
            self._current_cov_raw = {}
            # Reset rather than let a previous covariance block's epoch leak
            # into this one -- see the EPOCH-missing check below.
            self._current_cov_epoch = None
            return

        elif line == "COVARIANCE_STOP":
            # Perform reconstruction and hold in temp variable
            reconstructed = {
                epoch: reconstruct_covariance(lines)
                for epoch, lines in self._current_cov_raw.items()
            }
            self._reconstructed_cov = reconstructed
            self._state = "DATA" if self.message_type is OEMMessage else "HEADER"
            return

        # --- Data Parsing ---
        if "=" in line:
            key, val = [x.strip() for x in line.split("=", 1)]
            if self._state == "HEADER":
                self.header_raw[key] = val
            elif self._state == "META":
                self._current_meta[key] = val
            elif self._state == "COVARIANCE" and "EPOCH" in key:
                self._current_cov_epoch = val
                self._current_cov_raw[self._current_cov_epoch] = []
            elif self._state == "DATA" and self.message_type is TDMMessage:
                # TDM Measurement format: KEYWORD = TIMETAG VALUE
                parts = val.split()
                if len(parts) >= 2:
                    time, value = parts[0], parts[1]
                    if time not in self._current_data:
                        self._current_data[time] = {}
                    self._current_data[time][key] = float(value)

        elif self._state == "DATA" and self.message_type is OEMMessage:
            # OEM Ephemeris format: EPOCH [DATE] X Y Z VX VY VZ
            parts = line.split()
            epoch = parts[0] if len(parts) == 7 else f"{parts[0]} {parts[1]}"
            v = parts[1:] if len(parts) == 7 else parts[2:]
            self._current_data[epoch] = OEMStateVector(
                **{l: float(val) for l, val in zip(["x", "y", "z", "vx", "vy", "vz"], v)}
            )

        elif self._state == "COVARIANCE":
            if not self._current_cov_epoch:
                # Per CCSDS 502.0-B-3 5.2.5.3, a covariance matrix's EPOCH
                # must be provided; without it these rows can't be attached
                # to any epoch, so raise here instead of silently dropping
                # them (the previous behavior).
                raise ValueError(
                    "COVARIANCE block has matrix data before an 'EPOCH' "
                    "keyword -- OEM covariance data lines must be preceded "
                    "by 'EPOCH = <value>' (CCSDS 502.0-B-3 5.2.5.3)."
                )
            # Accumulate raw matrix lines for later reconstruction
            self._current_cov_raw[self._current_cov_epoch].append(line)

    def _process_omm_line(self, line: str):
        """
        Classifies a single OMM line by which CCSDS table its keyword
        belongs to and stashes it accordingly. OMM lines are always
        `KEYWORD = VALUE` (7.4.1) -- there are no bare marker lines like
        `META_START` to handle, since the OMM has none.

        Any keyword that is neither a header nor a metadata keyword is
        data by elimination: mean elements, spacecraft parameters,
        TLE-related parameters, the 21 named covariance keywords, or
        USER_DEFINED_x (Table 4-3) -- `_finalize_omm_segment` sorts those
        apart once the whole file has been read.
        """
        if "=" not in line:
            return
        key, val = [x.strip() for x in line.split("=", 1)]
        if key in OMM_HEADER_KEYWORDS:
            self.header_raw[key] = val
        elif key in OMM_METADATA_KEYWORDS:
            self._current_meta[key] = val
        else:
            self._current_data[key] = val
