import os
import numpy as np
from datetime import datetime, timezone
from typing import Union
from .message import (
    OEMMessage,
    TDMMessage,
    OMMMessage,
    NDMMessage,
)
from ..ccsds_utils import (
    make_lower_triangular,
    covariance_matrix_to_named_keywords,
)


class NDMWriter:
    def __init__(self):
        pass

    def validate_before_write(self, message: NDMMessage):
        """
        Runs Pydantic validation on the header and all segments.
        Raises ValidationError if mandatory fields are missing or logic is inconsistent.

        Paramters
        ---------
        message: NDMMessage
            message object to validate
        """
        # 1. Validate Header (Simple check as it's usually a CCSDSHeader model)
        # Assuming msg.header is already a Pydantic object, otherwise:
        # CCSDSHeader(**msg.header.model_dump())

        # 2. Validate Segments
        for i, segment in enumerate(message.segments):
            # This triggers Pydantic's internal validation (type checks, ranges, etc.)
            segment.metadata.model_validate(segment.metadata)

            # 3. Custom Logic Checks (Astro-logic)
            if hasattr(segment.metadata, "START_TIME") and hasattr(segment.metadata, "STOP_TIME"):
                if segment.metadata.STOP_TIME < segment.metadata.START_TIME:
                    raise ValueError(f"Segment {i}: STOP_TIME is earlier than START_TIME.")

    def write(
        self,
        message: Union[OEMMessage, TDMMessage, OMMMessage],
        output_path: str,
        input_message_in_meters: bool = True,
    ):
        """
        Validates the message and writes it to KVN format.

        Parameters
        ----------
        message : Union[OEMMessage, TDMMessage, OMMMessage]
            The CCSDS message object to be written.
        output_path : str
            The file system path where the output file will be created.
        input_message_in_meters : bool, optional
            If True, assumes input OEM coordinates are in meters and converts
            them to kilometers for the output file. Defaults to True. Not
            applicable to OMM, whose mean elements have no meters/kilometers
            ambiguity to resolve (ignored for OMMMessage).
        """
        # Step 0: Validate
        self.validate_before_write(message)

        # Step 1: Dispatch to specific writer
        if isinstance(message, OEMMessage):
            self._write_oem(message, output_path, input_message_in_meters)
        elif isinstance(message, TDMMessage):
            self._write_tdm(message, output_path)
        elif isinstance(message, OMMMessage):
            self._write_omm(message, output_path)

    def _write_oem(self, msg: OEMMessage, path: str, input_message_in_meters: bool):
        """
        Internal method to write an OEMMessage to a file in KVN format.

        Parameters
        ----------
        msg : OEMMessage
            The OEM message object.
        path : str
            The output file path.
        input_message_in_meters : bool
            Whether the input data is in meters (will be converted to km).
        """

        with open(path, "w") as f:
            # Write Header
            f.write(f"CCSDS_OEM_VERS = {msg.header.CCSDS_VERS}\n")
            if msg.header.CLASSIFICATION:
                f.write(f"CLASSIFICATION = {msg.header.CLASSIFICATION}\n")
            f.write(
                f"CREATION_DATE  = {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')}Z\n"
            )
            f.write(f"ORIGINATOR     = {msg.header.ORIGINATOR}\n")
            if msg.header.MESSAGE_ID:
                f.write(f"MESSAGE_ID = {msg.header.MESSAGE_ID}\n")
            f.write("\n")

            for seg in msg.segments:
                f.write("META_START\n")
                # model_dump(exclude_none=True) keeps the file clean
                for key, val in seg.metadata.model_dump(exclude_none=True).items():
                    f.write(f"{key.upper()} = {val}\n")
                f.write("META_STOP\n\n")

                if input_message_in_meters:
                    for epoch, state in sorted(seg.data.items()):
                        # Divide by 1000 to convert meters to kilometers for the file
                        f.write(
                            f"{epoch} {state['x']/1000.0} {state['y']/1000.0} {state['z']/1000.0} "
                            f"{state['vx']/1000.0} {state['vy']/1000.0} {state['vz']/1000.0}\n"
                        )

                else:
                    # Ephemeris Data
                    for epoch, state in sorted(seg.data.items()):
                        # Standard OEM order: X Y Z VX VY VZ
                        f.write(
                            f"{epoch} {state['x']} {state['y']} {state['z']} "
                            f"{state['vx']} {state['vy']} {state['vz']}\n"
                        )
                # Covariance Data
                if hasattr(seg, "covariance") and seg.covariance:
                    f.write("\nCOVARIANCE_START\n")
                    for epoch, matrix in sorted(seg.covariance.items()):
                        f.write(f"EPOCH = {epoch}\n")

                        # Core Fix: If the input is in meters, handle the 6x6 matrix block scaling
                        if input_message_in_meters:
                            scaled_matrix = np.array(matrix, dtype=float).copy()

                            # Position-Position block (Rows 0-2, Cols 0-2) -> scale from m^2 to km^2 (divide by 1e6)
                            scaled_matrix[0:3, 0:3] /= 1e6

                            # Position-Velocity cross terms (Rows 0-2/Cols 3-5 AND Rows 3-5/Cols 0-2) -> scale by (m * m/s) to (km * km/s) (divide by 1e6)
                            scaled_matrix[0:3, 3:6] /= 1e6
                            scaled_matrix[3:6, 0:3] /= 1e6

                            # Velocity-Velocity block (Rows 3-5, Cols 3-5) -> scale from (m/s)^2 to (km/s)^2 (divide by 1e6)
                            scaled_matrix[3:6, 3:6] /= 1e6

                            cov_lines = make_lower_triangular(scaled_matrix)
                        else:
                            cov_lines = make_lower_triangular(matrix)

                        for line in cov_lines:
                            f.write(f"{line}\n")

                    f.write("COVARIANCE_STOP\n")

    def _write_tdm(self, msg: TDMMessage, path: str):
        """
        Internal method to write a TDMMessage to a file in KVN format.

        Parameters
        ----------
        msg : TDMMessage
            The TDM message object.
        path : str
            The output file path.
        """
        with open(path, "w") as f:
            f.write(f"CCSDS_TDM_VERS = {msg.header.CCSDS_VERS}\n")
            if msg.header.CLASSIFICATION:
                f.write(f"CLASSIFICATION = {msg.header.CLASSIFICATION}\n")
            f.write(
                f"CREATION_DATE  = {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')}Z\n"
            )
            f.write(f"ORIGINATOR     = {msg.header.ORIGINATOR}\n")
            if msg.header.MESSAGE_ID:
                f.write(f"MESSAGE_ID = {msg.header.MESSAGE_ID}\n")
            f.write("\n")

            for seg in msg.segments:
                f.write("META_START\n")
                for key, val in seg.metadata.model_dump(exclude_none=True).items():
                    f.write(f"{key.upper()} = {val}\n")
                f.write("META_STOP\n\n")

                f.write("DATA_START\n")
                for epoch, obs in sorted(seg.data.items()):
                    for keyword, value in obs.items():
                        f.write(f"{keyword.upper()} = {epoch} {value}\n")
                f.write("DATA_STOP\n\n")

    def _write_omm(self, msg: OMMMessage, path: str):
        """
        Internal method to write an OMMMessage to a file in KVN format.

        Unlike OEM/TDM, the entire OMM -- header, metadata, and data -- is
        flat `KEYWORD = VALUE` lines with no `*_START`/`*_STOP` delimiters
        (7.4.1, 7.8.8; see `NDMParser._process_omm_line`).

        Parameters
        ----------
        msg : OMMMessage
            The OMM message object (always exactly one segment).
        path : str
            The output file path.
        """
        segment = msg.segments[0]

        with open(path, "w") as f:
            f.write(f"CCSDS_OMM_VERS = {msg.header.CCSDS_VERS}\n")
            if msg.header.CLASSIFICATION:
                f.write(f"CLASSIFICATION = {msg.header.CLASSIFICATION}\n")
            f.write(
                f"CREATION_DATE  = {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')}Z\n"
            )
            f.write(f"ORIGINATOR     = {msg.header.ORIGINATOR}\n")
            if msg.header.MESSAGE_ID:
                f.write(f"MESSAGE_ID = {msg.header.MESSAGE_ID}\n")
            f.write("\n")

            # Metadata (Table 4-2) -- exclude_none keeps the file clean.
            for key, val in segment.metadata.model_dump(exclude_none=True).items():
                f.write(f"{key.upper()} = {val}\n")
            f.write("\n")

            # Data (Table 4-3): mean elements / spacecraft / TLE-related fields.
            for key, val in segment.data.model_dump(exclude_none=True).items():
                f.write(f"{key.upper()} = {val}\n")

            # Covariance (also Table 4-3): a single {EPOCH: matrix} entry --
            # see NDMParser._finalize_omm_segment -- flattened back into the
            # 21 named keywords (no COVARIANCE_START/STOP; that's OEM/OCM only).
            if segment.covariance:
                for epoch, matrix in segment.covariance.items():
                    for key, val in covariance_matrix_to_named_keywords(matrix).items():
                        f.write(f"{key} = {val:.14e}\n")

            # User-Defined Parameters (Table 4-3): open-ended, so not part
            # of the OMMData model -- see Segment.user_defined.
            if segment.user_defined:
                for key, val in segment.user_defined.items():
                    f.write(f"{key} = {val}\n")
