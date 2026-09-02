from pydantic import BaseModel, Field, field_validator, ConfigDict
import re


class CCSDSHeader(BaseModel):
    """
    Global Header for all Navigation Data Messages (NDM) files.

    Attributes:
        CCSDS_VERS (str): The version of the CCSDS NDM standard. Defaults to "2.0".
            Populated from the type-specific version keyword (e.g.
            CCSDS_OEM_VERS, CCSDS_OMM_VERS, CCSDS_TDM_VERS) -- see
            `NDMParser._normalize_version_keyword`, which renames that key
            to CCSDS_VERS before this model is constructed, since the file
            format's version keyword name varies by message type (7.9.1,
            CCSDS 502.0-B-3) but this shared header model has one field.
        COMMENT (list[str] | None): Optional comments for the header.
        CLASSIFICATION (str | None): User-defined free-text classification/
            caveats for the message (e.g. OEM Table 5-2, OMM Table 4-1).
        CREATION_DATE (str): The date and time the file was created in UTC.
        ORIGINATOR (str): The organization creating the file. Defaults to
            "Tudat"; pass an explicit value to identify the actual
            producing organization/tool.
        MESSAGE_ID (str | None): ID uniquely identifying a message from a
            given originator (e.g. OEM Table 5-2, OMM Table 4-1).
    """

    model_config = ConfigDict(populate_by_name=True)

    CCSDS_VERS: str = Field(default="1.0")
    COMMENT: list[str] | None = []
    CLASSIFICATION: str | None = None
    CREATION_DATE: str
    ORIGINATOR: str = "TU Delft Astrodynamics Toolbox"
    MESSAGE_ID: str | None = None

    @field_validator("CREATION_DATE")
    @classmethod
    def validate_date(cls, v: str):
        """
        Validates that the CREATION_DATE follows the CCSDS ASCII time code formats.

        Supports both Calendar (YYYY-MM-DD) and Day of Year (YYYY-DDD) formats.

        Args:
            v (str): The date string to validate.

        Returns:
            str: The validated date string.
        """
        # Pattern 1: YYYY-MM-DD (Calendar)
        calendar_pattern = r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?Z?$"
        # Pattern 2: YYYY-DDD (Day of Year)
        doy_pattern = r"^\d{4}-\d{3}T\d{2}:\d{2}:\d{2}(\.\d+)?Z?$"

        if re.match(calendar_pattern, v) or re.match(doy_pattern, v):
            return v

        raise ValueError(
            f"Invalid CCSDS date format: {v}. " f"Expected YYYY-MM-DD or YYYY-DDD format."
        )
