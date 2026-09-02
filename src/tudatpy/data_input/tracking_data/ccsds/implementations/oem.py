from pydantic import BaseModel, model_validator


class OEMStateVector(BaseModel):
    """
    OEM Data section: a single Cartesian state vector at one epoch
    (CCSDS 502.0-B-3 Table 5-3).

    Mirrors `OMMData` in spirit (a validated model for the message type's
    data section) but is keyed by epoch rather than singular, since an OEM's
    data section is a time-indexed sequence of these -- see
    `Segment.data: dict[str, OEMStateVector]` for OEM segments.

    Supports dict-style item access (`state["x"]`) alongside normal
    attribute access (`state.x`), so existing code written against the
    plain-dict representation this replaces keeps working unchanged.
    """

    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float

    def __getitem__(self, key: str) -> float:
        return getattr(self, key)


class OEMMetadata(BaseModel):
    """
    Metadata fields specific to Orbit Ephemeris Messages (OEM).

    Attributes:
        OBJECT_NAME: Spacecraft name.
        OBJECT_ID: Spacecraft identifier.
        CENTER_NAME: Origin of the coordinate system (e.g., EARTH).
        REF_FRAME: Reference frame for the ephemeris data.
        TIME_SYSTEM: Time system used (e.g., UTC).
        START_TIME: Start time of the ephemeris data.
        STOP_TIME: Stop time of the ephemeris data.
        INTERPOLATION: Interpolation method used.
        INTERPOLATION_DEGREE: Degree of the interpolation.
        COMMENT: Optional comments.
    """

    COMMENT: str | None = None
    OBJECT_NAME: str
    OBJECT_ID: str
    CENTER_NAME: str
    REF_FRAME: str
    REF_FRAME_EPOCH: str | None = None
    TIME_SYSTEM: str
    START_TIME: str
    STOP_TIME: str
    USABLE_START_TIME: str | None = None
    USABLE_STOP_TIME: str | None = None
    INTERPOLATION: str | None = None
    INTERPOLATION_DEGREE: int | None = None

    @model_validator(mode="after")
    def validate_interpolation(self):
        if self.INTERPOLATION is not None and self.INTERPOLATION_DEGREE is None:
            raise ValueError(
                "INTERPOLATION_DEGREE must be provided when INTERPOLATION is specified"
            )
        return self
