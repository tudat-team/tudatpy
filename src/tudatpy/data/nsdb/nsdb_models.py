from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from uuid import UUID, uuid4

# ============================================================
# IDENTITY KEYS
# ============================================================


@dataclass(frozen=True)
class LinkKey:
    """Unique identity of a target-observatory link."""

    target_id: str
    observatory_id: str


@dataclass(frozen=True)
class ObservationGroupKey:
    """
    Unique identity of an observation group.

    One observable type + one link configuration.
    """

    observable_type: Any
    link_key: LinkKey


# ============================================================
# DOMAIN ENTITIES
# ============================================================


@dataclass
class Target:
    """A target represented in the NSDB data."""

    target_id: str
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class Observatory:
    """An observatory represented in the NSDB data."""

    observatory_id: str
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class NSDBDataset:
    """An imported NSDB source dataset."""

    dataset_id: UUID
    path: Path
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_path(cls, path: str | Path) -> "NSDBDataset":
        path = Path(path)

        return cls(
            dataset_id=uuid4(),
            path=path,
        )


# ============================================================
# LINK REGISTRY
# ============================================================


@dataclass
class LinkRecord:
    """
    Database record for one unique link.

    The LinkDefinition is constructed lazily.
    """

    key: LinkKey
    link_definition: Any = None


# ============================================================
# OBSERVATION GROUPS
# ============================================================


@dataclass
class ObservationGroup:
    """
    A homogeneous collection of observations.

    All observations belong to the same observable type
    and link configuration.
    """

    key: ObservationGroupKey

    observations: Any

    source_dataset_ids: set[UUID] = field(default_factory=set)

    metadata: dict[str, Any] = field(default_factory=dict)

    single_observation_set: Any = None


# ============================================================
# QUERY
# ============================================================


@dataclass(frozen=True)
class ObservationQuery:
    """Criteria for selecting observation groups."""

    target_id: str | None = None
    observatory_id: str | None = None
    observable_type: Any = None
    dataset_id: UUID | None = None
