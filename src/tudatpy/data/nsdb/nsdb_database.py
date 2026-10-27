from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any
from uuid import UUID

from tudatpy.data.nsdb.nsdb_models import (
    LinkKey,
    LinkRecord,
    NSDBDataset,
    ObservationGroup,
    ObservationGroupKey,
    ObservationQuery,
    Observatory,
    Target,
)


class NSDBDatabase:

    def __init__(self):

        # ----------------------------------------------------
        # Primary registries
        # ----------------------------------------------------

        self.targets: dict[str, Target] = {}

        self.observatories: dict[str, Observatory] = {}

        self.datasets: dict[UUID, NSDBDataset] = {}

        self.links: dict[LinkKey, LinkRecord] = {}

        self.observation_groups: dict[
            ObservationGroupKey,
            ObservationGroup,
        ] = {}

        # ----------------------------------------------------
        # Secondary indexes
        # ----------------------------------------------------

        self.groups_by_target = defaultdict(set)

        self.groups_by_observatory = defaultdict(set)

        self.groups_by_observable = defaultdict(set)

        self.groups_by_dataset = defaultdict(set)

    # ========================================================
    # TARGETS
    # ========================================================

    def register_target(self, target: Target) -> None:
        if target.target_id in self.targets:
            raise ValueError(f"Target already registered: {target.target_id}")

        self.targets[target.target_id] = target

    # ========================================================
    # OBSERVATORIES
    # ========================================================

    def register_observatory(
        self,
        observatory: Observatory,
    ) -> None:

        if observatory.observatory_id in self.observatories:
            raise ValueError("Observatory already registered: " f"{observatory.observatory_id}")

        self.observatories[observatory.observatory_id] = observatory

    # ========================================================
    # DATASETS
    # ========================================================

    def register_dataset(
        self,
        dataset: NSDBDataset,
    ) -> None:

        if dataset.dataset_id in self.datasets:
            raise ValueError(f"Dataset already registered: {dataset.dataset_id}")

        self.datasets[dataset.dataset_id] = dataset

    # ========================================================
    # LINKS
    # ========================================================

    def get_or_create_link(
        self,
        target_id: str,
        observatory_id: str,
    ) -> LinkRecord:

        if target_id not in self.targets:
            raise KeyError(f"Unknown target: {target_id}")

        if observatory_id not in self.observatories:
            raise KeyError(f"Unknown observatory: {observatory_id}")

        key = LinkKey(
            target_id=target_id,
            observatory_id=observatory_id,
        )

        if key not in self.links:
            self.links[key] = LinkRecord(key=key)

        return self.links[key]

    # ========================================================
    # OBSERVATION GROUPS
    # ========================================================

    def register_observation_group(
        self,
        group: ObservationGroup,
    ) -> None:

        key = group.key

        if key in self.observation_groups:
            raise ValueError(f"Observation group already registered: {key}")

        if key.link_key not in self.links:
            raise KeyError(f"Unknown link: {key.link_key}")

        for dataset_id in group.source_dataset_ids:
            if dataset_id not in self.datasets:
                raise KeyError(f"Unknown dataset: {dataset_id}")

        self.observation_groups[key] = group

        # Maintain secondary indexes.
        self.groups_by_target[key.link_key.target_id].add(key)

        self.groups_by_observatory[key.link_key.observatory_id].add(key)

        self.groups_by_observable[key.observable_type].add(key)

        for dataset_id in group.source_dataset_ids:
            self.groups_by_dataset[dataset_id].add(key)

    # ========================================================
    # QUERY
    # ========================================================

    def query(
        self,
        query: ObservationQuery,
    ) -> list[ObservationGroup]:

        candidate_sets = []

        if query.target_id is not None:
            candidate_sets.append(
                self.groups_by_target.get(
                    query.target_id,
                    set(),
                )
            )

        if query.observatory_id is not None:
            candidate_sets.append(
                self.groups_by_observatory.get(
                    query.observatory_id,
                    set(),
                )
            )

        if query.observable_type is not None:
            candidate_sets.append(
                self.groups_by_observable.get(
                    query.observable_type,
                    set(),
                )
            )

        if query.dataset_id is not None:
            candidate_sets.append(
                self.groups_by_dataset.get(
                    query.dataset_id,
                    set(),
                )
            )

        if not candidate_sets:
            keys = set(self.observation_groups)

        else:
            keys = set.intersection(*(set(candidate_set) for candidate_set in candidate_sets))

        return [self.observation_groups[key] for key in keys]
