"""Readable representations of observation selection conditions, observation datasets, link end ids and link definitions."""

from tudatpy.kernel.estimation.observable_models_setup.links import LinkDefinition, LinkEndId
from tudatpy.kernel.estimation.observations import ObservationDataset, ObservationSelectionCondition

_COMPARISONS = {
    "time_greater_equal": ">=",
    "time_greater_than": ">",
    "time_less_equal": "<=",
    "time_less_than": "<",
}


def _condition_repr(condition):
    kind = condition.condition_type.name
    if kind in ("and_condition", "or_condition"):
        first, second = condition.child_conditions
        operator = "&" if kind == "and_condition" else "|"
        return f"({_condition_repr(first)}) {operator} ({_condition_repr(second)})"
    if kind == "not_condition":
        return f"~({_condition_repr(condition.child_conditions[0])})"
    if kind == "observable_type":
        return f"observable_type == {condition.observable_type_value.name}"
    if kind == "link_definition":
        return f"link_definition == {condition.link_definition_value!r}"
    if kind == "link_end_type":
        return f"has link end {condition.link_end_type_value.name}"
    if kind == "link_end":
        return f"{condition.link_end_type_value.name} == {condition.link_end_id_value!r}"
    if kind == "set_id":
        return f"set_id == {condition.set_id_value}"
    if kind in _COMPARISONS:
        return f"time {_COMPARISONS[kind]} {float(condition.time_value)}"
    if kind in ("residual_absolute_value_greater_than", "observation_absolute_value_greater_than"):
        value = kind.split("_")[0]
        return f"|{value}| > {list(condition.vector_limit)}"
    if kind in ("dependent_variable_greater_than", "dependent_variable_less_than"):
        comparison = ">" if kind == "dependent_variable_greater_than" else "<"
        settings = type(condition.dependent_variable_settings).__name__
        return f"dependent_variable({settings}) {comparison} {list(condition.vector_limit)}"
    return kind


def _dataset_repr(dataset):
    metadata = dataset.get_metadata()
    observations_per_type = {}
    for set_id, set_metadata in metadata.items():
        name = set_metadata["observable_type"].name
        sets, observations = observations_per_type.get(name, (0, 0))
        observations_per_type[name] = (
            sets + 1,
            observations + dataset.number_of_observations_for_set(set_id),
        )
    if not observations_per_type:
        return "ObservationDataset (empty)"
    times = [float(time) for time in dataset.get_data(fields=["times"])["times"]]
    rejected = len(
        dataset.get_data(ObservationSelectionCondition.rejected(), fields=["times"])["times"]
    )
    lines = [
        f"ObservationDataset with {len(times)} observations ({rejected} rejected) in {len(metadata)} sets, "
        + (
            f"times {min(times)} to {max(times)} s since J2000 (TDB):"
            if times
            else "without observations:"
        )
    ]
    for name, (sets, observations) in sorted(observations_per_type.items()):
        lines.append(f"  {name}: {observations} observations in {sets} sets")
    return "\n".join(lines)


ObservationSelectionCondition.__repr__ = _condition_repr
ObservationDataset.__repr__ = _dataset_repr


def _link_end_id_repr(link_end_id):
    if link_end_id.reference_point:
        return f"LinkEndId({link_end_id.body_name!r}, {link_end_id.reference_point!r})"
    return f"LinkEndId({link_end_id.body_name!r})"


def _link_definition_repr(link_definition):
    link_ends = ", ".join(
        f"{link_end_type.name}: {link_end_id!r}"
        for link_end_type, link_end_id in link_definition.link_ends.items()
    )
    return f"LinkDefinition({{{link_ends}}})"


LinkEndId.__repr__ = _link_end_id_repr
LinkDefinition.__repr__ = _link_definition_repr
