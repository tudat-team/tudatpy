from tudatpy.kernel.estimation.observations import ObservationSelectionCondition
from tudatpy.estimation.observable_models_setup.links import (
    observed_body,
    observer,
    receiver,
    transmitter,
)

_INCOMPLETE_SELECTOR_MESSAGE = "Incomplete observation query selector. Compare it to a value or call a supported selector method."


class _ComparableSelector:
    def __init__(self, factory):
        self._factory = factory

    def __eq__(self, value):
        return self._factory(value)

    def __ne__(self, value):
        return ~self._factory(value)

    def __bool__(self):
        raise TypeError(_INCOMPLETE_SELECTOR_MESSAGE)


class _LinkEndSelector:
    def __init__(self, link_end_type):
        self._link_end_type = link_end_type

    def __eq__(self, link_end_id):
        return ObservationSelectionCondition.link_end(self._link_end_type, link_end_id)

    def __ne__(self, link_end_id):
        return ~ObservationSelectionCondition.link_end(self._link_end_type, link_end_id)

    def __bool__(self):
        raise TypeError(_INCOMPLETE_SELECTOR_MESSAGE)


class _TimeSelector:
    def between(self, start_time, end_time):
        return ObservationSelectionCondition.time_bounds(start_time, end_time)

    def __ge__(self, time):
        return ObservationSelectionCondition.time_greater_equal(time)

    def __gt__(self, time):
        return ObservationSelectionCondition.time_greater_than(time)

    def __le__(self, time):
        return ObservationSelectionCondition.time_less_equal(time)

    def __lt__(self, time):
        return ObservationSelectionCondition.time_less_than(time)

    def __bool__(self):
        raise TypeError(_INCOMPLETE_SELECTOR_MESSAGE)


class _VectorValueSelector:
    """Selector for row-level absolute-value comparisons.

    A vector observable row is selected when any scalar component exceeds the
    supplied absolute limit.
    """

    def __init__(self, abs_greater_than_factory):
        self._abs_greater_than_factory = abs_greater_than_factory

    def abs_greater_than(self, limit):
        """Select rows where any component has absolute value above ``limit``."""
        return self._abs_greater_than_factory(limit)

    def __bool__(self):
        raise TypeError(_INCOMPLETE_SELECTOR_MESSAGE)


class _DependentVariableSelector:
    """Selector for row-level dependent-variable comparisons.

    Unlike residual and observation selectors, ``greater_than`` compares signed
    dependent-variable values.
    """

    def __init__(self, settings, return_first_compatible_settings=False):
        self._settings = settings
        self._return_first_compatible_settings = return_first_compatible_settings

    def greater_than(self, limit):
        """Select rows where any compatible dependent-variable component is above ``limit``."""
        return ObservationSelectionCondition.dependent_variable_greater_than(
            self._settings,
            limit,
            self._return_first_compatible_settings,
        )

    def __bool__(self):
        raise TypeError(_INCOMPLETE_SELECTOR_MESSAGE)


class _ObservationQuery:
    """Query builder for creating ObservationSelectionCondition objects."""

    @property
    def observable_type(self):
        return _ComparableSelector(ObservationSelectionCondition.observable_type)

    @property
    def set_id(self):
        return _ComparableSelector(ObservationSelectionCondition.set_id)

    @property
    def link_definition(self):
        return _ComparableSelector(ObservationSelectionCondition.link_definition)

    def link_end(self, link_end_type):
        return _LinkEndSelector(link_end_type)

    @property
    def receiver(self):
        return self.link_end(receiver)

    @property
    def transmitter(self):
        return self.link_end(transmitter)

    @property
    def observer(self):
        return self.link_end(observer)

    @property
    def observed_body(self):
        return self.link_end(observed_body)

    @property
    def time(self):
        return _TimeSelector()

    @property
    def active(self):
        return ObservationSelectionCondition.active()

    @property
    def rejected(self):
        return ObservationSelectionCondition.rejected()

    @property
    def residual(self):
        return _VectorValueSelector(
            ObservationSelectionCondition.residual_absolute_value_greater_than
        )

    @property
    def observation(self):
        return _VectorValueSelector(
            ObservationSelectionCondition.observation_absolute_value_greater_than
        )

    def dependent_variable(self, settings, return_first_compatible_settings=False):
        return _DependentVariableSelector(settings, return_first_compatible_settings)

    def __bool__(self):
        raise TypeError(_INCOMPLETE_SELECTOR_MESSAGE)


observation_query = _ObservationQuery()
