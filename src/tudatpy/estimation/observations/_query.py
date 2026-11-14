import datetime
import operator

from tudatpy.astro.time_representation import DateTime
from tudatpy.estimation.observable_models_setup import links
from tudatpy.kernel.estimation.observations import ObservationSelectionCondition

# Readable representations of the dataset, conditions and link definitions
from . import _repr  # noqa: F401


def _as_time(time):
    """Convert a DateTime or Python datetime (interpreted as TDB) to a Time (in seconds since J2000)."""
    if isinstance(time, datetime.datetime):
        if time.tzinfo is not None:
            raise ValueError(
                "Times in observation queries must be given without time zone (a datetime is interpreted as TDB)."
            )
        time = DateTime.from_python_datetime(time)
    if isinstance(time, DateTime):
        return time.to_epoch_time_object()
    return time


def _any_of(conditions):
    """Combine conditions with OR as a balanced tree, so that the size and evaluation depth grow slowly with their number."""
    if len(conditions) == 1:
        return conditions[0]
    middle = len(conditions) // 2
    return _any_of(conditions[:middle]) | _any_of(conditions[middle:])


def _as_link_end_id(link_end_id):
    """Convert a (body, station) tuple to a LinkEndId; a string is rejected since it may be a body or a station."""
    if isinstance(link_end_id, tuple):
        return links.LinkEndId(link_end_id)
    if isinstance(link_end_id, str):
        raise TypeError(
            f"A link end must be given as a LinkEndId or a (body, station) tuple, not as the string "
            f"'{link_end_id}': use ('{link_end_id}', '') for a body, or (body, '{link_end_id}') for a station."
        )
    return link_end_id


class _Selector:
    """Base class of the selectors of :data:`observation_query`.

    A selector must be compared to a value (or one of its functions must be called) to create a
    condition. Using a selector itself as a boolean raises a ``TypeError``.
    """

    # Let numpy scalars on the left of a comparison (e.g. ``times[0] <= observation_query.time``) defer to the selector
    __array_ufunc__ = None

    def __bool__(self):
        raise TypeError(
            "Incomplete observation query selector. "
            "Compare it to a value or call a supported selector method."
        )


class _ComparableSelector(_Selector):
    def __init__(self, factory):
        self._factory = factory

    def __eq__(self, value):
        return self._factory(value)

    def __ne__(self, value):
        return ~self._factory(value)

    def isin(self, values):
        """Function to create a condition that selects the observations for which this attribute is equal to any of
        the given values.

        Parameters
        ----------
        values : list
            Values to compare to (an empty list selects no observations).

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations for which the attribute is equal to any of the values.
        """
        conditions = [self._factory(value) for value in values]
        if not conditions:
            return ~ObservationSelectionCondition.all()
        return _any_of(conditions)


class _TimeSelector(_ComparableSelector):
    def __init__(self):
        super().__init__(lambda time: self.between(time, time))

    def between(self, start_time, end_time):
        """Function to create a condition that selects the observations with a time between (and
        including) two bounds.

        Parameters
        ----------
        start_time : tudatpy.astro.time_representation.Time or tudatpy.astro.time_representation.DateTime or datetime.datetime
            Lower bound of the observation time, in seconds since J2000 (TDB); a ``DateTime`` or ``datetime`` is
            interpreted as a TDB epoch.
        end_time : tudatpy.astro.time_representation.Time or tudatpy.astro.time_representation.DateTime or datetime.datetime
            Upper bound of the observation time, in seconds since J2000 (TDB); a ``DateTime`` or ``datetime`` is
            interpreted as a TDB epoch.

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations within the time bounds.
        """
        return ObservationSelectionCondition.time_bounds(_as_time(start_time), _as_time(end_time))

    def __ge__(self, time):
        return ObservationSelectionCondition.time_greater_equal(_as_time(time))

    def __gt__(self, time):
        return ObservationSelectionCondition.time_greater_than(_as_time(time))

    def __le__(self, time):
        return ObservationSelectionCondition.time_less_equal(_as_time(time))

    def __lt__(self, time):
        return ObservationSelectionCondition.time_less_than(_as_time(time))


class _VectorValueSelector(_Selector):
    def __init__(self, abs_greater_than_factory):
        self._abs_greater_than_factory = abs_greater_than_factory

    def abs_greater_than(self, limit):
        """Function to create a condition that selects the observations for which the absolute
        value of any component is larger than a limit.

        Parameters
        ----------
        limit : float or numpy.ndarray
            Limit on the absolute value, either as a single value for all components, or with one
            value per observation component.

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations with a value larger than the limit.
        """
        return self._abs_greater_than_factory(limit)


class _DependentVariableSelector(_Selector):
    def __init__(self, settings, return_first_compatible_settings):
        self._settings = settings
        self._return_first_compatible_settings = return_first_compatible_settings

    def greater_than(self, limit):
        """Function to create a condition that selects the observations for which any component of
        the dependent variable is larger than a limit.

        Parameters
        ----------
        limit : float or numpy.ndarray
            Limit on the (signed) value of the dependent variable, either as a single value for all
            components, or with one value per dependent variable component.

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations with a dependent variable larger than the limit.
        """
        return ObservationSelectionCondition.dependent_variable_greater_than(
            self._settings, limit, self._return_first_compatible_settings
        )

    def less_than(self, limit):
        """Function to create a condition that selects the observations for which any component of
        the dependent variable is smaller than a limit.

        Unlike ``~greater_than(limit)``, this does not select observations in sets for which no
        compatible dependent variable is stored.

        Parameters
        ----------
        limit : float or numpy.ndarray
            Limit on the (signed) value of the dependent variable, either as a single value for all
            components, or with one value per dependent variable component.

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations with a dependent variable smaller than the limit.
        """
        return ObservationSelectionCondition.dependent_variable_less_than(
            self._settings, limit, self._return_first_compatible_settings
        )


class _ObservationIdSelector(_ComparableSelector):
    def __init__(self):
        super().__init__(lambda observation_id: self.isin([observation_id]))

    def isin(self, observation_ids):
        """Function to create a condition that selects the observations with the given ids.

        Parameters
        ----------
        observation_ids : list[int]
            Ids of the observations to select (see
            :meth:`~tudatpy.estimation.observations.ObservationDataset.get_observation_ids`); ids that
            do not exist in the dataset are ignored.

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations with the given ids.
        """
        return ObservationSelectionCondition.observation_ids(
            [operator.index(i) for i in observation_ids]
        )


def _link_end(link_end_type):
    return _ComparableSelector(
        lambda link_end_id: ObservationSelectionCondition.link_end(
            link_end_type, _as_link_end_id(link_end_id)
        )
    )


class _LinkEndSelector(_Selector):
    """Selector of a link end of a given type: ``observation_query.link_end(link_end_type) == link_end_id``."""

    def __call__(self, link_end_type):
        return _link_end(link_end_type)

    def __eq__(self, other):
        raise TypeError(
            "observation_query.link_end must first be called with a link end type, e.g. "
            "observation_query.link_end(links.reflector1) == ('Moon', 'Reflector')."
        )

    __ne__ = __eq__
    __hash__ = object.__hash__


class _ObservationQuery(_Selector):
    """Object used to create conditions that select observations in an
    :class:`~tudatpy.estimation.observations.ObservationDataset`.

    The attributes of ``observation_query`` are compared to a value to create an
    :class:`~tudatpy.estimation.observations.ObservationSelectionCondition`, which selects the
    individual observations for which the comparison is true:

    * ``observable_type == ...``, ``set_id == ...`` and ``link_definition == ...``: observations
      with a given observable type, observation set id or link definition (``link_definition``
      selects all observable types with these link ends, for instance both Doppler and range);
    * ``receiver == ...``, ``transmitter == ...``, ``retransmitter == ...`` (also used for
      ``reflector1``), ``receiver2 == ...``, ``transmitter2 == ...``, ``observer == ...`` and
      ``observed_body == ...`` (or ``link_end(link_end_type) == ...`` for other link end types):
      observations with a given link end
      (:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`, or a
      ``(body, station)`` tuple such as ``("Earth", "F51")``);
    * ``time == ...``, ``time < ...``, ``time <= ...``, ``time > ...``, ``time >= ...`` and
      ``time.between(start_time, end_time)``: observations with a given time (in seconds since
      J2000, TDB, or a ``DateTime``/``datetime`` interpreted as TDB). ``time == ...`` requires
      the exact time: for times converted from other time scales, which are stored with a higher
      precision than a float, use ``time.between(t - dt, t + dt)`` or the times returned by
      :meth:`~tudatpy.estimation.observations.ObservationDataset.get_times`;
    * ``active`` and ``rejected``: observations that are not rejected, or that are rejected;
    * ``rejection_reason == ...``: rejected observations with a given rejection reason, such as
      ``"carpino_outlier_rejection"`` or ``"simple_outlier_rejection"`` for the outlier rejection
      during an estimation (active observations have an empty reason, so ``rejection_reason == ""``
      also selects them);
    * ``observation_id == ...`` and ``observation_id.isin(observation_ids)``: observations with
      given ids (see :meth:`~tudatpy.estimation.observations.ObservationDataset.get_observation_ids`);
    * ``residual.abs_greater_than(limit)`` and ``observation.abs_greater_than(limit)``:
      observations with a residual or value larger than ``limit``;
    * ``dependent_variable(settings).greater_than(limit)`` and
      ``dependent_variable(settings).less_than(limit)``: observations with a dependent variable
      larger or smaller than ``limit`` (observations in sets without this dependent variable are
      not selected by either, but are selected by their negation with ``~``);
    * ``all``: all observations.

    ``!=`` creates the opposite condition of ``==``, and ``.isin(values)`` (for instance
    ``set_id.isin([0, 2])`` or ``receiver.isin([("Earth", "DSS-14"), ("Earth", "DSS-63")])``) selects
    the observations equal to any of the values. Conditions are combined with ``&`` (and) and
    ``|`` (or), and negated with ``~``; Python's ``and``, ``or`` and ``not`` cannot be used.
    Comparisons must be put in parentheses when they are combined, for example:

    .. code-block:: python

        from tudatpy.estimation.observable_models_setup import model_settings
        from tudatpy.estimation.observations import observation_query

        condition = (
            observation_query.observable_type == model_settings.angular_position_type
        ) & (observation_query.time >= start_time)
        observation_dataset.reject_observations(
            condition & observation_query.residual.abs_greater_than(1.0e-6)
        )
    """

    observable_type = _ComparableSelector(ObservationSelectionCondition.observable_type)
    set_id = _ComparableSelector(ObservationSelectionCondition.set_id)
    observation_id = _ObservationIdSelector()
    link_definition = _ComparableSelector(ObservationSelectionCondition.link_definition)
    link_end = _LinkEndSelector()
    receiver = _link_end(links.receiver)
    transmitter = _link_end(links.transmitter)
    retransmitter = _link_end(links.retransmitter)
    receiver2 = _link_end(links.receiver2)
    transmitter2 = _link_end(links.transmitter2)
    observer = _link_end(links.observer)
    observed_body = _link_end(links.observed_body)
    time = _TimeSelector()
    all = ObservationSelectionCondition.all()
    active = ObservationSelectionCondition.active()
    rejection_reason = _ComparableSelector(ObservationSelectionCondition.rejection_reason)
    rejected = ObservationSelectionCondition.rejected()
    residual = _VectorValueSelector(
        ObservationSelectionCondition.residual_absolute_value_greater_than
    )
    observation = _VectorValueSelector(
        ObservationSelectionCondition.observation_absolute_value_greater_than
    )

    def dependent_variable(self, settings, return_first_compatible_settings=False):
        """Function to select a dependent variable, to create a condition on its value with
        ``greater_than(limit)`` or ``less_than(limit)``.

        Parameters
        ----------
        settings : ObservationDependentVariableSettings
            Settings of the dependent variable. Observations in sets for which no compatible
            dependent variable is stored are not selected.
        return_first_compatible_settings : bool, default = False
            Whether to use the first compatible dependent variable if a set contains several
            dependent variables that are compatible with ``settings``. If ``False``, an error is
            raised in this case.

        Returns
        -------
        object
            Selector, of which the ``greater_than`` and ``less_than`` functions create the condition.
        """
        return _DependentVariableSelector(settings, return_first_compatible_settings)


observation_query = _ObservationQuery()
