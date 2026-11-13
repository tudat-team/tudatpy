from tudatpy.estimation.observable_models_setup import links
from tudatpy.kernel.estimation.observations import ObservationSelectionCondition


class _Selector:
    """Base class of the selectors of :data:`observation_query`.

    A selector must be compared to a value (or one of its functions must be called) to create a
    condition. Using a selector itself as a boolean raises a ``TypeError``.
    """

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


class _TimeSelector(_ComparableSelector):
    def __init__(self):
        super().__init__(lambda time: self.between(time, time))

    def between(self, start_time, end_time):
        """Function to create a condition that selects the observations with a time between (and
        including) two bounds.

        Parameters
        ----------
        start_time : tudatpy.astro.time_representation.Time
            Lower bound of the observation time, in seconds since J2000 (TDB).
        end_time : tudatpy.astro.time_representation.Time
            Upper bound of the observation time, in seconds since J2000 (TDB).

        Returns
        -------
        ObservationSelectionCondition
            Condition selecting the observations within the time bounds.
        """
        return ObservationSelectionCondition.time_bounds(start_time, end_time)

    def __ge__(self, time):
        return ObservationSelectionCondition.time_greater_equal(time)

    def __gt__(self, time):
        return ObservationSelectionCondition.time_greater_than(time)

    def __le__(self, time):
        return ObservationSelectionCondition.time_less_equal(time)

    def __lt__(self, time):
        return ObservationSelectionCondition.time_less_than(time)


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


def _link_end(link_end_type):
    return _ComparableSelector(
        lambda link_end_id: ObservationSelectionCondition.link_end(link_end_type, link_end_id)
    )


class _ObservationQuery(_Selector):
    """Object used to create conditions that select observations in an
    :class:`~tudatpy.estimation.observations.ObservationDataset`.

    The attributes of ``observation_query`` are compared to a value to create an
    :class:`~tudatpy.estimation.observations.ObservationSelectionCondition`, which selects the
    individual observations for which the comparison is true:

    * ``observable_type == ...``, ``set_id == ...`` and ``link_definition == ...``: observations
      with a given observable type, observation set id or link definition;
    * ``receiver == ...``, ``transmitter == ...``, ``observer == ...`` and
      ``observed_body == ...`` (or ``link_end(link_end_type) == ...`` for other link end types):
      observations with a given link end
      (:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`);
    * ``time == ...``, ``time < ...``, ``time <= ...``, ``time > ...``, ``time >= ...`` and
      ``time.between(start_time, end_time)``: observations with a given time (in seconds since
      J2000, TDB);
    * ``active`` and ``rejected``: observations that are not rejected, or that are rejected;
    * ``residual.abs_greater_than(limit)`` and ``observation.abs_greater_than(limit)``:
      observations with a residual or value larger than ``limit``;
    * ``dependent_variable(settings).greater_than(limit)``: observations with a dependent variable
      larger than ``limit``.

    ``!=`` creates the opposite condition of ``==``. Conditions are combined with ``&`` (and) and
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
    link_definition = _ComparableSelector(ObservationSelectionCondition.link_definition)
    link_end = staticmethod(_link_end)
    receiver = _link_end(links.receiver)
    transmitter = _link_end(links.transmitter)
    observer = _link_end(links.observer)
    observed_body = _link_end(links.observed_body)
    time = _TimeSelector()
    active = ObservationSelectionCondition.active()
    rejected = ObservationSelectionCondition.rejected()
    residual = _VectorValueSelector(
        ObservationSelectionCondition.residual_absolute_value_greater_than
    )
    observation = _VectorValueSelector(
        ObservationSelectionCondition.observation_absolute_value_greater_than
    )

    def dependent_variable(self, settings, return_first_compatible_settings=False):
        """Function to select a dependent variable, to create a condition on its value with
        ``greater_than(limit)``.

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
            Selector, of which the ``greater_than`` function creates the condition.
        """
        return _DependentVariableSelector(settings, return_first_compatible_settings)


observation_query = _ObservationQuery()
