from ....estimation.observable_models_setup import biases
from ....estimation.observable_models_setup import light_time_corrections
from ....estimation.observable_models_setup import links
import typing
__all__ = ['DopplerProperTimeRateSettings', 'NWayRangeObservationSettings', 'ObservableType', 'ObservationSettings', 'OneWayDopplerObservationSettings', 'angular_position', 'angular_position_type', 'cartesian_position', 'cartesian_velocity', 'doppler_measured_frequency', 'doppler_measured_frequency_type', 'dsn_n_way_Range', 'dsn_n_way_averaged_doppler', 'dsn_n_way_doppler_averaged', 'dsn_n_way_doppler_averaged_from_one_way_links', 'dsn_n_way_range', 'dsn_one_way_averaged_doppler', 'euler_angle_313_observable_type', 'euler_angles_313', 'n_way_averaged_doppler_type', 'n_way_doppler_averaged', 'n_way_doppler_averaged_from_one_way_links', 'n_way_range', 'n_way_range_from_one_way_links', 'n_way_range_type', 'one_way_averaged_doppler_type', 'one_way_closed_loop_doppler', 'one_way_doppler_averaged', 'one_way_doppler_instantaneous', 'one_way_instantaneous_doppler_type', 'one_way_open_loop_doppler', 'one_way_range', 'one_way_range_type', 'position_observable_type', 'relative_angular_position', 'relative_angular_position_type', 'relative_cartesian_position', 'relative_position_observable_type', 'two_doppler_instantaneous', 'two_way_doppler_averaged', 'two_way_doppler_averaged_from_one_way_links', 'two_way_doppler_instantaneous_from_one_way_links', 'two_way_instantaneous_doppler_type', 'two_way_open_loop_doppler', 'two_way_open_loop_doppler_from_one_way_links', 'two_way_range', 'two_way_range_from_one_way_links', 'velocity_observable_type']

class DopplerProperTimeRateSettings:
    """Base class to define proper time rate settings.
    
    Base class to define proper time rate settings (at a single link end) for instantaneous Doppler observation model settings.
    Specific proper time rate settings must be defined using an object derived from this class.
    The derived classes are made accessible via dedicated functions."""

class NWayRangeObservationSettings(ObservationSettings):
    """No documentation found."""

class ObservableType:
    """Enumeration of available observable types.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to print all available Observable Types
        from tudatpy.estimation.observable_models_setup import model_settings
    
        num_observable_types = len(model_settings.ObservableType.__members__)
        print(f'The length of all available Tudatpy Observable Types is: {num_observable_types}')
    
        # Print all available Observable Types using the "name" property
        for i in range(num_observable_types):
            print(i, model_settings.ObservableType(i).name)
    
    
    
    
          
    
    Members:
    
      one_way_range_type
    
      n_way_range_type
    
      angular_position_type
    
      relative_angular_position_type
    
      position_observable_type
    
      velocity_observable_type
    
      relative_position_observable_type
    
      one_way_instantaneous_doppler_type
    
      one_way_averaged_doppler_type
    
      two_way_instantaneous_doppler_type
    
      n_way_averaged_doppler_type
    
      euler_angle_313_observable_type
    
      dsn_one_way_averaged_doppler
    
      dsn_n_way_averaged_doppler
    
      doppler_measured_frequency_type
    
      dsn_n_way_range"""
    __members__: typing.ClassVar[dict[str, ObservableType]]
    angular_position_type: typing.ClassVar[ObservableType]
    doppler_measured_frequency_type: typing.ClassVar[ObservableType]
    dsn_n_way_averaged_doppler: typing.ClassVar[ObservableType]
    dsn_n_way_range: typing.ClassVar[ObservableType]
    dsn_one_way_averaged_doppler: typing.ClassVar[ObservableType]
    euler_angle_313_observable_type: typing.ClassVar[ObservableType]
    n_way_averaged_doppler_type: typing.ClassVar[ObservableType]
    n_way_range_type: typing.ClassVar[ObservableType]
    one_way_averaged_doppler_type: typing.ClassVar[ObservableType]
    one_way_instantaneous_doppler_type: typing.ClassVar[ObservableType]
    one_way_range_type: typing.ClassVar[ObservableType]
    position_observable_type: typing.ClassVar[ObservableType]
    relative_angular_position_type: typing.ClassVar[ObservableType]
    relative_position_observable_type: typing.ClassVar[ObservableType]
    two_way_instantaneous_doppler_type: typing.ClassVar[ObservableType]
    velocity_observable_type: typing.ClassVar[ObservableType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class ObservationSettings:
    """Base class to define settings of observation models.
    
    Base class to define settings of observation models.
    Observation model settings define at least the type and geometry of a given observation.
    They can furthermore set observation biases and/or light-time corrections.
    Simple observation models settings that are fully characterised by these elements can be managed by this base class.
    Instances of this class are typically created via functions, such as
    :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`, :func:`~tudatpy.estimation.observable_models_setup.model_settings.cartesian_position`,
    :func:`~tudatpy.estimation.observable_models_setup.model_settings.angular_position`, etc.
    Model settings for specific observation models that require additional information such as integration time, retransmission time, etc. must be defined using an object derived from this class.
    The derived classes are made accessible through further functions.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of an ObservationSettings object
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
        # Other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) are set by default
        observation_settings = model_settings.one_way_range(Link_Definition_Object)
    
        # Show that it is an ObservationSettings object.
        print(observation_settings)"""

class OneWayDopplerObservationSettings(ObservationSettings):
    """Derived Class for defining the settings of one-way instantaneous Doppler observation models.
    
    Derived Class for defining the settings of one-way instantaneous Doppler observation models.
    Settings object can account for additional observation model aspects such as light time corrections and proper time rate settings.
    Instances of this class can be created via the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_instantaneous` function.
    Associated base class: :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of a OneWayDopplerObservationSettings object
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # Use: model_settings.one_way_doppler_instantaneous to create a OneWayDopplerObservationSettings object (only required Link_Definition_Object argument is passed)
        # Other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings, proper time rate) are set by default
        doppler_observation_settings = model_settings.one_way_doppler_instantaneous(Link_Definition_Object)
    
        # Show that it is an OneWayDopplerObservationSettings object.
        print(doppler_observation_settings)"""

def angular_position(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for an angular position observable.
    
    Function for creating observation model settings of angular position type observables (as right ascension :math:`\\alpha` and declination :math:`\\delta`),
    for a single link definition. The associated observation model creates an observable :math:`\\mathbf{h}_{_{\\text{ang.pos.}}}` of type two as follows (in the unbiased case):
    
    .. math::
       \\Delta\\mathbf{r}=\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T}(t_{T})\\\\
       \\tan\\alpha=\\frac{\\Delta r_{y}}{\\Delta r_{x}}\\\\
       \\delta=\\frac{\\Delta r_{z}}{\\sqrt{\\Delta r_{x}^{2}+\\Delta r_{y}^{2}}}\\\\
       \\mathbf{h}_{_{\\text{ang.pos.}}} = [\\alpha;\\delta]
    
    The relative position vector :math:`\\Delta\\mathbf{r}` is computed identically as described for the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`
    The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
    environment.
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` class defining the settings for the angular position observable."""

def cartesian_position(link_ends: links.LinkDefinition, bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for a Cartesian position observable.
    
    Function for creating observation model settings of Cartesian position type observables.
    Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
    This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
    The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires that the
        `observed_body` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    Returns
    -------
    :class:`ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` class defining the settings for the cartesian position observable."""

def cartesian_velocity(link_ends: links.LinkDefinition, bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for a Cartesian velocity observable.
    
    Function for creating observation model settings of Cartesian position type observables.
    Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
    This observable provides the inertial (w.r.t. global frame origin) Cartesian velocity of the `observed_body` defined by the `link_ends` input.
    The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` velocity
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires that the
        `observed_body` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    Returns
    -------
    :class:`ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` class defining the settings for the cartesian velocity observable."""

def doppler_measured_frequency(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """No documentation found."""

def dsn_n_way_Range(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """No documentation found."""

def dsn_n_way_doppler_averaged(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=..., subtract_doppler_signature: bool=True) -> ObservationSettings:
    """No documentation found."""

def dsn_n_way_doppler_averaged_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=..., subtract_doppler_signature: bool=True) -> ObservationSettings:
    """No documentation found."""

def euler_angles_313(link_ends: links.LinkDefinition, bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """No documentation found."""

def n_way_doppler_averaged(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for an n-way averaged Doppler observable.
    
    Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implementation is
    analogous to the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` observable. But, in the present case
    the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\\Delta t`.
    
    The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
        as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
    
    light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable."""

def n_way_doppler_averaged_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for an n-way averaged Doppler observable.
    
    Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition.
    The implementation is the same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_doppler_averaged`, with the difference
    that the constituent one-way range observables may have different settings.
    
    
    Parameters
    ----------
    one_way_range_settings : List[ :class:`ObservationSettings` ]
        List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way averaged range rate observable.
        The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
        ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` (n-1) and ``receiver`` are defined by the
        ``transmitter`` and ``receiver`` of the :math:`\\text{n}^{th}` entry of this list.
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable."""

def n_way_range(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for a n-way range observable.
    
    Function for creating observation model settings of n-way range type observables, for a single link definition. The associated observation model creates
    a single-valued observable :math:`h_{_{\\text{N-range}}}` by combining together a series :math:`n` one-way range observations
    (see :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`). By default, the reception time of the :math:`i^{th}` one-way range is set as the
    transmission time of the :math:`(i+1)^{th}` one-way range. A retransmission delay may be defined by ancilliary settings (see :func:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings`) when creating observation
    simulation setings.
    
    For this function, the settings for each constituent one-way range (with the exception of the link end identifiers) are equal.
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
        as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
    
    light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
        Note that only one bias setting is applied to the n-way observable.
    
    light_time_convergence_settings : :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings` class.
    
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the n_way_range function
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # n_way_range() takes a Link Definition Object as input to the function.
        # Let's create it from link_ends
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
        # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
        observation_settings = model_settings.n_way_range(Link_Definition_Object)
    
        # Show that n_way_range() returns an NWayRangeObservationSettings object.
        print(observation_settings)"""

def n_way_range_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for a n-way range observable.
    
    Function for creating observation model settings of n-way range type observables, for a single link definition. The
    implementation is the same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_range`, with the difference
    that the constituent one-way ranges may have different settings.s
    
    
    Parameters
    ----------
    one_way_range_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` ]
        List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way range observable.
        The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
        ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` (n-1) and ``receiver`` are defined by the
        ``transmitter`` and ``receiver`` of the :math:`\\text{n}^{th}` entry of this list.
    
    bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
        Note that only one bias setting is applied to the n-way observable.
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings` class.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the n_way_range_from_one_way_links function
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # n_way_range_from_one_way_links() takes 1) a list of ObservationSettings objects and 2) bias as input (default is None)
        # Let's create it.
        Link_Definition_Object = links.LinkDefinition(link_ends) # define LinkDefinition object
        n_way_observation_settings_list = [model_settings.n_way_range(Link_Definition_Object)] # define (minimal) ObservationSettings object
    
        n_way_from_one_link_observation_settings = model_settings.n_way_range_from_one_way_links(n_way_observation_settings_list, bias_settings = None)
    
        # Show that n_way_range_from_one_way_links() returns an NWayRangeObservationSettings object.
        print(n_way_from_one_link_observation_settings)"""

@typing.overload
def one_way_closed_loop_doppler(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    ...

@typing.overload
def one_way_closed_loop_doppler(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    ...

def one_way_doppler_averaged(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for a one-way averaged Doppler observable.
    
    Function for creating observation model settings for one-way averaged Doppler observables, for a single link definition. The associated observation model creates
    a single-valued observable :math:`h_{_{\\text{1-\\bar{Dopp}}}}` as follows (in the unbiased case):
    
    .. math::
       h_{_{\\text{1-\\bar{Dopp}}}}&=c\\int_{t-\\Delta t}^{t+\\Delta t}\\frac{t_{T}}{dt_{R}}d\\bar{t}\\\\
                                 &=\\frac{h_{_{\\text{1-range}}}(t_{R}=t+\\Delta t,t_{T})-h_{_{\\text{1-range}}}(t_{R}=t,t_{T})}{\\Delta t}
    
    where, in the latter formulation (which is the one that is implemented), the observable is referenced to the receiver time. This averaged Doppler observable
    is computed as the difference of two one-way range observables (see :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`),
    with the reference time shifted by :math:`\\Delta t`. As such, it is sensitive to numerical errors for small :math:`\\Delta t`
    
    The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires that the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived `OneWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable."""

def one_way_doppler_instantaneous(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, transmitter_proper_time_rate_settings: DopplerProperTimeRateSettings=None, receiver_proper_time_rate_settings: DopplerProperTimeRateSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    """Function for creating settings for a one-way instantaneous Doppler observable.
    
    Function for creating settings for a one-way instantaneous Doppler observable for a single link definition. The associated observation model creates
    a single-valued observable :math:`h_{_{\\text{1-Dopp.}}}` as follows (in the unbiased case):
    
    .. math::
       h_{_{\\text{1-Dopp.}}}=c\\left(\\frac{d\\tau_{T}}{dt_{T}}\\frac{t_{T}}{dt_{R}}\\frac{dt_{R}}{d\\tau_{R}}-1\\right)
    
    where :math:`t` and :math:`\\tau` denote coordinate and proper time of the transmitter T and receiver R, respectively.
    The receiver and transmitter position and coordinate time are computed identically as described for the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`.
    The detailed mathematical implementation are described in: `Moyer, T.D. (2000) Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation. Monograph 2, Deep Space Communications and Navigation Series, JPL Publication 00-7 <https://www.scirp.org/reference/referencespapers?referenceid=1827210>`_.
    
    This observable represents the 'instantaneous (non-integrated)' Doppler observable, as obtained from open-loop observations.
    It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the
    :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` should be used)
    
    The coordinate
    time derivative :math:`\\frac{t_{A}}{dt_{B}}` is always computed when generating this observable. Settings for the proper time
    rates :math:`\\frac{d\\tau}{dt}` can be specified by the user through the ``transmitter_proper_time_rate_settings`` and ``receiver_proper_time_rate_settings``
    arguments (inputs, see Parameters). Whenever these are left empty, the proper time rates are omitted (set to 1.0).
    
    The observable may be non-dimensionalized by the speed of light :math:`c`, which results in the observable being equal to thee received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires that the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    transmitter_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
        Settings for computing the transmitter proper time rate :math:`\\frac{d\\tau}{dt}`, default is none (:math:`\\frac{d\\tau}{dt}=1`)
    
    receiver_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
        Settings for computing the receiver proper time rate :math:`\\frac{d\\tau}{dt}`, default is none (:math:`\\frac{d\\tau}{dt}=1`)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    normalized_with_speed_of_light : bool, default = false
        Option to non-dimensionalize the observable with speed of light :math:`c`
    
    Returns
    -------
    :class:`OneWayDopplerObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived :class:`OneWayDopplerObservationSettings` class defining the settings for the one-way open doppler observable observable."""

def one_way_open_loop_doppler(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, transmitter_proper_time_rate_settings: DopplerProperTimeRateSettings=None, receiver_proper_time_rate_settings: DopplerProperTimeRateSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    ...

def one_way_range(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for a one-way range observable.
    
    Function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
    a single-valued observable :math:`h_{_{\\text{1-range}}}` as follows (in the unbiased case):
    
    .. math::
       h_{_{\\text{1-range}}}(t_{R},t_{T})=|\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T}(t_{T})| + \\Delta s
    
    where :math:`\\mathbf{r}_{R}`, :math:`\\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
    of receiver and transmitter. The term :math:`\\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
    The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
    As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:
    
    .. math::
       t_{R}-t_{T}=c\\left(|\\mathbf{r}_{R}(t_{R})-\\mathbf{r}(t_{R})| + \\Delta s\\right)
    
    The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)
    
    light_time_convergence_settings : :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` class defining the settings for the one-way observable.
    
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the one_way_range function
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary. This will be the input to the function.
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
        # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
        observation_settings = model_settings.one_way_range(Link_Definition_Object)
    
        # Show that this returns an ObservationSettings object.
        print(observation_settings)"""

def relative_angular_position(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for a relative angular position observable.
    
    Function for creating observation model settings of relative angular position type observables (as right ascension difference :math:`\\Delta \\alpha` and declination difference :math:`\\Delta \\delta`),
    from two link definitions. The associated observation model creates an observable :math:`\\mathbf{h}_{_{\\text{ang.pos.}}}` of type two as follows (in the unbiased case):
    
    .. math::
    
       \\Delta\\mathbf{r}_1=\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T1}(t_{T1})             \\\\
       \\tan\\alpha_{1} =\\frac{\\Delta r_{1y}}{\\Delta r_{1x}}                          \\\\
       \\delta_{1} =\\frac{\\Delta r_{1z}}{\\sqrt{\\Delta r_{1x}^{2}+\\Delta r_{1y}^{2}}} \\\\
       \\Delta\\mathbf{r}_2=\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T2}(t_{T2})             \\\\
       \\tan\\alpha_{2} =\\frac{\\Delta r_{2y}}{\\Delta r_{2x}}                          \\\\
       \\delta_{2} =\\frac{\\Delta r_{2z}}{\\sqrt{\\Delta r_{2x}^{2}+\\Delta r_{2y}^{2}}} \\\\
       \\mathbf{h}_{_{\\text{rel.ang.pos.}}} = [\\alpha_{2}-\\alpha_{1};\\delta_{2}-\\delta_{1}]
    
    The relative position vectors :math:`\\Delta\\mathbf{r}_1` and :math:`\\Delta\\mathbf{r}_2` are computed identically as described for the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`
    The relative angular position observable can be used for optical astrometry, optical navigation, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
    environment.
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter`, `transmitter2` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` class defining the settings for the relative angular position observable."""

def relative_cartesian_position(link_ends: links.LinkDefinition, bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for a relative Cartesian position observable.
    
    Function for creating observation model settings of relative Cartesian position type observables.
    Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
    This observable provides the inertial Cartesian position of the `observed_body`, w.r.t. the `observer` defined by the `link_ends` input.
    The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires that the
        `observed_body` and `observer` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    Returns
    -------
    :class:`ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` class defining the settings for the relative Cartesian position observable."""

def two_doppler_instantaneous(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    """No documentation found."""

def two_way_doppler_averaged(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for an n-way averaged Doppler observable.
    
    Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implementation is
    analogous to the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` observable. But, in the present case
    the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\\Delta t`.
    
    The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
        as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
    
    light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable."""

def two_way_doppler_averaged_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for an n-way averaged Doppler observable.
    
    Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
    analogous to the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` observable. But, in the present case
    the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\\Delta t`.
    
    The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
        as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
    
    light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable."""

def two_way_doppler_instantaneous_from_one_way_links(uplink_doppler_settings: OneWayDopplerObservationSettings, downlink_doppler_settings: OneWayDopplerObservationSettings, bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for a two-way instantaneous Doppler observable.
    
    
    Function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The
    implementation is the same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous`, with the difference
    that the constituent one-way ranges may have different settings.
    
    The observable may be non-dimensionalized by the speed of light :math:`c` (in the constituent one-way Doppler observable settings),
    which results in the observable being equal to the received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.
    
    
    Parameters
    ----------
    uplink_doppler_settings : :class:`OneWayDopplerObservationSettings`
        Settings for uplink leg of one-way observable, created using :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_open_loop_doppler`
    
    downlink_doppler_settings : :class:`OneWayDopplerObservationSettings`
        Settings for downlink leg of one-way observable, created using :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_open_loop_doppler`
    
    bias_settings : :class:`ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the full observation, default is none (unbiased observation). Note that,
        even if no bias is applied to the two-way observable, the constituent one-way observables may still be biased.
    
    light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`TwoWayDopplerObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived :class:`TwoWayDopplerObservationSettings` class defining the settings for the two-way open doppler observable."""

def two_way_open_loop_doppler(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    ...

def two_way_open_loop_doppler_from_one_way_links(uplink_doppler_settings: OneWayDopplerObservationSettings, downlink_doppler_settings: OneWayDopplerObservationSettings, bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    ...

def two_way_range(link_ends: links.LinkDefinition, light_time_correction_settings: list[light_time_corrections.LightTimeCorrectionSettings]=[], bias_settings: biases.ObservationBiasSettings=None, light_time_convergence_settings: light_time_corrections.LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Function for creating settings for a two-way range observable.
    
    Same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_range`, with :math:`n=2`. This function is provided
    for convenience.
    
    
    Parameters
    ----------
    link_ends : LinkDefinition
        Set of link ends that define the geometry of the observation. This observable requires the
        `transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined
    
    light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
        List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
        in the signal being modelled as moving in a straight line with the speed of light
    
    bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
        Note that only one bias setting is applied to the n-way observable.
    
    light_time_convergence_settings : :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
        Settings for convergence of the light-time
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings` class.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the two_way_range function
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # two_way_range() takes a Link Definition Object as input to the function.
        # Note: as for this case, transmitter, retransmitter and receiver are required to define the Link Ends dictionary
        link_ends = dict()
        link_ends[links.transmitter] = links.body_origin_link_end_id("Earth")
        link_ends[links.retransmitter] = links.body_origin_link_end_id("Delfi-C3")
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
    
        # Create the LinkDefinition object
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
        # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
        observation_settings = model_settings.two_way_range(Link_Definition_Object)
    
        # Show that two_way_range() returns an NWayRangeObservationSettings object.
        print(observation_settings)"""

def two_way_range_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: biases.ObservationBiasSettings=None) -> ObservationSettings:
    """Function for creating settings for a two-way range observable.
    
    Same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_range_from_one_way_links`, with :math:`n=2`. This function is provided
    for convenience.
    
    
    Parameters
    ----------
    one_way_range_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` ]
        List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
        The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
        ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
        ``transmitter`` and ``receiver`` of the second entry of this list.
    
    bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
        Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
        Note that only one bias setting is applied to the n-way observable.
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationSettings` class.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the two_way_range_from_one_way_links function
        from tudatpy.estimation.observable_models_setup import links, model_settings
    
        # two_way_range_from_one_way_links() takes a list of ObservationSettings objects
        # Note: as for this case, transmitter, retransmitter and receiver are required to define the Link Ends dictionary
        link_ends = dict()
        link_ends[links.transmitter] = links.body_origin_link_end_id("Earth")
        link_ends[links.retransmitter] = links.body_origin_link_end_id("Delfi-C3")
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
    
        # Create the LinkDefinition object to be used as input
        Link_Definition_Object = links.LinkDefinition(link_ends) # define LinkDefinition object
        two_way_range_observation_settings_list = [model_settings.two_way_range(Link_Definition_Object)] # define (minimal) NWayRangeObservationSettings object
        two_way_range_one_way_link_settings = model_settings.two_way_range_from_one_way_links(two_way_range_observation_settings_list)
    
        # Show that two_way_range_from_one_way_links() returns an NWayRangeObservationSettings object.
        print(two_way_range_one_way_link_settings)"""
angular_position_type: ObservableType
doppler_measured_frequency_type: ObservableType
dsn_n_way_averaged_doppler: ObservableType
dsn_n_way_range: ObservableType
dsn_one_way_averaged_doppler: ObservableType
euler_angle_313_observable_type: ObservableType
n_way_averaged_doppler_type: ObservableType
n_way_range_type: ObservableType
one_way_averaged_doppler_type: ObservableType
one_way_instantaneous_doppler_type: ObservableType
one_way_range_type: ObservableType
position_observable_type: ObservableType
relative_angular_position_type: ObservableType
relative_position_observable_type: ObservableType
two_way_instantaneous_doppler_type: ObservableType
velocity_observable_type: ObservableType