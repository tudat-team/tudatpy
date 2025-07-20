from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
import typing
__all__ = ['FrequencyBands', 'ObservationAncilliarySimulationSettings', 'ObservationAncilliarySimulationVariable', 'ObservationIntermediateSimulationVariable', 'add_ancilliary_settings_to_observable', 'add_ancilliary_settings_to_observable_for_link_ends', 'cassini_turnaround_ratios', 'doppler_ancilliary_settings', 'doppler_integration_time', 'doppler_measured_frequency_ancillary_settings', 'doppler_reference_frequency', 'dsn_default_turnaround_ratios', 'dsn_n_way_doppler_ancilliary_settings', 'dsn_n_way_range_ancilliary_settings', 'frequency_bands', 'link_ends_delays', 'n_way_doppler_ancilliary_settings', 'n_way_range_ancilliary_settings', 'range_conversion_factor', 'received_frequency_intermediate', 'reception_reference_frequency_band', 'sequential_range_lowest_ranging_component', 'transmitter_frequency_intermediate', 'two_way_doppler_ancilliary_settings', 'two_way_range_ancilliary_settings']

class FrequencyBands:
    """No documentation found.
    
    Members:
    
      s_band
    
      x_band
    
      ka_band
    
      ku_band"""
    __members__: typing.ClassVar[dict[str, FrequencyBands]]
    ka_band: typing.ClassVar[FrequencyBands]
    ku_band: typing.ClassVar[FrequencyBands]
    s_band: typing.ClassVar[FrequencyBands]
    x_band: typing.ClassVar[FrequencyBands]

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

class ObservationAncilliarySimulationSettings:
    """Base class for holding ancilliary settings for observation simulation.
    
    Base class for holding ancilliary settings for observation simulation.
    The user can create instances of this class via the :func:`~estimation.observations_setup.observations_dependent_variables.elevation_angle_dependent_variable` function.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of an ObservationAncillarySimulationSettings object
        from tudatpy.estimation.observations_setup import ancillary_settings
    
        # Example 1: Create ObservationAncillarySimulationSettings object using ancillary_settings.n_way_range_ancilliary_settings function
        # In this case the frequency bands of the retransmitter - we set it to x band.
        n_way_range_ancillary_settings = ancillary_settings.n_way_range_ancilliary_settings(frequency_bands=[ancillary_settings.FrequencyBands.x_band])
    
        # Show that this is indeed an ObservationAncillarySimulationSettings object
        print(n_way_range_ancillary_settings)
    
        # Example 2: Create ObservationAncillarySimulationSettings object using ancillary_settings.doppler_ancilliary_settings function
        # In this case the integration time (in seconds) has to be given as input - we set it to 60s
        doppler_ancillary_settings = ancillary_settings.doppler_ancilliary_settings(60)
    
        # Show that this is indeed an ObservationAncillarySimulationSettings object
        print(doppler_ancillary_settings)
    
        # [OPTIONAL] Verify that we indeed added Frequency Bands as Ancillary Simulation Variables for the n_way_range_ancillary_settings.
        list_num = n_way_range_ancillary_settings.get_float_list_settings(ancillary_settings.ObservationAncilliarySimulationVariable.frequency_bands)
        for num in list_num:
            name = ancillary_settings.ObservationAncilliarySimulationVariable(int(num)).name
            print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncilliarySimulationVariable Enumeration')"""

    def get_float_list_settings(self, setting_type: ObservationAncilliarySimulationVariable, throw_exception: bool=True) -> list[float]:
        """
                 Parameters
                 ----------
                 setting_type : ObservationAncilliarySimulationVariable
                     Type of the setting for which the value is to be returned
        
                 throw_exception : bool, default = false
                     Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.
        
                 Returns
                 -------
                 list[ float ]
                     Value of the requested ancilliary variable (or empty list if it does not exist)
        
                 Examples
                 --------
                 .. code-block:: python
        
                     # Code snippet to show how to retrieve ObservationAncillarySimulationSettings info
                     # using the ObservationAncilliarySimulationSettings.get_float_settings function
        
                     from tudatpy.estimation.observations_setup import ancillary_settings
        
                     # Create Ancillary Settings
                     n_way_range_ancillary_settings = ancillary_settings.n_way_range_ancilliary_settings(frequency_bands=[ancillary_settings.FrequencyBands.x_band])
        
                     # Verify that we indeed added Frequency Bands as Ancillary Simulation Variables, using n_way_range_ancillary_settings.get_float_list_settings
                     list_num = n_way_range_ancillary_settings.get_float_list_settings(ancillary_settings.ObservationAncilliarySimulationVariable.frequency_bands)
                     for num in list_num:
                         name = ancillary_settings.ObservationAncilliarySimulationVariable(int(num)).name
                         print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncilliarySimulationVariable Enumeration')
        """

    def get_float_settings(self, setting_type: ObservationAncilliarySimulationVariable, throw_exception: bool=True) -> float:
        """
                 Parameters
                 ----------
                 setting_type : ObservationAncilliarySimulationVariable
                     Type of the setting for which the value is to be returned
        
                 throw_exception : bool, default = false
                     Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.
        
                 Returns
                 -------
                 float
                     Value of the requested ancilliary variable (or NaN if it does not exist)
        """

    def set_intermediate_double_data(self, variable: ..., value: float) -> None:
        ...

class ObservationAncilliarySimulationVariable:
    """Enumeration of observation ancillary variable types.
    
    Examples
    --------
    
    .. code-block:: python
    
        # Code snippet to print all available Observation Ancillary Variable Types
        from tudatpy.estimation.observations_setup import ancillary_settings
    
        num_observation_ancillary_variable_types = len(ancillary_settings.ObservationAncilliarySimulationVariable.__members__)
        print(f'The length of all available Tudatpy  Observation Ancillary Variable Types is: {num_observation_ancillary_variable_types}')
    
        # Print all Observation Ancillary Variable Types using the "name" property
        for i in range(num_observation_ancillary_variable_types):
            print(i, ancillary_settings.ObservationAncilliarySimulationVariable(i).name)
    
    
    
          
    
    Members:
    
      link_ends_delays
    
      doppler_integration_time
    
      doppler_reference_frequency
    
      frequency_bands
    
      reception_reference_frequency_band
    
      sequential_range_lowest_ranging_component
    
      range_conversion_factor"""
    __members__: typing.ClassVar[dict[str, ObservationAncilliarySimulationVariable]]
    doppler_integration_time: typing.ClassVar[ObservationAncilliarySimulationVariable]
    doppler_reference_frequency: typing.ClassVar[ObservationAncilliarySimulationVariable]
    frequency_bands: typing.ClassVar[ObservationAncilliarySimulationVariable]
    link_ends_delays: typing.ClassVar[ObservationAncilliarySimulationVariable]
    range_conversion_factor: typing.ClassVar[ObservationAncilliarySimulationVariable]
    reception_reference_frequency_band: typing.ClassVar[ObservationAncilliarySimulationVariable]
    sequential_range_lowest_ranging_component: typing.ClassVar[ObservationAncilliarySimulationVariable]

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

class ObservationIntermediateSimulationVariable:
    """Members:
    
      transmitter_frequency_intermediate
    
      received_frequency_intermediate"""
    __members__: typing.ClassVar[dict[str, ObservationIntermediateSimulationVariable]]
    received_frequency_intermediate: typing.ClassVar[ObservationIntermediateSimulationVariable]
    transmitter_frequency_intermediate: typing.ClassVar[ObservationIntermediateSimulationVariable]

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

def add_ancilliary_settings_to_observable(observation_simulation_settings_list: list[...], ancilliary_settings: ObservationAncilliarySimulationSettings, observable_type: model_settings.ObservableType) -> None:
    """No documentation found."""

def add_ancilliary_settings_to_observable_for_link_ends(observation_simulation_settings_list: list[...], ancilliary_settings: ObservationAncilliarySimulationSettings, observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition) -> None:
    """No documentation found."""

def cassini_turnaround_ratios(uplink_band: FrequencyBands, downlink_band: FrequencyBands) -> float:
    """No documentation found."""

def doppler_ancilliary_settings(integration_time: float=60.0) -> ObservationAncilliarySimulationSettings:
    """Function for creating ancilliary settings for averaged Doppler observable.
    
    Function for creating ancilliary settings for an averaged Doppler observable. Specifically, this
    function can be used to create settings for the integration time of the observable. Note: in case no retransmission
    delays (or other additional ancilliary settings) are to be defined, this setting may be used for one-, two-, or N-way
    averaged Doppler.
    
    
    Parameters
    ----------
    integration_time : float, default = 60.0
        Integration time that is to be used for the averaged Doppler observable
    Returns
    -------
    :class:`ObservationAncilliarySimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings."""

def doppler_measured_frequency_ancillary_settings(frequency_bands: list[...]) -> ObservationAncilliarySimulationSettings:
    """No documentation found."""

def dsn_default_turnaround_ratios(uplink_band: FrequencyBands, downlink_band: FrequencyBands) -> float:
    """No documentation found."""

def dsn_n_way_doppler_ancilliary_settings(frequency_bands: list[...], reference_frequency_band: ..., reference_frequency: float, integration_time: float=60.0, link_end_delays: list[float]=[]) -> ObservationAncilliarySimulationSettings:
    """No documentation found."""

def dsn_n_way_range_ancilliary_settings(frequency_bands: list[...], lowest_ranging_component: float, link_end_delays: list[float]=[]) -> ObservationAncilliarySimulationSettings:
    """No documentation found."""

def n_way_doppler_ancilliary_settings(integration_time: float=60.0, link_end_delays: list[float]=[], frequency_bands: list[...]=[]) -> ObservationAncilliarySimulationSettings:
    """Function for creating ancilliary settings for n-way averaged Doppler observable.
    
    Function for creating ancilliary settings for a n-way averaged Doppler observable. Specifically, this
    function can be used to create settings for the integration time of the observable, and the  retransmission delays for each of the retransmitters.
    
    
    Parameters
    ----------
    integration_time : float, default = 60.0
        Integration time that is to be used for the averaged Doppler observable
    retransmission_delays : list[ float ], default = None
        Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
    Returns
    -------
    :class:`ObservationAncilliarySimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings."""

def n_way_range_ancilliary_settings(link_end_delays: list[float]=[], frequency_bands: list[...]=[]) -> ObservationAncilliarySimulationSettings:
    """Function for creating ancilliary settings for n-way range observable.
    
    Function for creating ancilliary settings for a n-way range observable. Specifically, this
    function can be used to create settings for the retransmission delays of the observable, for each of the retransmitters.
    
    
    Parameters
    ----------
    retransmission_delays : list[ float ], default = None
        Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
    Returns
    -------
    :class:`ObservationAncilliarySimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings."""

def two_way_doppler_ancilliary_settings(integration_time: float=60.0, retransmission_delay: float=0.0) -> ObservationAncilliarySimulationSettings:
    """Function for creating ancilliary settings for two-way averaged Doppler observable.
    
    Function for creating ancilliary settings for a two-way range observable. Specifically, this
    function can be used to create settings for the retransmission delay of the observable.  NOTE:
    this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.estimation.observations_setup.ancillary_settings.n_way_doppler_ancilliary_settings`
    with a single retransmission delay.
    
    
    Parameters
    ----------
    integration_time : float, default = 60.0
        Integration time that is to be used for the averaged Doppler observable
    retransmission_delay : float, default = 0.0
        Retransmission delay that is to be applied to the simulation of the two-way observable
    Returns
    -------
    :class:`ObservationAncilliarySimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings."""

def two_way_range_ancilliary_settings(retransmission_delay: float=0.0) -> ObservationAncilliarySimulationSettings:
    """Function for creating ancilliary settings for two-way range observable.
    
    Function for creating ancilliary settings for a two-way range observable. Specifically, this
    function can be used to create settings for the retransmission delay of the observable. NOTE:
    this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.estimation.observations_setup.ancillary_settings.n_way_range_ancilliary_settings`
    with a single retransmission delay.
    
    
    Parameters
    ----------
    retransmission_delay : float, default = 0.0
        Retransmission delay that is to be applied to the simulation of the two-way observable
    Returns
    -------
    :class:`ObservationAncilliarySimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings` with the required settings."""
doppler_integration_time: ObservationAncilliarySimulationVariable
doppler_reference_frequency: ObservationAncilliarySimulationVariable
frequency_bands: ObservationAncilliarySimulationVariable
link_ends_delays: ObservationAncilliarySimulationVariable
range_conversion_factor: ObservationAncilliarySimulationVariable
received_frequency_intermediate: ObservationIntermediateSimulationVariable
reception_reference_frequency_band: ObservationAncilliarySimulationVariable
sequential_range_lowest_ranging_component: ObservationAncilliarySimulationVariable
transmitter_frequency_intermediate: ObservationIntermediateSimulationVariable