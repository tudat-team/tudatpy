import numpy
import tudatpy.data.expose_data
import tudatpy.numerical_simulation.estimation.expose_estimation
import typing
__all__ = ['DopplerProperTimeRateSettings', 'FrequencyBands', 'LightTimeConvergenceCriteria', 'LightTimeCorrectionSettings', 'LightTimeFailureHandling', 'LinkDefinition', 'LinkEndId', 'LinkEndType', 'NWayRangeObservationSettings', 'ObservableType', 'ObservationAncilliarySimulationSettings', 'ObservationAncilliarySimulationVariable', 'ObservationBiasSettings', 'ObservationDependentVariableSettings', 'ObservationSettings', 'ObservationSimulationSettings', 'ObservationViabilitySettings', 'ObservationViabilityType', 'OneWayDopplerObservationSettings', 'ProcessedOdfFileContents', 'TabulatedObservationSimulationSettings', 'TroposphericMappingModel', 'WaterVaporPartialPressureModel', 'absolute_bias', 'accept_without_warning', 'add_ancilliary_settings_to_observable', 'add_ancilliary_settings_to_observable_for_link_ends', 'add_dependent_variables_to_all', 'add_dependent_variables_to_observable', 'add_dependent_variables_to_observable_for_link_ends', 'add_gaussian_noise_to_all', 'add_gaussian_noise_to_observable', 'add_gaussian_noise_to_observable_for_link_ends', 'add_noise_function_to_all', 'add_noise_function_to_observable', 'add_noise_function_to_observable_for_link_ends', 'add_viability_check_to_all', 'add_viability_check_to_observable', 'add_viability_check_to_observable_for_link_ends', 'angular_position', 'angular_position_type', 'arc_wise_time_drift_bias', 'arcwise_absolute_bias', 'arcwise_relative_bias', 'bean_and_dutton', 'body_avoidance_angle', 'body_avoidance_viability', 'body_avoidance_viability_list', 'body_occultation', 'body_occultation_viability', 'body_occultation_viability_list', 'body_origin_link_end_id', 'body_reference_point_link_end_id', 'cartesian_position', 'cartesian_velocity', 'cassini_turnaround_ratios', 'change_simulation_settings_observable_types', 'combined_bias', 'continuous_arc_simulation_settings', 'continuous_arc_simulation_settings_list', 'create_odf_observed_observation_collection', 'create_tracking_txtfile_observation_collection', 'doppler_ancilliary_settings', 'doppler_integration_time', 'doppler_reference_frequency', 'dsn_default_turnaround_ratios', 'dsn_n_way_averaged_doppler', 'dsn_n_way_doppler_ancilliary_settings', 'dsn_n_way_doppler_averaged', 'dsn_n_way_doppler_averaged_from_one_way_links', 'dsn_one_way_averaged_doppler', 'dsn_tabulated_ionospheric_light_time_correction', 'dsn_tabulated_tropospheric_light_time_correction', 'elevation_angle_viability', 'elevation_angle_viability_list', 'euler_angle_313_observable_type', 'euler_angles_313', 'first_order_relativistic_light_time_correction', 'frequency_bands', 'get_default_reference_link_end', 'inverse_power_series_solar_corona_light_time_correction', 'jakowski_ionospheric_light_time_correction', 'light_time_convergence_settings', 'link_definition', 'link_ends_delays', 'minimum_elevation_angle', 'n_way_averaged_doppler_type', 'n_way_doppler_ancilliary_settings', 'n_way_doppler_averaged', 'n_way_doppler_averaged_from_one_way_links', 'n_way_range', 'n_way_range_ancilliary_settings', 'n_way_range_from_one_way_links', 'n_way_range_type', 'niell', 'observation_settings_from_collection', 'observed_body', 'one_way_averaged_doppler_type', 'one_way_closed_loop_doppler', 'one_way_doppler_averaged', 'one_way_doppler_instantaneous', 'one_way_downlink_link_ends', 'one_way_instantaneous_doppler_type', 'one_way_open_loop_doppler', 'one_way_range', 'one_way_range_type', 'one_way_uplink_link_ends', 'position_observable_type', 'print_warning_and_accept', 'process_odf_data_multiple_files', 'process_odf_data_single_file', 'receiver', 'reception_reference_frequency_band', 'reflector1', 'reflector2', 'reflector3', 'reflector4', 'relative_angular_position', 'relative_angular_position_type', 'relative_bias', 'relative_cartesian_position', 'retransmitter', 'saastamoinen_tropospheric_light_time_correction', 'set_odf_information_in_bodies', 'simplified_chao', 'tabulated', 'tabulated_simulation_settings', 'tabulated_simulation_settings_list', 'throw_exception', 'time_drift_bias', 'transmitter', 'two_doppler_instantaneous', 'two_way_doppler_ancilliary_settings', 'two_way_doppler_averaged', 'two_way_doppler_averaged_from_one_way_links', 'two_way_doppler_instantaneous_from_one_way_links', 'two_way_instantaneous_doppler_type', 'two_way_open_loop_doppler', 'two_way_open_loop_doppler_from_one_way_links', 'two_way_range', 'two_way_range_ancilliary_settings', 'two_way_range_from_one_way_links', 'unidentified_link_end', 'velocity_observable_type']

class DopplerProperTimeRateSettings:
    """Base class to defining proper time rate settings.
	
	Functional (base) class for settings of proper time rate (at a single link end) for instantaneous Doppler observation model settings.
	Specific proper time rate settings must be defined using an object derived from this class.
	The derived classes are made accessible via dedicated factory functions.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class FrequencyBands:
    """Members:
	
	s_band
	
	x_band
	
	ka_band
	
	ku_band
	"""
    __members__: typing.ClassVar[dict[str, FrequencyBands]]
    ka_band: typing.ClassVar[FrequencyBands]
    ku_band: typing.ClassVar[FrequencyBands]
    s_band: typing.ClassVar[FrequencyBands]
    x_band: typing.ClassVar[FrequencyBands]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class LightTimeConvergenceCriteria:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class LightTimeCorrectionSettings:
    """Base class to defining light time correction settings.
	
	Functional (base) class for settings of light time corrections.
	This class is not used for calculations of corrections, but is used for the purpose of defining the light time correction properties.
	Specific light time correction settings must be defined using an object derived from this class.
	The derived classes are made accessible via dedicated factory functions, such as e.g.
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction`
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class LightTimeFailureHandling:
    """Enumeration of behaviour when failing to converge light-time with required settings.
	
	
	:member accept_without_warning:
	:member print_warning_and_accept:
	:member throw_exception:
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, LightTimeFailureHandling]]
    accept_without_warning: typing.ClassVar[LightTimeFailureHandling]
    print_warning_and_accept: typing.ClassVar[LightTimeFailureHandling]
    throw_exception: typing.ClassVar[LightTimeFailureHandling]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class LinkDefinition:
    """Object storing the link ends involved in a given observation.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, link_ends: dict[LinkEndType, LinkEndId]) -> None:
        ...

class LinkEndId:
    """Object serving as identifier of a specific link end.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def body_name(self) -> str:
        ...

    @property
    def reference_point(self) -> str:
        ...

class LinkEndType:
    """Enumeration of available link end types.
	
	
	:member unidentified_link_end:
	:member transmitter:
	:member reflector1:
	:member retransmitter:
	:member reflector2:
	:member reflector3:
	:member reflector4:
	:member receiver:
	:member observed_body:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, LinkEndType]]
    observed_body: typing.ClassVar[LinkEndType]
    receiver: typing.ClassVar[LinkEndType]
    reflector1: typing.ClassVar[LinkEndType]
    reflector2: typing.ClassVar[LinkEndType]
    reflector3: typing.ClassVar[LinkEndType]
    reflector4: typing.ClassVar[LinkEndType]
    retransmitter: typing.ClassVar[LinkEndType]
    transmitter: typing.ClassVar[LinkEndType]
    unidentified_link_end: typing.ClassVar[LinkEndType]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class NWayRangeObservationSettings(ObservationSettings):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservableType:
    """Enumeration of available observable types.
	
	
	:member one_way_range_type:
	:member n_way_range_type:
	:member angular_position_type:
	:member relative_angular_position_type:
	:member position_observable_type:
	:member velocity_observable_type:
	:member one_way_instantaneous_doppler_type:
	:member one_way_averaged_doppler_type:
	:member two_way_instantaneous_doppler_type:
	:member n_way_averaged_doppler_type:
	:member euler_angle_313_observable_type:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, ObservableType]]
    angular_position_type: typing.ClassVar[ObservableType]
    dsn_n_way_averaged_doppler: typing.ClassVar[ObservableType]
    dsn_one_way_averaged_doppler: typing.ClassVar[ObservableType]
    euler_angle_313_observable_type: typing.ClassVar[ObservableType]
    n_way_averaged_doppler_type: typing.ClassVar[ObservableType]
    n_way_range_type: typing.ClassVar[ObservableType]
    one_way_averaged_doppler_type: typing.ClassVar[ObservableType]
    one_way_instantaneous_doppler_type: typing.ClassVar[ObservableType]
    one_way_range_type: typing.ClassVar[ObservableType]
    position_observable_type: typing.ClassVar[ObservableType]
    relative_angular_position_type: typing.ClassVar[ObservableType]
    two_way_instantaneous_doppler_type: typing.ClassVar[ObservableType]
    velocity_observable_type: typing.ClassVar[ObservableType]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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
    """Class for holding ancilliary settings for observation simulation.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def get_float_list_settings(self, setting_type: ObservationAncilliarySimulationVariable, throw_exception: bool=True) -> list[float]:
        """
        	:param setting_type:
        		Type of the setting for which the value is to be returned
        
        	:param throw_exception:
        		Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.
        
        	:return:
        		Value of the requested ancilliary variable (or empty list if it does not exist)
        """

    def get_float_settings(self, setting_type: ObservationAncilliarySimulationVariable, throw_exception: bool=True) -> float:
        """
        	:param setting_type:
        		Type of the setting for which the value is to be returned
        
        	:param throw_exception:
        		Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.
        
        	:return:
        		Value of the requested ancilliary variable (or NaN if it does not exist)
        """

class ObservationAncilliarySimulationVariable:
    """Members:
	
	link_ends_delays
	
	doppler_integration_time
	
	doppler_reference_frequency
	
	frequency_bands
	
	reception_reference_frequency_band
	"""
    __members__: typing.ClassVar[dict[str, ObservationAncilliarySimulationVariable]]
    doppler_integration_time: typing.ClassVar[ObservationAncilliarySimulationVariable]
    doppler_reference_frequency: typing.ClassVar[ObservationAncilliarySimulationVariable]
    frequency_bands: typing.ClassVar[ObservationAncilliarySimulationVariable]
    link_ends_delays: typing.ClassVar[ObservationAncilliarySimulationVariable]
    reception_reference_frequency_band: typing.ClassVar[ObservationAncilliarySimulationVariable]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class ObservationBiasSettings:
    """Base class to defining observation bias settings.
	
	Functional (base) class for settings of observation bias.
	Specific observation bias settings must be defined using an object derived from this class.
	The derived classes are made accessible via dedicated factory functions.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationDependentVariableSettings:
    """Base class for setting observation dependent variables.
	
	Functional (base) class for setting observation dependent variables as part of the observation output.
	Note: The associated functionality is not yet mature enough for the end user. Class is exposed for development purposes only.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationSettings:
    """Base class for defining observation settings.
	
	Functional (base) class for settings of observation models.
	Observation model settings define at least the type and geometry of a given observation.
	They can furthermore set observation biases and/or light-time corrections.
	Simple observation models settings that are fully characterised by these elements can be managed by this base class, which can be instantiated through dedicated factory functions, such as
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`, :func:`~tudatpy.numerical_simulation.estimation_setup.observation.cartesian_position`,
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.angular_position`, etc.
	Model settings for specific observation models that require additional information such as integration time, retransmission time, etc. must be defined using an object derived from this class.
	The derived classes are made accessible through further factory functions.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationSimulationSettings:
    """Base class for defining settings for simulating observations.
	
	Base class for defining settings for simulating observations.
	This simulation settings object defines observation times, noise and viability criteria, *etc.* at which observations are to be simulated.
	Therefore, one simulation settings object of this type can only refer to one combination of observable type and link geometry (LinkDefinition).
	The user does not interact with this class directly, but defines specific observation simulation settings using an object derived from this class (created through the associated factory function).
	"""
    observable_type: ObservableType
    viability_settings_list: list[ObservationViabilitySettings]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def link_ends(self) -> LinkDefinition:
        ...

    @property
    def noise_function(self) -> typing.Callable[[float], numpy.ndarray]:
        ...

    @noise_function.setter
    def noise_function(self, arg1: typing.Callable[[float], float]) -> None:
        ...

class ObservationViabilitySettings:
    """Enumeration of observation viability criterion types.
	
	
	:member minimum_elevation_angle:
	:member body_avoidance_angle:
	:member body_occultation:
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationViabilityType:
    """Enumeration of observation viability criterion types.
	
	
	:member minimum_elevation_angle:
	:member body_avoidance_angle:
	:member body_occultation:
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, ObservationViabilityType]]
    body_avoidance_angle: typing.ClassVar[ObservationViabilityType]
    body_occultation: typing.ClassVar[ObservationViabilityType]
    minimum_elevation_angle: typing.ClassVar[ObservationViabilityType]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class OneWayDopplerObservationSettings(ObservationSettings):
    """Class for defining the settings of one-way instantaneous Doppler observation models.
	
	:class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived class for one-way instantaneous Doppler observation model settings.
	Settings object can account for additional observation model aspects such as light time corrections and proper time rate settings.
	Instances of this class can be created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_instantaneous` factory function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ProcessedOdfFileContents:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def ground_station_names(self) -> list[str]:
        ...

    @property
    def ignored_ground_stations(self) -> list[str]:
        ...

    @property
    def ignored_odf_observable_types(self) -> list[int]:
        ...

    @property
    def processed_observable_types(self) -> list[ObservableType]:
        ...

    @property
    def raw_odf_data(self) -> list[tudatpy.data.expose_data.OdfRawFileContents]:
        ...

    @property
    def start_and_end_time(self) -> tuple[float, float]:
        ...

class TabulatedObservationSimulationSettings(ObservationSimulationSettings):
    """Class for defining settings for simulating observations at a predefined set of times.
	
	:class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived class for defining settings for simulating observations
	at a predefined set of times
	This type defines predefined time epochs at which applicable observations are to be simulated, stored in a rigid, "tabulated" form.
	Some observation times may be discarded due to the use of viability settings.
	Instances of this class are typicall created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
	and :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings_list` factory functions.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class TroposphericMappingModel:
    """Members:
	
	simplified_chao
	
	niell
	"""
    __members__: typing.ClassVar[dict[str, TroposphericMappingModel]]
    niell: typing.ClassVar[TroposphericMappingModel]
    simplified_chao: typing.ClassVar[TroposphericMappingModel]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class WaterVaporPartialPressureModel:
    """Members:
	
	tabulated
	
	bean_and_dutton
	"""
    __members__: typing.ClassVar[dict[str, WaterVaporPartialPressureModel]]
    bean_and_dutton: typing.ClassVar[WaterVaporPartialPressureModel]
    tabulated: typing.ClassVar[WaterVaporPartialPressureModel]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

def absolute_bias(bias_value: numpy.ndarray) -> ObservationBiasSettings:
    """Factory function for creating settings for an absolute observation bias.
	
	Factory function for creating settings for an absolute observation bias. When calculating the observable value, applying this setting
	will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`	ilde{h}` as follows:
	
	.. math::
			ilde{h}=h+K
	
	where :math:`K` is the `bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the addition is component-wise.
	
	
	:param bias_value:
			A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
			applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class, defining the settings for a constant, absolute observation bias.
	"""

def add_ancilliary_settings_to_observable(observation_simulation_settings_list: list[ObservationSimulationSettings], ancilliary_settings: ObservationAncilliarySimulationSettings, observable_type: ObservableType) -> None:
    ...

def add_ancilliary_settings_to_observable_for_link_ends(observation_simulation_settings_list: list[ObservationSimulationSettings], ancilliary_settings: ObservationAncilliarySimulationSettings, observable_type: ObservableType, link_ends: LinkDefinition) -> None:
    ...

def add_dependent_variables_to_all(observation_simulation_settings: list[ObservationSimulationSettings], dependent_variable_settings: list[ObservationDependentVariableSettings], bodies: ...) -> None:
    ...

def add_dependent_variables_to_observable(observation_simulation_settings: list[ObservationSimulationSettings], dependent_variable_settings: list[ObservationDependentVariableSettings], bodies: ..., observable_type: ObservableType) -> None:
    ...

def add_dependent_variables_to_observable_for_link_ends(observation_simulation_settings: list[ObservationSimulationSettings], dependent_variable_settings: list[ObservationDependentVariableSettings], bodies: ..., observable_type: ObservableType, link_ends: LinkDefinition) -> None:
    ...

def add_gaussian_noise_to_all(observation_simulation_settings_list: list[ObservationSimulationSettings], noise_amplitude: float) -> None:
    """Function for adding gaussian noise function to all existing observation simulation settings.
	
	Function for including simple time-independent and time-uncorrelated Gaussian noise function to the simulation settings of one or more observable(s).
	The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
	list.
	
	Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
	and thus the function does not return anything.
	
	
	:param observation_simulation_settings:
			Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param noise_amplitude:
			Standard deviation defining the un-biased Gaussian distribution for the noise.
	:return:
			The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.
	"""

def add_gaussian_noise_to_observable(observation_simulation_settings_list: list[ObservationSimulationSettings], noise_amplitude: float, observable_type: ObservableType) -> None:
    """Function for adding gaussian noise function to existing observation simulation settings of a given observable type.
	
	As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
	`observation_simulation_settings` list that matches the specified `observable_type`.
	
	
	:param observation_simulation_settings:
			Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param noise_amplitude:
			Standard deviation defining the un-biased Gaussian distribution for the noise.
	:param observable_type:
			Identifies the observable type in the observation simulation settings to which the noise is to be added.
	
	:return:
			The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.
	"""

def add_gaussian_noise_to_observable_for_link_ends(observation_simulation_settings_list: list[ObservationSimulationSettings], noise_amplitude: float, observable_type: ObservableType, link_definition: LinkDefinition) -> None:
    ...

def add_noise_function_to_all(observation_simulation_settings_list: list[ObservationSimulationSettings], noise_amplitude: typing.Callable[[float], numpy.ndarray]) -> None:
    ...

def add_noise_function_to_observable(observation_simulation_settings_list: list[ObservationSimulationSettings], noise_amplitude: typing.Callable[[float], numpy.ndarray], observable_type: ObservableType) -> None:
    ...

def add_noise_function_to_observable_for_link_ends(observation_simulation_settings_list: list[ObservationSimulationSettings], noise_amplitude: typing.Callable[[float], numpy.ndarray], observable_type: ObservableType, link_ends: LinkDefinition) -> None:
    ...

def add_viability_check_to_all(observation_simulation_settings_list: list[ObservationSimulationSettings], viability_settings: list[ObservationViabilitySettings]) -> None:
    """Function for including viability checks into existing observation simulation settings.
	
	Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
	The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
	list.
	Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
	and thus the function does not return anything.
	
	
	:param observation_simulation_settings:
			Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param viability_settings:
			List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.
	
	:return:
			The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.
	"""

def add_viability_check_to_observable(observation_simulation_settings_list: list[ObservationSimulationSettings], viability_settings: list[ObservationViabilitySettings], observable_type: ObservableType) -> None:
    ...

def add_viability_check_to_observable_for_link_ends(observation_simulation_settings_list: list[ObservationSimulationSettings], viability_settings: list[ObservationViabilitySettings], observable_type: ObservableType, link_ends: LinkDefinition) -> None:
    """Function for including viability checks into existing observation simulation settings.
	
	As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all`, except that the function only adds noise to entries of the
	`observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.
	
	
	:param observation_simulation_settings:
			Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param viability_settings:
			List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.
	
	:param observable_type:
			Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.
	
	:param link_definition:
			Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.
	
	:return:
			The :class
	"""

def angular_position(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for an angular position observable.
	
	Factory function for creating observation model settings of angular position type observables (as right ascension :math:`\x07lpha` and declination :math:`\\delta`),
	for a single link definition. The associated observation model creates an observable :math:`\\mathbf{h}_{_{	  ext{ang.pos.}}}` of type two as follows (in the unbiased case):
	
	.. math::
	   \\Delta\\mathbf{r}=\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T}(t_{T})\\
			an\x07lpha=
	
	   \\delta=
	
	   \\mathbf{h}_{_{	   ext{ang.pos.}}} = [\x07lpha;\\delta]
	
	The relative position vector :math:`\\Delta\\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
	The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
	environment.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the angular position observable.
	"""

@typing.overload
def arc_wise_time_drift_bias(bias_value: list[numpy.ndarray], arc_start_times: list[float], time_link_end: LinkEndType, ref_epochs: list[float]) -> ObservationBiasSettings:
    ...

@typing.overload
def arc_wise_time_drift_bias(bias_value_per_start_time: dict[float, numpy.ndarray], time_link_end: LinkEndType, ref_epochs: list[float]) -> ObservationBiasSettings:
    ...

@typing.overload
def arcwise_absolute_bias(arc_start_times: list[float], bias_values: list[numpy.ndarray], reference_link_end_type: LinkEndType) -> ObservationBiasSettings:
    """Factory function for creating settings for arc-wise absolute observation biases.
	
	Factory function for creating settings for arc-wise absolute observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.
	
	
	:param bias_values_per_start_time:
			Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
			The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)
	
	:param reference_link_end_type:
			Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
	"""

@typing.overload
def arcwise_absolute_bias(bias_values_per_start_time: dict[float, numpy.ndarray], reference_link_end_type: LinkEndType) -> ObservationBiasSettings:
    """Factory function for creating settings for arc-wise absolute observation biases.
	
	Factory function for creating settings for arc-wise absolute observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.
	
	
	:param bias_values_per_start_time:
			Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
			The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)
	
	:param reference_link_end_type:
			Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
	"""

@typing.overload
def arcwise_relative_bias(arc_start_times: list[float], bias_values: list[numpy.ndarray], reference_link_end_type: LinkEndType) -> ObservationBiasSettings:
    """Factory function for creating settings for arc-wise relative observation biases.
	
	Factory function for creating settings for arc-wise relative observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.
	
	
	:param bias_values_per_start_time:
			Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
			The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)
	
	:param reference_link_end_type:
			Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
	"""

@typing.overload
def arcwise_relative_bias(bias_values_per_start_time: dict[float, numpy.ndarray], reference_link_end_type: LinkEndType) -> ObservationBiasSettings:
    """Factory function for creating settings for arc-wise relative observation biases.
	
	Factory function for creating settings for arc-wise relative observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.
	
	
	:param bias_values_per_start_time:
			Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
			The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)
	
	:param reference_link_end_type:
			Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
	"""

def body_avoidance_viability(link_end_id: tuple[str, str], body_to_avoid: str, avoidance_angle: float) -> ObservationViabilitySettings:
    """Factory function for defining body avoidance observation viability settings.
	
	Factory function for defining body avoidance observation viability settings for single link ends.
	When simulating observations, this settings ensures that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
	The definition of 'too close' is computed as the angle between:
	
	* The line-of-sight vector from a link end to a given third body
	* The line-of-sight between two link ends
	
	This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
	a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.
	
	
	:param link_end_id:
			Link end (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
			To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
			For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
			link end passes too close to the specified body.
	
	:param body_to_avoid:
			Name of the body which the signal path should not pass 'too close' to.
	
	:param avoidance_angle:
			Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
			value must be in radians.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.
	"""

def body_avoidance_viability_list(link_end_ids: list[tuple[str, str]], body_to_avoid: str, avoidance_angle: float) -> list[ObservationViabilitySettings]:
    """Factory function for defining list of body avoidance viability settings.
	
	Factory function for defining body avoidance viability settings for multiple link ends.
	Each entry in the returned list contains the observation viability settings for one link end.
	When simulating observations, these settings ensure that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
	The definition of 'too close' is computed as the angle between:
	
	* The line-of-sight vector from a link end to a given third body
	* The line-of-sight between two link ends
	
	This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
	a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.
	
	
	:param link_end_ids:
			List of individual link ends (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
			To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
	
	:param body_to_avoid:
			Name of the body which the signal path should not pass 'too close' to.
	
	:param avoidance_angle:
			Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
			value must be in radians.
	
	:return:
			List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.
	"""

def body_occultation_viability(link_end_id: tuple[str, str], occulting_body: str) -> ObservationViabilitySettings:
    """Factory function for defining body occultation viability settings.
	
	Factory function for defining body occultation viability settings for single link ends.
	When simulating observations, this setting ensures that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
	The occultation is computed using the shape model of the specified body.
	
	
	:param link_end_id:
			Link end (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
			To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
	
	:param body_to_avoid:
			Name of the body which the signal path should not be occulted by.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.
	"""

def body_occultation_viability_list(link_end_id: tuple[str, str], occulting_body: str) -> ObservationViabilitySettings:
    """Factory function for defining body occultation viability settings.
	
	Factory function for defining body occultation viability settings for multiple link ends.
	Each entry in the returned list contains the observation viability settings for one link end.
	When simulating observations, these settings ensure that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
	The occultation is computed using the shape model of the specified body.
	
	
	:param link_end_ids:
			List of individual link ends (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
			To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
			For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
			link end is occulted by the specified body.
	
	:param body_to_avoid:
			Name of the body which the signal path should not be occulted by.
	
	:return:
			List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.
	"""

def body_origin_link_end_id(body_name: str) -> LinkEndId:
    """Function to create a link end identifier for the origin (typically center of mass) of a body.
	
	Function to create a link end identifier for the origin (typically center of mass) of a body.
	Using this option will simulate the origin of a body transmitter, receiving, etc. the observation.
	Although this is typically not physically realistic, it can be a useful approximation, in particular for simulation studies.
	
	
	:param body_name:
			Name of the body
	
	:return:
			A LinkEndId object representing the center of mass of a body
	"""

def body_reference_point_link_end_id(body_name: str, reference_point_id: str) -> LinkEndId:
    """Function to create a link end identifier for a reference point on a body.
	
	Function to create a link end identifier for a reference point on a body, where the reference point
	is typically the identifier of a ground stations
	
	
	:param body_name:
			Name of the body on which the reference point is located
	
	:param body_name:
			Name of the reference point on the body.
	
	:return:
			A LinkEndId object representing a reference point on a body
	"""

def cartesian_position(link_ends: LinkDefinition, bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for a Cartesian position observable.
	
	Factory function for creating observation model settings of Cartesian position type observables.
	Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
	This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
	The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires that the
			`observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian position observable.
	"""

def cartesian_velocity(link_ends: LinkDefinition, bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for a Cartesian velocity observable.
	
	Factory function for creating observation model settings of Cartesian position type observables.
	Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
	This observable provides the inertial (w.r.t. global frame origin) Cartesian velocity of the `observed_body` defined by the `link_ends` input.
	The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` velocity
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires that the
			`observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian velocity observable.
	"""

def cassini_turnaround_ratios(uplink_band: FrequencyBands, downlink_band: FrequencyBands) -> float:
    ...

def change_simulation_settings_observable_types(observation_simulation_settings: list[ObservationSimulationSettings], replacement_observable_types: dict[ObservableType, ObservableType]=...) -> None:
    ...

def combined_bias(bias_list: list[ObservationBiasSettings]) -> ObservationBiasSettings:
    """Factory function for creating settings for a combined observation bias.
			
				Factory function for creating settings for a combined observation bias, calculating by combining any number of bias types.
				Each contribution of the combined bias is computed from the unbiased observable, so when applying both a relative and absolute bias, we get:
			
				.. math::
						ilde{h}=h+K_{a}+hK_{r}
			
				And, crucially:
			
				.. math::
						ilde{h}
		eq (h+K_{a})(1+K_{r})
			
				where :math:`K_{r}` and :math:`K_{a}` is the relative and absolute bias, respectively.
			
			
				:param bias_list:
						A list containing the bias the bias settings that are to be applied to the observable.
			
				:return:
						Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`MultipleObservationBiasSettings` class, combining the settings for multiple observation biases.
			
		"""

def continuous_arc_simulation_settings(observable_type: ObservableType, link_ends: LinkDefinition, start_time: float, end_time: float, interval_between_observations: float, arc_limiting_constraints: ObservationViabilitySettings, minimum_arc_duration: float, maximum_arc_duration: float, minimum_time_between_arcs: float, reference_link_end_type: LinkEndType=..., additional_viability_settings: list[ObservationViabilitySettings]=[], noise_function: typing.Callable[[float], numpy.ndarray]=None) -> ObservationSimulationSettings:
    """Factory function for creating settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.
	
	Factory function for creating settings object for observation simulation. Unlike the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
	function, the resulting settings do not define the observation times explicitly. Instead, this settings object determines the observation times adaptively during the
	simulation of the observation, with the requirement that observations should be simulated over a set of contiguous arcs (if possible). The exact algorithm meets the following conditions:
	
	* Observations are only simulated within the time span of ``start_time`` and ``end_time``
	* A contiguous tracking arc has simulated observations separated by ``interval_between_observations``
	* Starting from ``start_time``, an observation is simulated each ``interval_between_observations``. Once an observation is unviable, as defined by
	  the ``arc_limiting_constraints`` input, it is checked whether the arc up until that point
	  is longer in duration than ``minimum_arc_duration``. If it is, the arc is added to the simulated observations. If not, the arc is discarded. In either case, a new arc is started once a
	  viable is observation is encountered
	* If the current arc reaching a duration greater than ``maximum_arc_duration``, the arc is added to the existing observations, and a new arc is started
	* If defined (e.g. if not NaN), the current observation time is incremented by ``minimum_time_between_arcs`` when an arc has been added to the observations.
	
	Nominally, this algorithm ensures that any arc of observations has a minimum and maximum duration. In addition, it ensures that (if desired) there is a minimum time interval
	between two tracking arcs. This behaviour can be modified by adding ``additional_viability_settings``, which are *not* used when computing the tracking arcs, but which are instead only used
	to reduce the set of simulated observations afterwards.
	
	
	:param observable_type:
			Observable type of which observations are to be simulated.
	:param link_ends:
			Link ends for which observations are to be simulated.
	:param start_time:
			First time at which an observation is to be simulated (and checked for viability).
	:param end_time:
			Maximum time at which an observation is to be simulated (and checked for viability).
	:param interval_between_observations:
			Cadence (in seconds) of subsequent observations in an arc
	:param arc_limiting_constraints:
			List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
			whether an arc should be terminated.
	
	:param minimum_arc_duration:
			Minimum permissible time for a tracking arc
	:param maximum_arc_duration:
			Maximum permissible time for a tracking arc
	:param minimum_time_between_arc:
			Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
	:param additional_viability_settings:
			Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
			These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.
	
	:param noise_function:
			Function providing the observation noise factors as a function of observation time.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.
	"""

def continuous_arc_simulation_settings_list(link_ends_per_observable: dict[ObservableType, list[LinkDefinition]], start_time: float, end_time: float, interval_between_observations: float, arc_limiting_constraints: ObservationViabilitySettings, minimum_arc_duration: float, maximum_arc_duration: float, minimum_time_between_arcs: float, reference_link_end_type: LinkEndType=..., additional_viability_settings: list[ObservationViabilitySettings]=[]) -> list[ObservationSimulationSettings]:
    ...

def create_odf_observed_observation_collection(processed_odf_file: ProcessedOdfFileContents, observable_types_to_process: list[ObservableType]=[], start_and_end_times_to_process: tuple[float, float]=...) -> tudatpy.numerical_simulation.estimation.expose_estimation.ObservationCollection:
    ...

def create_tracking_txtfile_observation_collection(raw_tracking_txtfile_contents: tudatpy.data.expose_data.TrackingTxtFileContents, spacecraft_name: str, observable_types_to_process: list[ObservableType]=[], earth_fixed_ground_station_positions: dict[str, numpy.ndarray]=..., ancillary_settings: ObservationAncilliarySimulationSettings=..., start_and_end_times_to_process: tuple[float, float]=...) -> tudatpy.numerical_simulation.estimation.expose_estimation.ObservationCollection:
    ...

def doppler_ancilliary_settings(integration_time: float=60.0) -> ObservationAncilliarySimulationSettings:
    ...

def dsn_default_turnaround_ratios(uplink_band: FrequencyBands, downlink_band: FrequencyBands) -> float:
    ...

def dsn_n_way_doppler_ancilliary_settings(frequency_bands: list[...], reference_frequency_band: ..., reference_frequency: float, integration_time: float=60.0, link_end_delays: list[float]=[]) -> ObservationAncilliarySimulationSettings:
    ...

def dsn_n_way_doppler_averaged(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    ...

def dsn_n_way_doppler_averaged_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    ...

def dsn_tabulated_ionospheric_light_time_correction(file_names: list[str], spacecraft_name_per_id: dict[int, str]={}, quasar_name_per_id: dict[int, str]={}, reference_frequency: float=2295000000.0, body_with_atmosphere_name: str='Earth') -> LightTimeCorrectionSettings:
    ...

def dsn_tabulated_tropospheric_light_time_correction(file_names: list[str], body_with_atmosphere_name: str='Earth', mapping_model: TroposphericMappingModel=...) -> LightTimeCorrectionSettings:
    ...

def elevation_angle_viability(link_end_id: tuple[str, str], elevation_angle: float) -> ObservationViabilitySettings:
    """Factory function for defining single elevation angle viability setting.
	
	Factory function for defining elevation angle viability settings for single link end.
	When simulating observations, this setting ensures that any applicable observations, for which the local elevation angle at link end is less than some limit value, will be omitted.
	
	
	:param link_end_id:
			Link end (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
			To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
	
	:param elevation_angle:
			Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
			value must be in radians.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` class, defining the settings for observation viability
	"""

def elevation_angle_viability_list(link_end_ids: list[tuple[str, str]], elevation_angle: float) -> list[ObservationViabilitySettings]:
    """Factory function for defining list of elevation angle viability settings.
	
	Factory function for defining elevation angle viability settings for multiple link ends.
	Each entry in the returned list contains the observation viability settings for one link end.
	When simulating observations, these settings ensure that any applicable observations, for which the local elevation angle at a link end is less than some limit value, will be omitted.
	
	
	:param link_end_ids:
			List of individual link ends (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
			To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
			For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
			link end violates the minimum elevation angle constraint.
	
	:param elevation_angle:
			Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
			value must be in radians.
	
	:return:
			List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.
	"""

def euler_angles_313(link_ends: LinkDefinition, bias_settings: ...=None) -> ObservationSettings:
    ...

def first_order_relativistic_light_time_correction(perturbing_bodies: list[str]) -> LightTimeCorrectionSettings:
    """Factory function for creating settings for first-order relativistic light-time corrections.
	
	Factory function for creating settings for first-order relativistic light-time corrections: the correction to
	the light time of a (set of) stationary point masses, computed up to c−2 according to general relativity as formulated by e.g. Moyer (2000).
	One ambiguity in the model is the time at which the states of the perturbing bodies are evaluated. We distinguish two cases:
	
	* In the case where the perturbing body contains a link end of the observation (for instance perturbation due to Earth gravity field,
	  with one of the link ends being an Earth-based station), the time at which the Earth’s state is evaluated equals the transmission time if Earth acts as transmitter, and reception time if
	  Earth acts as receiver.
	* In other cases, where the perturbing body is not involved in the link ends, its state is evaluated at the midpoint time between transmitter and receiver.
	
	
	:param perturbing_bodies:
			A list containing the names of the bodies due to which the light-time correction is to be taken into account.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` derived :class:`FirstOrderRelativisticLightTimeCorrectionSettings` class,
			defining the settings for the light-time corrections
	"""

def get_default_reference_link_end(observabl_type: ...) -> LinkEndType:
    ...

def inverse_power_series_solar_corona_light_time_correction(coefficients: list[float]=[7.8207e-06], positive_exponents: list[float]=[2.0], delay_coefficient: float=40.3, sun_body_name: str='Sun') -> LightTimeCorrectionSettings:
    ...

def jakowski_ionospheric_light_time_correction(ionosphere_height: float=400000.0, first_order_delay_coefficient: float=40.3, solar_activity_data: dict[float, tudatpy.data.expose_data.SolarActivityData]=..., geomagnetic_pole_latitude: float=1.4119713648634127, geomagnetic_pole_longitude: float=-1.2671090369478832, use_utc_for_local_time_computation: bool=False, body_with_atmosphere_name: str='Earth') -> LightTimeCorrectionSettings:
    ...

def light_time_convergence_settings(iterate_corrections: bool=False, maximum_number_of_iterations: int=50, absolute_tolerance: float=..., failure_handling: LightTimeFailureHandling=...) -> LightTimeConvergenceCriteria:
    """Factory function for creating settings for a one-way range observable.
	
	Factory function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{ ext{1-range}}}` as follows (in the unbiased case):
	
	.. math::
	   h_{_{		ext{1-range}}}(t_{R},t_{T})=|\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T}(t_{T})| + \\Delta s
	
	where :math:`\\mathbf{r}_{R}`, :math:`\\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
	of receiver and transmitter. The term :math:`\\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
	The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
	As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:
	
	.. math::
	   t_{R}-t_{T}=c\\left(|\\mathbf{r}_{R}(t_{R})-\\mathbf{r}(t_{R})| + \\Delta s
	
	
	 The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the one-way observable.
	"""

def link_definition(link_ends: dict[LinkEndType, LinkEndId]) -> LinkDefinition:
    """Function to create a link definition object.
	
	:param link_ends:
			Dictionary of link ends, with the key denoting the role in the observaton, and the associated value the identifier for the link end.
	:return:
			The ``LinkDefinition`` object storing the link ends of the observation
	"""

def n_way_doppler_ancilliary_settings(integration_time: float=60.0, link_end_delays: list[float]=[], frequency_bands: list[...]=[]) -> ObservationAncilliarySimulationSettings:
    """Factory function for creating ancilliary settings for n-way averaged Doppler observable.
	
	Factory function for creating ancilliary settings for a n-way averaged Doppler observable. Specifically, this
	function can be used to create settings for the integration time of the observable, and the  retransmission delays for each of the retransmitters.
	
	
	:param integration_time:
			Integration time that is to be used for the averaged Doppler observable
	:param retransmission_delays:
			Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
	"""

def n_way_doppler_averaged(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for an n-way averaged Doppler observable.
	
	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
	analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
	the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\\Delta t`.
	
	The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
			as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
	"""

def n_way_doppler_averaged_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for an n-way averaged Doppler observable.
	
	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition.
	The implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged`, with the difference
	that the constituent one-way range observables may have different settings.
	
	
	:param one_way_range_settings:
			List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way averaged range rate observable.
			The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
			``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter``(n-1) and ``receiver`` are defined by the
			``transmitter`` and ``receiver`` of the :math:`n`^{th} entry of this list.
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
	"""

def n_way_range(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for a n-way range observable.
	
	Factory function for creating observation model settings of n-way range type observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{ ext{N-range}}}` by combining together a series :math:`n` one-way range observations
	(see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`). By default, the reception time of the :math:`i^{th}` one-way range is set as the
	transmission time of the :math:`(i+1)^{th}` one-way range. A retransmission delay may be defined by ancilliary settings (see TODO) when creating observation
	simulation setings.
	
	For this factory function, the settings for each constituent one-way range (with the exception of the link end identifiers) are equal.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
			as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
			Note that only one bias setting is applied to the n-way observable.
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
	"""

def n_way_range_ancilliary_settings(link_end_delays: list[float]=[], frequency_bands: list[...]=[]) -> ObservationAncilliarySimulationSettings:
    """Factory function for creating ancilliary settings for n-way range observable.
	
	Factory function for creating ancilliary settings for a n-way range observable. Specifically, this
	function can be used to create settings for the retransmission delays of the observable, for each of the retransmitters.
	
	
	:param retransmission_delays:
			Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
	"""

def n_way_range_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for a n-way range observable.
	
	Factory function for creating observation model settings of n-way range type observables, for a single link definition. The
	implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with the difference
	that the constituent one-way ranges may have different settings.s
	
	
	:param one_way_range_settings:
			List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way range observable.
			The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
			``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter``(n-1) and ``receiver`` are defined by the
			``transmitter`` and ``receiver`` of the :math:`n`^{th} entry of this list.
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
			Note that only one bias setting is applied to the n-way observable.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
	"""

def observation_settings_from_collection(observed_observation_collection: tudatpy.numerical_simulation.estimation.expose_estimation.ObservationCollection) -> list[ObservationSimulationSettings]:
    ...

@typing.overload
def one_way_closed_loop_doppler(link_ends: LinkDefinition, light_time_correction_settings: list[LightTimeCorrectionSettings]=[], bias_settings: ObservationBiasSettings=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    ...

@typing.overload
def one_way_closed_loop_doppler(link_ends: LinkDefinition, light_time_correction_settings: list[LightTimeCorrectionSettings]=[], bias_settings: ObservationBiasSettings=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    ...

def one_way_doppler_averaged(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for a one-way averaged Doppler observable.
	
	Factory function for creating observation model settings for one-way averaged Doppler observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{ ext{1-\x08ar{Dopp}}}}` as follows (in the unbiased case):
	.. math::
	   h_{_{		ext{1-\x08ar{Dopp}}}}&=c\\int_{t-\\Delta t}^{t+\\Delta t}
	
								 &=
	   ext{1-range}}}(t_{R}=t+\\Delta t,t_{T})-h_{_{	ext{1-range}}}(t_{R}=t,t_{T})}{\\Delta t}
	
	where, in the latter formulation (which is the one that is implemented), the observable is referenced to the receiver time. This averaged Doppler observable
	is computed as the difference of two one-way range observables (see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`),
	with the reference time shifted by :math:`\\Delta t`. As such, it is sensitive to numerical errors for small :math:`\\Delta t`
	
	The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires that the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `OneWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
	"""

def one_way_doppler_instantaneous(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, transmitter_proper_time_rate_settings: DopplerProperTimeRateSettings=None, receiver_proper_time_rate_settings: DopplerProperTimeRateSettings=None, light_time_convergence_settings: LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    """Factory function for creating settings for a one-way instantaneous Doppler observable.
	
	Factory function for creating settings for a one-way instantaneous Doppler observable for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{ ext{1-Dopp.}}}` as follows (in the unbiased case):
	
	.. math::
	   h_{_{		ext{1-Dopp.}}}=c\\left(
	au_{T}}{dt_{T}}
	
		au_{R}}-1
	
	
	where :math:`t` and :math:`	 au` denote coordinate and proper time of the transmitter T and receiver R, respectively.
	The receiver and transmitter position and coordinate time are computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`.
	The detailed mathematical implementation are described on TODO.
	
	This observable represents the 'instantaneous (non-integrated)' Doppler observable, as obtained from open-loop observations.
	It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` should be used)
	
	The coordinate
	time derivative :math:`
	
	rates :math:`
	 au}{dt}` can be specified by the user through the ``transmitter_proper_time_rate_settings`` and ``receiver_proper_time_rate_settings``
	inputs. Whene these are left empty, the proper time rates are omitted (set to 1.0).
	
	The observable may be non-dimensionalized by the speed of light :math:`c`, which results in the observable being equal to thee received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires that the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param transmitter_proper_time_rate_settings:
			Settings for computing the transmitter proper time rate :math:`
	
		au}{dt}=1`)
	
	:param receiver_proper_time_rate_settings:
			Settings for computing the receiver proper time rate :math:`
	  au}{dt}`, default is none (:math:`
		au}{dt}=1`)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:param normalized_with_speed_of_light:
			Option to non-dimensionalize the observable with speed of light :math:`c`
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`OneWayDopplerObservationSettings` class defining the settings for the one-way open doppler observable observable.
	"""

def one_way_downlink_link_ends(transmitter: ..., receivers: list[...]) -> list[dict[LinkEndType, ...]]:
    """Function for defining one-way downlinks via LinkDefinition types.
	
	Function for defining single or multiple one-way downlinks.
	Multiple downlinks share the same transmitters, but may each have a different receiver.
	For each downlink, the returned list will contain an additional `LinkDefinition` type.
	
	
	:param transmitter:
			`LinkEndId` type (tuple of strings), where the first entrance identifies the body and the second entry the reference point of the single transmitter link end.
	
	:param receivers:
			List of `LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the receiver link end(s).
	
	:return:
			List of one or more `LinkDefinition` types, each defining the geometry for one one-way downlink.
			A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a `LinkEndId` type is mapped.
	"""

def one_way_open_loop_doppler(link_ends: LinkDefinition, light_time_correction_settings: list[LightTimeCorrectionSettings]=[], bias_settings: ObservationBiasSettings=None, transmitter_proper_time_rate_settings: DopplerProperTimeRateSettings=None, receiver_proper_time_rate_settings: DopplerProperTimeRateSettings=None, light_time_convergence_settings: LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    ...

def one_way_range(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for a one-way range observable.
	
	Factory function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{ ext{1-range}}}` as follows (in the unbiased case):
	
	.. math::
	   h_{_{		ext{1-range}}}(t_{R},t_{T})=|\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T}(t_{T})| + \\Delta s
	
	where :math:`\\mathbf{r}_{R}`, :math:`\\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
	of receiver and transmitter. The term :math:`\\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
	The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
	As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:
	
	.. math::
	   t_{R}-t_{T}=c\\left(|\\mathbf{r}_{R}(t_{R})-\\mathbf{r}(t_{R})| + \\Delta s
	
	
	 The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the one-way observable.
	"""

def one_way_uplink_link_ends(transmitters: list[...], receiver: ...) -> list[dict[LinkEndType, ...]]:
    """Function for defining one-way uplinks via LinkDefinition types.
	
	Function for defining single or multiple one-way uplinks.
	Multiple uplinks share the same receiver, but may each have a different transmitter.
	For each uplink, the returned list will contain an additional `LinkDefinition` type.
	
	
	:param transmitters:
			List of `LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the transmitter link end(s).
	
	:param receiver:
			`LinkEndId` type (tuple of strings), where the first entrance identifies the body and the second entry the reference point of the single receiver link end.
	
	:return:
			List of one or more `LinkDefinition` types, each defining the geometry for one one-way uplink.
			A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a `LinkEndId` type is mapped.
	"""

def process_odf_data_multiple_files(file_names: list[str], spacecraft_name: str, verbose: bool=True, earth_fixed_ground_station_positions: dict[str, numpy.ndarray]=...) -> ProcessedOdfFileContents:
    ...

def process_odf_data_single_file(file_name: str, spacecraft_name: str, verbose: bool=True, earth_fixed_ground_station_positions: dict[str, numpy.ndarray]=...) -> ProcessedOdfFileContents:
    ...

def relative_angular_position(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for an angular position observable.
	
	Factory function for creating observation model settings of angular position type observables (as right ascension :math:`\x07lpha` and declination :math:`\\delta`),
	for a single link definition. The associated observation model creates an observable :math:`\\mathbf{h}_{_{	  ext{ang.pos.}}}` of type two as follows (in the unbiased case):
	
	.. math::
	   \\Delta\\mathbf{r}=\\mathbf{r}_{R}(t_{R})-\\mathbf{r}_{T}(t_{T})\\
			an\x07lpha=
	
	   \\delta=
	
	   \\mathbf{h}_{_{	   ext{ang.pos.}}} = [\x07lpha;\\delta]
	
	The relative position vector :math:`\\Delta\\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
	The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
	environment.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the angular position observable.
	"""

def relative_bias(bias_value: numpy.ndarray) -> ObservationBiasSettings:
    """Factory function for creating settings for a relative observation bias.
	
	Factory function for creating settings for a relative observation bias. When calculating the observable value, applying this setting
	will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`	ilde{h}` as follows:
	
	.. math::
			ilde{h}=h(1+K)
	
	where :math:`K` is the`bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the multiplication is component-wise.
	
	
	:param bias_value:
			A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
			applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class,
			defining the settings for a constant, relative observation bias.
	"""

def relative_cartesian_position(link_ends: LinkDefinition, bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for a Cartesian position observable.
	
	Factory function for creating observation model settings of Cartesian position type observables.
	Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
	This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
	The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires that the
			`observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian position observable.
	"""

def saastamoinen_tropospheric_light_time_correction(body_with_atmosphere_name: str='Earth', mapping_model: TroposphericMappingModel=..., water_vapor_partial_pressure_model: WaterVaporPartialPressureModel=...) -> LightTimeCorrectionSettings:
    ...

def set_odf_information_in_bodies(processed_odf_file: ProcessedOdfFileContents, bodies: ..., body_with_ground_stations_name: str='Earth', turnaround_ratio_function: typing.Callable[[FrequencyBands, FrequencyBands], float]=...) -> None:
    ...

def tabulated_simulation_settings(observable_type: ObservableType, link_ends: LinkDefinition, simulation_times: list[float], reference_link_end_type: LinkEndType=..., viability_settings: list[ObservationViabilitySettings]=[], noise_function: typing.Callable[[float], numpy.ndarray]=None, ancilliary_settings: ObservationAncilliarySimulationSettings=None) -> ObservationSimulationSettings:
    """Factory function for creating settings object for observation simulation, using a predefined list of observation times.
	
	Factory function for creating single simulation settings object, using a predefined list of observation times.
	The list of resulting observations may be reduced compared to the ``simulation_times`` provided here, as
	only observations that meet the viability settings are retained during observation simulation (these may be
	provide directly here through the ``viability_settings`` input, or added later to the resulting settings object).
	
	
	:param observable_type:
			Observable type of which observations are to be simulated.
	:param link_ends:
			Link ends for which observations are to be simulated.
	:param simulation_times:
			List of times at which to perform the observation simulation.
	:param reference_link_end_type:
			Defines the link end (via the :class:`LinkEndType`) which is used as a reference time for the observation.
	:param viability_settings:
			Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
	
	:param noise_function:
			Function providing the observation noise factors as a function of observation time.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.
	"""

def tabulated_simulation_settings_list(link_ends_per_observable: dict[ObservableType, list[LinkDefinition]], simulation_times: list[float], reference_link_end_type: LinkEndType=..., viability_settings: list[ObservationViabilitySettings]=[]) -> list[ObservationSimulationSettings]:
    """Factory function for creating a list of settings object for observation simulation, using a predefined list of observation times.
	
	Factory function for creating multiple tabulated observation simulation settings objects in a list. This function is
	equivalent to calling the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings` repeatly, with the different
	observables and link definition provided here through `link_ends_per_observable`.
	During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.
	
	
	:param link_ends_per_observable:
			Link geometry per observable type of which observations are to be simulated.
	:param simulation_times:
			List of times at which to perform the observation simulation.
	:param reference_link_end_type:
			Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
			The single link end specified here will be considered as the reference link end for all simulation settings object created in the function call.
	
	:param viability_settings:
			Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
			The single settings list given here will be considered as potential viability settings for all simulation settings object created in the function call.
	
	:return:
			List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` objects.
	"""

def time_drift_bias(bias_value: numpy.ndarray, time_link_end: LinkEndType, ref_epoch: float) -> ObservationBiasSettings:
    ...

def two_doppler_instantaneous(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    ...

def two_way_doppler_ancilliary_settings(integration_time: float=60.0, retransmission_delay: float=0.0) -> ObservationAncilliarySimulationSettings:
    """Factory function for creating ancilliary settings for two-way averaged Doppler observable.
	
	Factory function for creating ancilliary settings for a two-way range observable. Specifically, this
	function can be used to create settings for the retransmission delay of the observable.  NOTE:
	this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_ancilliary_settings`
	with a single retransmission delay.
	
	
	:param integration_time:
			Integration time that is to be used for the averaged Doppler observable
	:param retransmission_delay:
			Retransmission delay that is to be applied to the simulation of the two-way observable
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
	"""

def two_way_doppler_averaged(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for an n-way averaged Doppler observable.
	
	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
	analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
	the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\\Delta t`.
	
	The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
			as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
	"""

def two_way_doppler_averaged_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for an n-way averaged Doppler observable.
	
	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
	analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
	the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\\Delta t`.
	
	The integration time :math:`\\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
			as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
	"""

def two_way_doppler_instantaneous_from_one_way_links(uplink_doppler_settings: OneWayDopplerObservationSettings, downlink_doppler_settings: OneWayDopplerObservationSettings, bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for a two-way instantaneous Doppler observable.
	
	
	Factory function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The
	implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_instantaneous`, with the difference
	that the constituent one-way ranges may have different settings.
	
	The observable may be non-dimensionalized by the speed of light :math:`c` (in the constituent one-way Doppler observable settings),
	which results in the observable being equal to the received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.
	
	
	:param uplink_doppler_settings:
			Settings for uplink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`
	
	:param downlink_doppler_settings:
			Settings for downlink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the full observation, default is none (unbiased observation). Note that,
			even if no bias is applied to the two-way observable, the constituent one-way observables may still be biased.
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`TwoWayDopplerObservationSettings` class defining the settings for the two-way open doppler observable.
	"""

def two_way_open_loop_doppler(link_ends: LinkDefinition, light_time_correction_settings: list[LightTimeCorrectionSettings]=[], bias_settings: ObservationBiasSettings=None, light_time_convergence_settings: LightTimeConvergenceCriteria=..., normalized_with_speed_of_light: bool=False) -> ObservationSettings:
    ...

def two_way_open_loop_doppler_from_one_way_links(uplink_doppler_settings: OneWayDopplerObservationSettings, downlink_doppler_settings: OneWayDopplerObservationSettings, bias_settings: ObservationBiasSettings=None) -> ObservationSettings:
    ...

def two_way_range(link_ends: LinkDefinition, light_time_correction_settings: list[...]=[], bias_settings: ...=None, light_time_convergence_settings: LightTimeConvergenceCriteria=...) -> ObservationSettings:
    """Factory function for creating settings for a two-way range observable.
	
	Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with :math:`n=2`. This function is provided
	for convenience.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the observation. This observable requires the
			`transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined
	
	:param light_time_correction_settings:
			List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
			in the signal being modelled as moving in a straight line with the speed of light
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
			Note that only one bias setting is applied to the n-way observable.
	
	:param light_time_convergence_settings:
			Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
	"""

def two_way_range_ancilliary_settings(retransmission_delay: float=0.0) -> ObservationAncilliarySimulationSettings:
    """Factory function for creating ancilliary settings for two-way range observable.
	
	Factory function for creating ancilliary settings for a two-way range observable. Specifically, this
	function can be used to create settings for the retransmission delay of the observable. NOTE:
	this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_ancilliary_settings`
	with a single retransmission delay.
	
	
	:param retransmission_delay:
			Retransmission delay that is to be applied to the simulation of the two-way observable
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
	"""

def two_way_range_from_one_way_links(one_way_range_settings: list[ObservationSettings], bias_settings: ...=None) -> ObservationSettings:
    """Factory function for creating settings for a two-way range observable.
	
	Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_from_one_way_links`, with :math:`n=2`. This function is provided
	for convenience.
	
	
	:param one_way_range_settings:
			List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
			The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
			``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
			``transmitter`` and ``receiver`` of the second entry of this list.
	
	:param bias_settings:
			Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
			Note that only one bias setting is applied to the n-way observable.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
	"""
accept_without_warning: LightTimeFailureHandling
angular_position_type: ObservableType
bean_and_dutton: WaterVaporPartialPressureModel
body_avoidance_angle: ObservationViabilityType
body_occultation: ObservationViabilityType
doppler_integration_time: ObservationAncilliarySimulationVariable
doppler_reference_frequency: ObservationAncilliarySimulationVariable
dsn_n_way_averaged_doppler: ObservableType
dsn_one_way_averaged_doppler: ObservableType
euler_angle_313_observable_type: ObservableType
frequency_bands: ObservationAncilliarySimulationVariable
link_ends_delays: ObservationAncilliarySimulationVariable
minimum_elevation_angle: ObservationViabilityType
n_way_averaged_doppler_type: ObservableType
n_way_range_type: ObservableType
niell: TroposphericMappingModel
observed_body: LinkEndType
one_way_averaged_doppler_type: ObservableType
one_way_instantaneous_doppler_type: ObservableType
one_way_range_type: ObservableType
position_observable_type: ObservableType
print_warning_and_accept: LightTimeFailureHandling
receiver: LinkEndType
reception_reference_frequency_band: ObservationAncilliarySimulationVariable
reflector1: LinkEndType
reflector2: LinkEndType
reflector3: LinkEndType
reflector4: LinkEndType
relative_angular_position_type: ObservableType
retransmitter: LinkEndType
simplified_chao: TroposphericMappingModel
tabulated: WaterVaporPartialPressureModel
throw_exception: LightTimeFailureHandling
transmitter: LinkEndType
two_way_instantaneous_doppler_type: ObservableType
unidentified_link_end: LinkEndType
velocity_observable_type: ObservableType