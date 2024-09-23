import numpy
import tudatpy.astro.element_conversion.expose_element_conversion
import tudatpy.math.interpolators.expose_interpolators
import tudatpy.numerical_simulation.environment_setup.gravity_field_variation.expose_gravity_field_variation
import tudatpy.numerical_simulation.environment_setup.ground_station.expose_ground_station
import typing
__all__ = ['BodyListSettings', 'BodySettings', 'add_aerodynamic_coefficient_interface', 'add_empty_tabulated_ephemeris', 'add_engine_model', 'add_flight_conditions', 'add_gravity_field_model', 'add_ground_station', 'add_mass_properties_model', 'add_radiation_pressure_interface', 'add_radiation_pressure_target_model', 'add_rigid_body_properties', 'add_rotation_model', 'add_variable_direction_engine_model', 'convert_ground_station_state_between_itrf_frames', 'create_aerodynamic_coefficient_interface', 'create_body_ephemeris', 'create_ground_station_ephemeris', 'create_radiation_pressure_interface', 'create_simplified_system_of_bodies', 'create_system_of_bodies', 'create_tabulated_ephemeris_from_spice', 'get_default_body_settings', 'get_default_body_settings_time_limited', 'get_default_single_alternate_body_settings', 'get_default_single_alternate_body_settings_time_limited', 'get_default_single_body_settings', 'get_default_single_body_settings_time_limited', 'get_ground_station_list', 'get_safe_interpolation_interval', 'set_aerodynamic_guidance', 'set_aerodynamic_orientation_functions', 'set_constant_aerodynamic_orientation']

class BodyListSettings:
    """Class for defining settings for the creation of a system of bodies.
	
	Class for defining settings for the creation of a system of bodies. This object is typically created from default settings, and
	then adapted to the user's specific needs.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, frame_origin: str, frame_orientation: str) -> None:
        ...

    def add_empty_settings(self, body_name: str) -> None:
        ...

    def add_settings(self, settings_to_add: BodySettings, body_name: str) -> None:
        ...

    def get(self, arg0: str) -> BodySettings:
        """
        This function extracts a single BodySettings object .
        
        	:param body_name:
        		Name of the body for which settings are to be retrieved
        """

    @property
    def frame_orientation(self) -> str:
        """
        Definition of the global frame orientation for the bodies
        	
        """

    @property
    def frame_origin(self) -> str:
        """
        Definition of the global frame origin for the bodies
        	
        """

class BodySettings:
    """Class for defining settings for the creation of a single body.
	
	Class for defining settings for the creation of a single body, this object is typically stored inside a
	:class:`BodyListSettings`, object.
	"""
    aerodynamic_coefficient_settings: ...
    gravity_field_variation_settings: list[tudatpy.numerical_simulation.environment_setup.gravity_field_variation.expose_gravity_field_variation.GravityFieldVariationSettings]
    ground_station_settings: list[tudatpy.numerical_simulation.environment_setup.ground_station.expose_ground_station.GroundStationSettings]
    radiation_pressure_settings: dict[str, ...]
    radiation_pressure_target_settings: ...
    radiation_source_settings: ...
    rigid_body_settings: ...
    shape_deformation_settings: list[...]
    vehicle_shape_settings: ...

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def atmosphere_settings(self) -> ...:
        """
        Object that defines the settings of the atmosphere model that is to be created. Note that wind model settings
        may be defined inside this object. A variable of this type is typically assigned by using a factory function from the
        :ref:`\\`\\`atmosphere\\`\\`` module.
        """

    @atmosphere_settings.setter
    def atmosphere_settings(self, arg0: ...) -> None:
        ...

    @property
    def constant_mass(self) -> float:
        """
        Mass that gets assigned to the vehicle. Note that this mass does *not* automatically define a gravity field
        model, but is instead used for the calculation of non-conservative forces only. When creating a body with a gravity field,
        leave this entry empty.
        """

    @constant_mass.setter
    def constant_mass(self, arg0: float) -> None:
        ...

    @property
    def ephemeris_settings(self) -> ...:
        """
        Object that defines the settings of the ephemeris model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\\`\\`ephemeris\\`\\`` module.
        """

    @ephemeris_settings.setter
    def ephemeris_settings(self, arg0: ...) -> None:
        ...

    @property
    def gravity_field_settings(self) -> ...:
        """
        Object that defines the settings of the gravity field model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\\`\\`gravity_field\\`\\`` module.
        """

    @gravity_field_settings.setter
    def gravity_field_settings(self, arg0: ...) -> None:
        ...

    @property
    def rotation_model_settings(self) -> ...:
        """
        Object that defines the settings of the rotation model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\\`\\`rotation_model\\`\\`` module.
        """

    @rotation_model_settings.setter
    def rotation_model_settings(self, arg0: ...) -> None:
        ...

    @property
    def shape_settings(self) -> ...:
        """
        Object that defines the settings of the shape model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\\`\\`shape\\`\\`` module.
        """

    @shape_settings.setter
    def shape_settings(self, arg0: ...) -> None:
        ...

def add_aerodynamic_coefficient_interface(bodies: ..., body_name: str, coefficient_settings: ...) -> None:
    """Function that creates an aerodynamic coefficient interface from settings, and adds it to an existing body.
	
	This function can be used to add an aerodynamic coefficient interface to an existing body. It requires
	settings for the aerodynamic coefficients, created using one of the factory functions from the `~tudatpy.numerical_simulation_environment_setup.aerodynamic_coefficient` module.
	This function creates the actual coefficient interface from these settings, and assigns it to the
	selected body. In addition to the identifier for the body to which it is assigned, this function
	requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
	inter-body dependencies in the coefficient interface
	
	
	:param bodies:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
			Name of the body to which the aerodynamic coefficients are to be assigned
	:param coefficient_settings:
			Settings defining the coefficient interface that is to be created.
	"""

def add_empty_tabulated_ephemeris(bodies: ..., body_name: str, ephemeris_origin: str='', is_part_of_multi_arc: bool=False) -> None:
    ...

def add_engine_model(body_name: str, engine_name: str, thrust_magnitude_settings: ..., bodies: ..., body_fixed_thrust_direction: numpy.ndarray=...) -> None:
    """Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.
	
	Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body. It creates and
	object of class :class:`~tudatpy.numerical_simulation.environment.EngineModel`, and adds it to an existing body. Properties
	assigned to this engine model are:
	* The (constant) direction in body-fixed frame in which the engine is pointing (e.g. the body-fixed thrust direction when the engine is on)
	* Settings for computing the thrust magnitude (as a function of time and/or other parameters), using a suitable function from the :ref:`\\`\\`thrust\\`\\`` submodule
	
	
	:param body_name:
			Name of the body to which the engine is to be added.
	:param engine_name:
			Name (e.g. unique identifier) of the engine that is to be added to the body
	:param thrust_magnitude_settings:
			Settings for computing the thrust magnitude (and specific impulse) as a function of time
	:param bodies:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
			Name of the body to which the rotation model is to be assigned
	:param body_fixed_thrust_direction:
			Unit vector along which the thrust from the engine will point in a body-fixed frame
	"""

def add_flight_conditions(bodies: ..., body_name: str, central_body_name: str) -> None:
    """Function that creates a flight conditions, and adds it to an existing body.
	
	This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.FlightConditions` object to an existing body.
	Typically, the ``FlightConditions`` are created automatically when they are required (for the calulcation of an
	aerodynamic acceleration, or the saving of certain dependent variables). However, in some cases it may be useful
	to manually trigger their creation, which is done through this function. If the ``central_body_name`` input
	denotes a body that is endowed with an :class:`~tudatpy.numerical_simulation.environment.AtmosphereModel`, this function
	automically creates an :class:`~tudatpy.numerical_simulation.environment.AtmosphericFlightConditions` object (capable of
	calculating density, speed of sound, etc.), instead of the more basic :class:`~tudatpy.numerical_simulation.environment.FlightConditions`
	(which is limited to properties such as altitude, latitude, etc.)
	
	
	:param bodies:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
			Name of the body for which the flight conditions are to be created
	:param central_body_name:
			Name of the cenral body w.r.t. which the flight conditions are to be created (typically, but not necesarilly, the central body of propagation)/
	"""

def add_gravity_field_model(bodies: ..., body_name: str, gravity_field_settings: ..., gravity_field_variation_settings: list[tudatpy.numerical_simulation.environment_setup.gravity_field_variation.expose_gravity_field_variation.GravityFieldVariationSettings]=[]) -> None:
    ...

@typing.overload
def add_ground_station(body: ..., ground_station_name: str, ground_station_position: numpy.ndarray, position_type: tudatpy.astro.element_conversion.expose_element_conversion.PositionElementTypes=..., station_motion_settings: list[tudatpy.numerical_simulation.environment_setup.ground_station.expose_ground_station.GroundStationMotionSettings]=[]) -> None:
    ...

@typing.overload
def add_ground_station(body: ..., ground_station_settings: tudatpy.numerical_simulation.environment_setup.ground_station.expose_ground_station.GroundStationSettings) -> None:
    ...

def add_mass_properties_model(bodies: ..., body_name: str, mass_property_settings: ...) -> None:
    ...

def add_radiation_pressure_interface(bodies: ..., body_name: str, radiation_pressure_settings: ...) -> None:
    ...

def add_radiation_pressure_target_model(bodies: ..., body_name: str, radiation_pressure_target_settings: ...) -> None:
    """Function that creates an radiation pressure interface from settings, and adds it to an existing body.
	
	This function can be used to add an radiation pressure interface to an existing body. It requires
	settings for the radiation pressure interface, created using one of the factory functions from the :ref:`\\`\\`radiation_pressure\\`\\`` module.
	This function creates the actual coefficient interface from these settings, and assigns it to the
	selected body. In addition to the identifier for the body to which it is assigned, this function
	requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
	inter-body dependencies in the radiation pressure interface
	
	
	:param bodies:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
			Name of the body to which the radiation pressure interface is to be assigned
	:param radiation_pressure_settings:
			Settings defining the radiation pressure interface that is to be created.
	"""

def add_rigid_body_properties(bodies: ..., body_name: str, rigid_body_property_settings: ...) -> None:
    ...

def add_rotation_model(bodies: ..., body_name: str, rotation_model_settings: ...) -> None:
    """Function that creates a rotation model, and adds it to an existing body.
	
	This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.RotationalEphemeris` object to an existing body.
	Typically, the ``RotationalEphemeris`` is created along with the `~tudatpy.numerical_simulation.environment.Body` itself However, in some cases it may be useful
	to create a rotation model after the Body objects have been created. This function requires
	settings for the rotation model, created using one of the factory functions from the :ref:`~tudatpy.numerical_simulation_environment_setup.rotation_model` module.
	This function creates the actual coefficient interface from these settings, and assigns it to the
	selected body. In addition to the identifier for the body to which it is assigned, this function
	requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
	inter-body dependencies in the radiation model
	
	
	:param bodies:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
			Name of the body to which the rotation model is to be assigned
	:param rotation_model_settings:
			Settings defining the rotation model that is to be created.
	"""

def add_variable_direction_engine_model(body_name: str, engine_name: str, thrust_magnitude_settings: ..., bodies: ..., body_fixed_thrust_direction_function: typing.Callable[[float], numpy.ndarray]) -> None:
    """Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.
	
	Same as :func:`add_engine_model`, but with a time-variable body-fixed thrust direction
	
	
	:param body_name:
			Name of the body to which the engine is to be added.
	:param engine_name:
			Name (e.g. unique identifier) of the engine that is to be added to the body
	:param thrust_magnitude_settings:
			Settings for computing the thrust magnitude (and specific impulse) as a function of time
	:param bodies:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
			Name of the body to which the rotation model is to be assigned
	:param body_fixed_thrust_direction_function:
			Function returning a unit vector, as a function of time, along which the thrust from the engine will point in a body-fixed frame
	"""

def convert_ground_station_state_between_itrf_frames(ground_station_state: numpy.ndarray, epoch: float, base_frame: str, target_frame: str) -> numpy.ndarray:
    ...

@typing.overload
def create_aerodynamic_coefficient_interface(coefficient_settings: ..., body: str) -> ...:
    ...

@typing.overload
def create_aerodynamic_coefficient_interface(coefficient_settings: ..., body: str, bodies: ...) -> ...:
    ...

def create_body_ephemeris(ephemeris_settings: ..., body_name: str) -> ...:
    ...

def create_ground_station_ephemeris(arg0: ..., arg1: str, arg2: ...) -> ...:
    """station_name
	"""

def create_radiation_pressure_interface(radiationPressureInterfaceSettings: ..., body_name: str, body_dict: ...) -> ...:
    ...

def create_simplified_system_of_bodies(initial_time: float=0) -> ...:
    """Function that creates a simplified System of bodies.
	
	Function that creates a simplified system of bodies. The following bodies are created in this system: the Sun, all planets of the Solar system, and Pluto.
	All bodies in this system use Gtop ephemerides and point mass gravity. The Earth is setup with a spherical shape model and a simple rotation model.
	The reference frame used to setup this simplified system of bodies has its origin at the SSB, and has an ECLIPJ2000 orientation.
	
	
	:param initial_time:
			Initial system time in seconds since J2000.
	:return:
			Object containing the objects for bodies and environment models constituting the physical environment
	"""

def create_system_of_bodies(body_settings: BodyListSettings) -> ...:
    """Function that creates a System of bodies from associated settings.
	
	Function that creates a System of bodies from associated settings. This function creates the separate :class:`~tudatpy.numerical_simulation.Body`
	objects and stores them in a :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. This object represents the full
	physical environment in the simulation.
	
	
	:param body_settings:
			Object defining the physical environment, with all properties of artificial and natural bodies.
	:return:
			Object containing the objects for bodies and environment models constituting the physical environment
	"""

def create_tabulated_ephemeris_from_spice(body: str, initial_time: float, end_time: float, time_step: float, observer_name: str, reference_frame_name: str, interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=...) -> ...:
    ...

def get_default_body_settings(bodies: list[str], base_frame_origin: str='SSB', base_frame_orientation: str='ECLIPJ2000') -> BodyListSettings:
    """Function that retrieves the default settings for the given set of input bodies.
	
	Function that retrieves the default settings for the given set of input bodies. Default settings are described in
	detail `here <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/create_bodies/default_settings.html>`_ .
	Note that if a body is provided as input for which default settings do not exist, an exception is thrown. In addition
	to settings for each separate body, this function returns an object that defines the global frame origin and orientation,
	
	
	:param bodies:
			List of name of bodies for which default settings are to be retrieved.
	:param base_frame_origin:
			Base frame origin of the set of bodies that is to be created. It defaults to the solar system barycenter (SSB), but it can by any of the bodies in `bodies_to_create` (provided it has an ephemeris defined).
	:param base_frame_orientation:
			Base frame orientation of the set of bodies that is to be created. It can be either ECLIPJ2000 (default) or J2000.
	:return:
			Object containing the settings for the SystemOfBodies that are to be created
	"""

def get_default_body_settings_time_limited(bodies: list[str], initial_time: float, final_time: float, base_frame_origin: str='SSB', base_frame_orientation: str='ECLIPJ2000', time_step: float=300.0) -> BodyListSettings:
    """Function that retrieves the default settings for the given set of input bodies, with a limited valid time interval.
	
	Same as :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but with body settings valid over a limited time interval. This makes the
	the extraction of states from ephemerides more computationally efficient, at the expense of more RAM usage, and a
	constrained time interval over which the ephemerides are valid. See `this page <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/valid_time_range.html>`_ for more details.
	
	
	:param bodies:
			List of name of bodies for which default settings are to be retrieved.
	:param initial_time:
			Start time from which the environment settings should be created.
	:param final_time:
			End time up to which the environment settings should be created.
	:param base_frame_origin:
			Base frame origin of the set of bodies that is to be created.
	:param base_frame_orientation:
			Base frame orientation of the set of bodies that is to be created.
	:param time_step:
			Time step to be used for the tabulated ephemeris.
	:return:
			Object containing the settings for the SystemOfBodies that are to be created
	"""

def get_default_single_alternate_body_settings(body_name: str, source_body_name: str, base_frame_orientation: str='ECLIPJ2000') -> BodySettings:
    ...

def get_default_single_alternate_body_settings_time_limited(body_name: str, source_body_name: str, initial_time: float, final_time: float, base_frame_orientation: str='ECLIPJ2000', time_step: float=300.0) -> BodySettings:
    ...

def get_default_single_body_settings(body_name: str, base_frame_orientation: str='ECLIPJ2000') -> BodySettings:
    ...

def get_default_single_body_settings_time_limited(body_name: str, initial_time: float, final_time: float, base_frame_orientation: str='ECLIPJ2000', time_step: float=300.0) -> BodySettings:
    ...

def get_ground_station_list(body: ...) -> list[tuple[str, str]]:
    ...

def get_safe_interpolation_interval(ephemeris_model: ...) -> tuple[float, float]:
    ...

def set_aerodynamic_guidance(aerodynamic_guidance: ..., body: ..., silence_warnings: bool=False) -> None:
    ...

def set_aerodynamic_orientation_functions(body: ..., angle_of_attack_function: typing.Callable[[], float]=None, sideslip_angle_function: typing.Callable[[], float]=None, bank_angle_function: typing.Callable[[], float]=None, update_function: typing.Callable[[float], None]=None) -> None:
    ...

def set_constant_aerodynamic_orientation(body: ..., angle_of_attack: float, sideslip_angle: float, bank_angle: float, silence_warnings: bool=False) -> None:
    ...