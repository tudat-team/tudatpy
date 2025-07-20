import numpy
import pybind11_stubgen.typing_ext
from ...astro import element_conversion
from ...astro import time_representation
from ...math import interpolators
import typing
from . import aerodynamic_coefficients
from . import atmosphere
from . import ephemeris
from . import gravity_field
from . import gravity_field_variation
from . import ground_station
from . import radiation_pressure
from . import rigid_body
from . import rotation_model
from . import shape
from . import shape_deformation
from . import vehicle_systems
__all__ = ['BodyListSettings', 'BodySettings', 'add_aerodynamic_coefficient_interface', 'add_empty_tabulated_ephemeris', 'add_engine_model', 'add_flight_conditions', 'add_gravity_field_model', 'add_ground_station', 'add_mass_properties_model', 'add_radiation_pressure_interface', 'add_radiation_pressure_target_model', 'add_rigid_body_properties', 'add_rotation_model', 'add_variable_direction_engine_model', 'aerodynamic_coefficients', 'atmosphere', 'convert_ground_station_state_between_itrf_frames', 'create_aerodynamic_coefficient_interface', 'create_body_ephemeris', 'create_ground_station_ephemeris', 'create_radiation_pressure_interface', 'create_simplified_system_of_bodies', 'create_system_of_bodies', 'create_tabulated_ephemeris_from_spice', 'ephemeris', 'get_default_body_settings', 'get_default_body_settings_time_limited', 'get_default_single_alternate_body_settings', 'get_default_single_alternate_body_settings_time_limited', 'get_default_single_body_settings', 'get_default_single_body_settings_time_limited', 'get_ground_station_list', 'get_safe_interpolation_interval', 'gravity_field', 'gravity_field_variation', 'ground_station', 'radiation_pressure', 'rigid_body', 'rotation_model', 'set_aerodynamic_guidance', 'set_aerodynamic_orientation_functions', 'set_constant_aerodynamic_orientation', 'shape', 'shape_deformation', 'vehicle_systems']

class BodyListSettings:
    """Class for defining settings for the creation of a system of bodies.
    
    Class for defining settings for the creation of a system of bodies. This object is typically created from default settings, and
    then adapted to the user's specific needs."""

    def __init__(self, frame_origin: str, frame_orientation: str) -> None:
        """
                 Class initialization method.
        
                 Class method to initialize an empty BodyListSettings object.
        
                 .. note::
        
                     When creating BodyListSettings from this method, the settings for each body will have to be added manually.
                     It is typically more convenient to use the :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings` function to create a BodyListSettings object with default settings for all bodies, and then modify the settings for specific bodies as needed.
        
        
                 Parameters
                 ----------
                 frame_origin : str
                     Definition of the global frame origin for the bodies.
                 frame_orientation : str
                     Definition of the global frame orientation for the bodies.
        """

    def add_empty_settings(self, body_name: str) -> None:
        """
                 This method adds empty settings to the :class:`BodyListSettings` instance.
        
                 Adds empty settings to the :class:`BodyListSettings` instance. This is typically used to add settings for custom bodies, for which no default settings are available.
                 See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/creation_celestial_body_settings.html>`_ for more information.
        
                 Parameters
                 ----------
                 body_name : str
                     Name of the body for which settings are added
        """

    def add_settings(self, settings_to_add: BodySettings, body_name: str) -> None:
        """
                 Add a single :class:`BodySettings` object to the :class:`BodyListSettings` instance.
        
                 .. warning::
        
                     This method is rarely called by the user, as :class:`BodySettings` objects cannot be created directly but only be extracted from a BodyListSettings instance.
                     Instead, users are recommended to use the :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings` to create settings for major celestial bodies, and the :func:`~tudatpy.dynamics.environment_setup.BodyListSettings.add_empty_settings` function to create settings for custom bodies.
                     See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/creation_celestial_body_settings.html>`_ for more information.
        
        
                 Parameters
                 ----------
                 settings_to_add : BodySettings
                     Settings to be added
                 body_name : str
                     Name of the body for which settings are added
        """

    def get(self, body_name: str) -> BodySettings:
        """
                 This function extracts a single BodySettings object.
        
        
                 Parameters
                 ----------
                 body_name : str
                     Name of the body for which settings are to be retrieved
        
        
                 Returns
                 -------
                 BodySettings
                     Settings for the requested body
        """

    @property
    def frame_orientation(self) -> str:
        """
                 **read-only**
        
                 Definition of the global frame orientation for the bodies
        
                 :type: str
        """

    @property
    def frame_origin(self) -> str:
        """
                 **read-only**
        
                 Definition of the global frame origin for the bodies
        
                 :type: str
        """

class BodySettings:
    """Class for defining settings for the creation of a single body.
    
    Class for defining settings for the creation of a single body, this object is typically stored inside a
    :class:`BodyListSettings` object."""

    @property
    def aerodynamic_coefficient_settings(self) -> aerodynamic_coefficients.AerodynamicCoefficientSettings:
        """
                 Object that defines the settings of the aerodynamic coefficient model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`aerodynamic_coefficients` module.
        
        
                 :type: AerodynamicCoefficientSettings
        """

    @aerodynamic_coefficient_settings.setter
    def aerodynamic_coefficient_settings(self, arg0: aerodynamic_coefficients.AerodynamicCoefficientSettings) -> None:
        ...

    @property
    def atmosphere_settings(self) -> atmosphere.AtmosphereSettings:
        """
                 Object that defines the settings of the atmosphere model that is to be created. Note that wind model settings
                 may be defined inside this object. A variable of this type is typically assigned by using a function from the
                 :ref:`atmosphere` module.
        
        
                 :type: AtmosphereSettings
        """

    @atmosphere_settings.setter
    def atmosphere_settings(self, arg0: atmosphere.AtmosphereSettings) -> None:
        ...

    @property
    def constant_mass(self) -> float:
        """
                 Mass that gets assigned to the vehicle. This mass does *not* automatically define a gravity field
                 model, but is instead used for the calculation of non-conservative forces only. When creating a body with a gravity field,
                 leave this entry empty. NOTE: this option is a shorthand for assigning a mass-only
                 :func:`~tudatpy.dynamics.environment_setup.rigid_body.constant_rigid_body_properties` to ``mass_property_settings``, and will be deprecated.
        
        
                 :type: float
        """

    @constant_mass.setter
    def constant_mass(self, arg0: float) -> None:
        ...

    @property
    def ephemeris_settings(self) -> ephemeris.EphemerisSettings:
        """
                 Object that defines the settings of the ephemeris model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`ephemeris` module.
        
        
                 :type: EphemerisSettings
        """

    @ephemeris_settings.setter
    def ephemeris_settings(self, arg0: ephemeris.EphemerisSettings) -> None:
        ...

    @property
    def gravity_field_settings(self) -> gravity_field.GravityFieldSettings:
        """
                 Object that defines the settings of the gravity field model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`gravity_field` module.
        
        
                 :type: GravityFieldSettings
        """

    @gravity_field_settings.setter
    def gravity_field_settings(self, arg0: gravity_field.GravityFieldSettings) -> None:
        ...

    @property
    def gravity_field_variation_settings(self) -> list[gravity_field_variation.GravityFieldVariationSettings]:
        """
                 List of objects that define the settings of time variations of the gravity field variation models that are to be created. Variables in this list are typically
                 assigned by using a function from the :ref:`gravity_field_variation` module.
        
        
                 :type: list[GravityFieldVariationSettings]
        """

    @gravity_field_variation_settings.setter
    def gravity_field_variation_settings(self, arg0: list[gravity_field_variation.GravityFieldVariationSettings]) -> None:
        ...

    @property
    def ground_station_settings(self) -> list[ground_station.GroundStationSettings]:
        """
        No documentation found.
        """

    @ground_station_settings.setter
    def ground_station_settings(self, arg0: list[ground_station.GroundStationSettings]) -> None:
        ...

    @property
    def radiation_pressure_settings(self) -> dict[str, radiation_pressure.RadiationPressureInterfaceSettings]:
        """
                 .. warning::
        
                     This interface is deprecated and will be removed in a future release. Use :attr:`~tudatpy.dynamics.environment_setup.BodySettings.radiation_source_settings` and :attr:`~tudatpy.dynamics.environment_setup.BodySettings.radiation_pressure_target_settings` instead.
        """

    @radiation_pressure_settings.setter
    def radiation_pressure_settings(self, arg0: dict[str, radiation_pressure.RadiationPressureInterfaceSettings]) -> None:
        ...

    @property
    def radiation_pressure_target_settings(self) -> radiation_pressure.RadiationPressureTargetModelSettings:
        """
                 Object that defines the settings of the radiation pressure target model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`radiation_pressure` module.
        
        
                 :type: RadiationPressureTargetModelSettings
        """

    @radiation_pressure_target_settings.setter
    def radiation_pressure_target_settings(self, arg0: radiation_pressure.RadiationPressureTargetModelSettings) -> None:
        ...

    @property
    def radiation_source_settings(self) -> radiation_pressure.RadiationSourceModelSettings:
        """
                 Object that defines the settings of the radiation source model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`radiation_pressure` module.
        
        
                 :type: RadiationSourceModelSettings
        """

    @radiation_source_settings.setter
    def radiation_source_settings(self, arg0: radiation_pressure.RadiationSourceModelSettings) -> None:
        ...

    @property
    def rigid_body_settings(self) -> rigid_body.RigidBodyPropertiesSettings:
        """
                 Object that defines the settings of the body rigid body (mass, center of mass, inertia) properties that are to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`rigid_body` module. Note that this setting does *not* define
                 the gravity field, but rather only the mass, center of mass and inertia tensor.
        
        
                 :type: RigidBodyPropertiesSettings
        """

    @rigid_body_settings.setter
    def rigid_body_settings(self, arg0: rigid_body.RigidBodyPropertiesSettings) -> None:
        ...

    @property
    def rotation_model_settings(self) -> rotation_model.RotationModelSettings:
        """
                 Object that defines the settings of the rotation model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`rotation_model` module.
        
        
                 :type: RotationModelSettings
        """

    @rotation_model_settings.setter
    def rotation_model_settings(self, arg0: rotation_model.RotationModelSettings) -> None:
        ...

    @property
    def shape_deformation_settings(self) -> list[shape_deformation.BodyDeformationSettings]:
        """
                 List of objects that define the settings of time variations of the exterior shape of natural bodies are to be created. Variables in this list are typically
                 assigned by using a function from the :ref:`shape_deformation` module.
        
        
                 :type: list[BodyDeformationSettings]
        """

    @shape_deformation_settings.setter
    def shape_deformation_settings(self, arg0: list[shape_deformation.BodyDeformationSettings]) -> None:
        ...

    @property
    def shape_settings(self) -> shape.BodyShapeSettings:
        """
                 Object that defines the settings of the shape model that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`shape` module.
        
        
                 :type: BodyShapeSettings
        """

    @shape_settings.setter
    def shape_settings(self, arg0: shape.BodyShapeSettings) -> None:
        ...

    @property
    def vehicle_shape_settings(self) -> vehicle_systems.FullPanelledBodySettings:
        """
                 Object that defines the settings of an exterior panelled vehicle shape that is to be created. A variable of this type is typically
                 assigned by using a function from the :ref:`vehicle_systems` module.
        
        
                 :type: FullPanelledBodySettings
        """

    @vehicle_shape_settings.setter
    def vehicle_shape_settings(self, arg0: vehicle_systems.FullPanelledBodySettings) -> None:
        ...

def add_aerodynamic_coefficient_interface(bodies: ..., body_name: str, coefficient_settings: aerodynamic_coefficients.AerodynamicCoefficientSettings) -> None:
    """Function that creates an aerodynamic coefficient interface from settings, and adds it to an existing body.
    
    This function can be used to add an aerodynamic coefficient interface to an existing body. It requires
    settings for the aerodynamic coefficients, created using one of the functions from the `~tudatpy.dynamics.environment_setup.aerodynamic_coefficients` module.
    This function creates the actual coefficient interface from these settings, and assigns it to the
    selected body. In addition to the identifier for the body to which it is assigned, this function
    requires the full :class:`~tudatpy.dynamics.environment.SystemOfBodies` as input, to facilitate
    inter-body dependencies in the coefficient interface
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_name : str
        Name of the body to which the aerodynamic coefficients are to be assigned
    coefficient_settings : AerodynamicCoefficientSettings
        Settings defining the coefficient interface that is to be created."""

def add_empty_tabulated_ephemeris(bodies: ..., body_name: str, ephemeris_origin: str='', is_part_of_multi_arc: bool=False) -> None:
    """No documentation found."""

def add_engine_model(body_name: str, engine_name: str, thrust_magnitude_settings: ..., bodies: ..., body_fixed_thrust_direction: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]=...) -> None:
    """Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.
    
    Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body. It creates and
    object of class :class:`~tudatpy.dynamics.environment.EngineModel`, and adds it to an existing body. Properties
    assigned to this engine model are:
    
    * The (constant) direction in body-fixed frame in which the engine is pointing (e.g. the body-fixed thrust direction when the engine is on)
    * Settings for computing the thrust magnitude (as a function of time and/or other parameters), using a suitable function from the :ref:`thrust` submodule
    
    
    Parameters
    ----------
    body_name : str
        Name of the body to which the engine is to be added.
    engine_name : str
        Name (e.g. unique identifier) of the engine that is to be added to the body
    thrust_magnitude_settings : ThrustMagnitudeSettings
        Settings for computing the thrust magnitude (and specific impulse) as a function of time
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_fixed_thrust_direction : numpy.ndarray[numpy.float64[3, 1]], default = [1,0,0]
        Unit vector along which the thrust from the engine will point in a body-fixed frame"""

def add_flight_conditions(bodies: ..., body_name: str, central_body_name: str) -> None:
    """Function that creates a flight conditions, and adds it to an existing body.
    
    This function can be used to add  a :class:`~tudatpy.dynamics.environment.FlightConditions` object to an existing body.
    Typically, the ``FlightConditions`` are created automatically when they are required (for the calculation of an
    aerodynamic acceleration, or the saving of certain dependent variables). However, in some cases it may be useful
    to manually trigger their creation, which is done through this function. If the ``central_body_name`` input
    denotes a body that is endowed with an :class:`~tudatpy.dynamics.environment.AtmosphereModel`, this function
    automatically creates an :class:`~tudatpy.dynamics.environment.AtmosphericFlightConditions` object (capable of
    calculating density, speed of sound, etc.), instead of the more basic :class:`~tudatpy.dynamics.environment.FlightConditions`
    (which is limited to properties such as altitude, latitude, etc.)
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_name : str
        Name of the body for which the flight conditions are to be created
    central_body_name : str
        Name of the central body w.r.t. which the flight conditions are to be created (typically, but not necessarily, the central body of propagation)/"""

def add_gravity_field_model(bodies: ..., body_name: str, gravity_field_settings: gravity_field.GravityFieldSettings, gravity_field_variation_settings: list[gravity_field_variation.GravityFieldVariationSettings]=[]) -> None:
    """No documentation found."""

@typing.overload
def add_ground_station(body: ..., ground_station_name: str, ground_station_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], position_type: element_conversion.PositionElementTypes=..., station_motion_settings: list[ground_station.GroundStationMotionSettings]=[]) -> None:
    ...

@typing.overload
def add_ground_station(body: ..., ground_station_settings: ground_station.GroundStationSettings) -> None:
    """No documentation found."""

def add_mass_properties_model(bodies: ..., body_name: str, mass_property_settings: rigid_body.RigidBodyPropertiesSettings) -> None:
    ...

def add_radiation_pressure_interface(bodies: ..., body_name: str, radiation_pressure_settings: radiation_pressure.RadiationPressureInterfaceSettings) -> None:
    ...

def add_radiation_pressure_target_model(bodies: ..., body_name: str, radiation_pressure_target_settings: radiation_pressure.RadiationPressureTargetModelSettings) -> None:
    """Function that creates a radiation pressure target model from settings, and adds it to an existing body.
    
    This function can be used to add a radiation pressure target model to an existing body. It requires
    settings for the radiation pressure target model, created using one of the functions from the :ref:`radiation_pressure` module.
    This function creates the actual target model from these settings, and assigns it to the
    selected body. In addition to the identifier for the body to which it is assigned, this function
    requires the full :class:`~tudatpy.dynamics.environment.SystemOfBodies` as input, to facilitate
    inter-body dependencies in the radiation pressure interface.
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_name : str
        Name of the body to which the radiation pressure interface is to be assigned
    radiation_pressure_target_settings : RadiationPressureTargetModelSettings
       Settings defining the radiation pressure target model that is to be created."""

def add_rigid_body_properties(bodies: ..., body_name: str, rigid_body_property_settings: rigid_body.RigidBodyPropertiesSettings) -> None:
    """Function that creates a rigid body property model, and adds it to an existing body.
    
    This function can be used to add a :class:`~tudatpy.dynamics.environment.RigidBodyProperties` object to an existing body.
    Typically, the ``RigidBodyProperties`` are created along with the :class:`~tudatpy.dynamics.environment.Body` itself. However, in some cases it may be useful
    to create body mass properties after the Body objects have been created. This function requires
    settings for the rigid body properties, created using one of the functions from the :ref:`rigid_body` module.
    This function creates the actual rigid body properties from these settings, and assigns it to the
    selected body.
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_name : str
        Name of the body to which the model is to be assigned
    rigid_body_property_settings : RigidBodyPropertiesSettings
        Settings defining the rigid body properties model that is to be created."""

def add_rotation_model(bodies: ..., body_name: str, rotation_model_settings: rotation_model.RotationModelSettings) -> None:
    """Function that creates a rotation model, and adds it to an existing body.
    
    This function can be used to add  a :class:`~tudatpy.dynamics.environment.RotationalEphemeris` object to an existing body.
    Typically, the ``RotationalEphemeris`` is created along with the :class:`~tudatpy.dynamics.environment.Body` itself. However, in some cases it may be useful
    to create a rotation model after the Body objects have been created. This function requires
    settings for the rotation model, created using one of the functions from the :ref:`rotation_model` module.
    This function creates the actual coefficient interface from these settings, and assigns it to the
    selected body. In addition to the identifier for the body to which it is assigned, this function
    requires the full :class:`~tudatpy.dynamics.environment.SystemOfBodies` as input, to facilitate
    inter-body dependencies in the radiation model
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_name : str
        Name of the body to which the rotation model is to be assigned
    rotation_model_settings
        Settings defining the rotation model that is to be created."""

def add_variable_direction_engine_model(body_name: str, engine_name: str, thrust_magnitude_settings: ..., bodies: ..., body_fixed_thrust_direction_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]) -> None:
    """Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.
    
    Same as :func:`add_engine_model`, but with a time-variable body-fixed thrust direction
    
    
    Parameters
    ----------
    body_name : str
        Name of the body to which the engine is to be added.
    engine_name : str
        Name (e.g. unique identifier) of the engine that is to be added to the body
    thrust_magnitude_settings : ThrustMagnitudeSettings
        Settings for computing the thrust magnitude (and specific impulse) as a function of time
    bodies : SystemOfBodies
        Object defining the physical environment, with all properties of artificial and natural bodies.
    body_fixed_thrust_direction_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
        Function returning a unit vector, as a function of time, along which the thrust from the engine will point in a body-fixed frame"""

def convert_ground_station_state_between_itrf_frames(ground_station_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], epoch: float, base_frame: str, target_frame: str) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """No documentation found."""

@typing.overload
def create_aerodynamic_coefficient_interface(coefficient_settings: aerodynamic_coefficients.AerodynamicCoefficientSettings, body: str) -> ...:
    ...

@typing.overload
def create_aerodynamic_coefficient_interface(coefficient_settings: aerodynamic_coefficients.AerodynamicCoefficientSettings, body: str, bodies: ...) -> ...:
    """No documentation found."""

def create_body_ephemeris(ephemeris_settings: ephemeris.EphemerisSettings, body_name: str) -> ...:
    """Function that creates an Ephemeris object.
    
    Function that creates an :class:`~tudatpy.dynamics.environment.Ephemeris` object, but does *not*
    associate it with any specific body (e.g., it does not go into the environment, but can be used independently of it)
    
    
    Parameters
    ----------
    ephemeris_settings : EphemerisSettings
        Object defining the ephemeris settings.
    body_name : str
        Name of body for which the ephemeris is created. Note that this input is only relevant for some ephemeris settings (for instance, a spice ephemeris setting), and it does *not* imply that the ephemeris object is associated with a Body object of this name.
    Returns
    -------
    :class:`~tudatpy.dynamics.environment.Ephemeris`
        Ephemeris object, created according to the provided settings"""

def create_ground_station_ephemeris(arg0: ..., arg1: str, arg2: ...) -> ...:
    """station_name"""

def create_radiation_pressure_interface(radiationPressureInterfaceSettings: radiation_pressure.RadiationPressureInterfaceSettings, body_name: str, body_dict: ...) -> ...:
    ...

def create_simplified_system_of_bodies(initial_time: float=0) -> ...:
    """Function that creates a simplified System of bodies.
    
    Function that creates a simplified system of bodies. The following bodies are created in this system: the Sun, all planets of the Solar system, and Pluto.
    All bodies in this system use Gtop ephemerides and point mass gravity. The Earth is setup with a spherical shape model and a simple rotation model.
    The reference frame used to setup this simplified system of bodies has its origin at the SSB, and has an ECLIPJ2000 orientation.
    
    
    Parameters
    ----------
    initial_time : float, optional, default=0
        Initial system time in seconds since J2000.
    Returns
    -------
    :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object containing the objects for bodies and environment models constituting the physical environment"""

def create_system_of_bodies(body_settings: BodyListSettings) -> ...:
    """Function that creates a System of bodies from associated settings.
    
    Function that creates a System of bodies from associated settings. This function creates the separate :class:`~tudatpy.dynamics.environment.Body`
    objects and stores them in a :class:`~tudatpy.dynamics.environment.SystemOfBodies` object. This object represents the full
    physical environment in the simulation.
    
    
    Parameters
    ----------
    body_settings : BodyListSettings
        Object defining the physical environment, with all properties of artificial and natural bodies.
    Returns
    -------
    :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object containing the objects for bodies and environment models constituting the physical environment"""

def create_tabulated_ephemeris_from_spice(body: str, initial_time: time_representation.Time, end_time: time_representation.Time, time_step: time_representation.Time, observer_name: str, reference_frame_name: str, interpolator_settings: interpolators.InterpolatorSettings=...) -> ...:
    ...

def get_default_body_settings(bodies: list[str], base_frame_origin: str='SSB', base_frame_orientation: str='ECLIPJ2000') -> BodyListSettings:
    """Function that retrieves the default settings for the given set of input bodies.
    
    Function that retrieves the default settings for the given set of input bodies. Default settings are described in
    detail `in the user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models.html>`_ .
    Note that if a body is provided as input for which default settings do not exist, an exception is thrown. In addition
    to settings for each separate body, this function returns an object that defines the global frame origin and orientation.
    
    .. note::
    
        Before using this function, make sure to have the appropriate set of SPICE kernels loaded.
        Typically, this is done through the :func:`~tudatpy.interface.spice.load_standard_kernels` function.
    
    
    Parameters
    ----------
    bodies : list[str]
        List of name of bodies for which default settings are to be retrieved and created.
    base_frame_origin : str, default = 'SSB'
        Base frame origin of the set of bodies that is to be created. It defaults to the solar system barycenter (SSB), but it can by any of the bodies in `bodies` (provided it has an ephemeris defined).
    base_frame_orientation : str, default = 'ECLIPJ2000'
        Base frame orientation of the set of bodies that is to be created. It can be either ECLIPJ2000 (default) or J2000.
    Returns
    -------
    BodyListSettings
        Object containing the settings for the SystemOfBodies that are to be created"""

def get_default_body_settings_time_limited(bodies: list[str], initial_time: float, final_time: float, base_frame_origin: str='SSB', base_frame_orientation: str='ECLIPJ2000', time_step: float=300.0) -> BodyListSettings:
    """Function that retrieves the default settings for the given set of input bodies, with a limited valid time interval.
    
    Same as :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings`, but with body settings valid over a limited time interval. This makes the
    the extraction of states from ephemerides more computationally efficient, at the expense of more RAM usage, and a
    constrained time interval over which the ephemerides are valid. See `this page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models/default_bodies_limited_time_range.html>`_ for more details.
    
    
    Parameters
    ----------
    bodies : list[str]
        List of name of bodies for which default settings are to be retrieved and created.
    initial_time : float
        Start time from which the environment settings should be created.
    final_time : float
        End time up to which the environment settings should be created.
    base_frame_origin : str
        Base frame origin of the set of bodies that is to be created.
    base_frame_orientation : str
        Base frame orientation of the set of bodies that is to be created.
    time_step : float, default = 300.0
        Time step to be used for the tabulated ephemeris.
    Returns
    -------
    BodyListSettings
        Object containing the settings for the SystemOfBodies that are to be created"""

def get_default_single_alternate_body_settings(body_name: str, source_body_name: str, base_frame_orientation: str='ECLIPJ2000') -> BodySettings:
    """Function that retrieves the default settings for a single body, and assigns them to another body.
    
    As :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings`, but for retrieving default settings of only a single body,
    where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
    For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
    of Mars
    
    
    Parameters
    ----------
    body_name : str
        Name of body for which default settings are to be created.
    source_body_name : str
        Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
    base_frame_orientation : str, default = 'ECLIPJ2000'
        Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
    Returns
    -------
    BodySettings
        Object containing the settings for the body that is to be created"""

def get_default_single_alternate_body_settings_time_limited(body_name: str, source_body_name: str, initial_time: float, final_time: float, base_frame_orientation: str='ECLIPJ2000', time_step: float=300.0) -> BodySettings:
    """Function that retrieves the default settings for a single body, with a limited valid time interval.
    
    As :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body,
    where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
    For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
    of Mars
    
    
    Parameters
    ----------
    body_name : str
        Name of body for which default settings are to be retrieved.
    source_body_name : str
        Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
    initial_time : float
        Start time from which the environment settings should be created.
    final_time : float
        End time up to which the environment settings should be created.
    base_frame_orientation : str, default = 'ECLIPJ2000'
        Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
    time_step : float, default = 300.0
        Time step to be used for the tabulated ephemeris.
    Returns
    -------
    BodySettings
        Object containing the settings for the body that is to be created"""

def get_default_single_body_settings(body_name: str, base_frame_orientation: str='ECLIPJ2000') -> BodySettings:
    """Function that retrieves the default settings for a single body.
    
    As :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings`, but for retrieving default settings of only a single body
    
    
    Parameters
    ----------
    body_name : str
        Name of body for which default settings are to be retrieved and created.
    base_frame_orientation : str, default = 'ECLIPJ2000'
        Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
    Returns
    -------
    BodySettings
        Object containing the settings for the body that is to be created"""

def get_default_single_body_settings_time_limited(body_name: str, initial_time: float, final_time: float, base_frame_orientation: str='ECLIPJ2000', time_step: float=300.0) -> BodySettings:
    """Function that retrieves the default settings for a single body, with a limited valid time interval.
    
    As :func:`~tudatpy.dynamics.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body
    
    
    Parameters
    ----------
    body_name : str
        Name of body for which default settings are to be retrieved and created.
    initial_time : float
        Start time from which the environment settings should be created.
    final_time : float
        End time up to which the environment settings should be created.
    base_frame_orientation : str, default = 'ECLIPJ2000'
        Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
    time_step : float, default = 300.0
        Time step to be used for the tabulated ephemeris.
    Returns
    -------
    BodySettings
        Object containing the settings for the body that is to be created"""

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