import numpy
import pybind11_stubgen.typing_ext
from ....dynamics import environment
from ....dynamics.environment_setup import aerodynamic_coefficients
from ....dynamics.environment_setup import gravity_field_variation
from ....dynamics.propagation_setup import acceleration
from ....dynamics.propagation_setup import torque
import typing
__all__ = ['PropagationDependentVariables', 'SingleAccelerationDependentVariableSaveSettings', 'SingleDependentVariableSaveSettings', 'VariableSettings', 'acceleration_partial_wrt_body_translational_state_type', 'actual_cross_section', 'aerodynamic_coefficients', 'aerodynamic_force_coefficients', 'aerodynamic_force_coefficients_control_surface_free', 'aerodynamic_force_coefficients_control_surface_increment', 'aerodynamic_force_coefficients_type', 'aerodynamic_moment_coefficients', 'aerodynamic_moment_coefficients_control_surface_free', 'aerodynamic_moment_coefficients_control_surface_increment', 'aerodynamic_moment_coefficients_type', 'airspeed', 'airspeed_type', 'altitude', 'altitude_type', 'angle_of_attack', 'apoapsis_altitude', 'apoapsis_altitude_type', 'bank_angle', 'body_fixed_airspeed_based_velocity_type', 'body_fixed_airspeed_velocity', 'body_fixed_groundspeed_based_velocity_type', 'body_fixed_groundspeed_velocity', 'body_fixed_relative_cartesian_position_type', 'body_fixed_relative_spherical_position_type', 'body_mass', 'center_of_mass', 'central_body_fixed_cartesian_position', 'central_body_fixed_spherical_position', 'control_surface_deflection', 'control_surface_deflection_type', 'cross_section_change', 'current_body_mass_type', 'custom', 'custom_dependent_variable', 'custom_type', 'density', 'dynamic_pressure', 'euler_angles_to_body_fixed_type', 'flight_path_angle', 'full_body_paneled_geometry', 'geodetic_latitude', 'geodetic_latitude_type', 'get_dependent_variable_id', 'get_dependent_variable_shape', 'get_dependent_variable_size', 'gravity_field_laplacian_of_potential', 'gravity_field_laplacian_of_potential_type', 'gravity_field_potential', 'gravity_field_potential_type', 'heading_angle', 'illuminated_panel_fraction', 'inertia_tensor', 'inertial_to_body_fixed_313_euler_angles', 'inertial_to_body_fixed_rotation_frame', 'intermediate_aerodynamic_rotation_matrix_type', 'intermediate_aerodynamic_rotation_matrix_variable', 'keplerian_state', 'keplerian_state_type', 'latitude', 'local_aerodynamic_g_load', 'local_density_type', 'local_dynamic_pressure_type', 'local_temperature_type', 'longitude', 'mach_number', 'mach_number_type', 'minimum_body_distance', 'minimum_visible_station_body_distances', 'modified_equinoctial_state', 'modified_equinoctial_state_type', 'per_target_panel_radiation_pressure_force', 'periapsis_altitude', 'periapsis_altitude_type', 'radiation_pressure', 'radiation_pressure_coefficient', 'radiation_pressure_coefficient_type', 'radiation_pressure_source_panel_geometry', 'radiation_pressure_source_panel_irradiance', 'radiation_pressure_type', 'received_irradiance', 'received_irradiance_shadow_function', 'relative_body_aerodynamic_orientation_angle_type', 'relative_distance', 'relative_distance_type', 'relative_position', 'relative_position_type', 'relative_speed', 'relative_speed_type', 'relative_velocity', 'relative_velocity_type', 'rotation_matrix_to_body_fixed_frame_type', 'rsw_to_inertial_frame_rotation_type', 'rsw_to_inertial_rotation_matrix', 'sideslip_angle', 'single_acceleration', 'single_acceleration_norm', 'single_acceleration_norm_type', 'single_acceleration_type', 'single_gravity_field_variation_acceleration', 'single_gravity_field_variation_acceleration_terms_type', 'single_gravity_field_variation_acceleration_type', 'single_per_term_gravity_field_variation_acceleration', 'single_torque', 'single_torque_norm', 'single_torque_norm_type', 'single_torque_type', 'spherical_harmonic_acceleration_norm_terms_type', 'spherical_harmonic_acceleration_terms_type', 'spherical_harmonic_terms_acceleration', 'spherical_harmonic_terms_acceleration_norm', 'stagnation_point_heat_flux_type', 'temperature', 'tnw_to_inertial_frame_rotation_type', 'tnw_to_inertial_rotation_matrix', 'total_acceleration', 'total_acceleration_norm', 'total_acceleration_norm_type', 'total_acceleration_type', 'total_aerodynamic_g_load_type', 'total_gravity_field_variation_acceleration', 'total_gravity_field_variation_acceleration_type', 'total_mass_rate', 'total_mass_rate_type', 'total_spherical_harmonic_cosine_coefficien_variations', 'total_spherical_harmonic_cosine_coefficien_variations_from_indices', 'total_spherical_harmonic_sine_coefficien_variations', 'total_spherical_harmonic_sine_coefficien_variations_from_indices', 'total_torque', 'total_torque_norm', 'total_torque_norm_type', 'total_torque_type', 'vehicle_panel_surface_normals_body_fixed_frame', 'vehicle_panel_surface_normals_inertial_frame', 'visible_radiation_source_area']

class PropagationDependentVariables:
    """Enumeration of available propagation dependent variables.
    
             Enumeration of propagation dependent variables supported by tudat.
    
    
    
    
    
          
    
    Members:
    
      mach_number_type : 
          
    
      altitude_type : 
          
    
      airspeed_type : 
          
    
      local_density_type : 
          
    
      relative_speed_type : 
          
    
      relative_position_type : 
          
    
      relative_distance_type : 
          
    
      relative_velocity_type : 
          
    
      radiation_pressure_type : 
          
    
      total_acceleration_norm_type : 
          
    
      single_acceleration_norm_type : 
          
    
      total_acceleration_type : 
          
    
      single_acceleration_type : 
          
    
      aerodynamic_force_coefficients_type : 
          
    
      aerodynamic_moment_coefficients_type : 
          
    
      rotation_matrix_to_body_fixed_frame_type : 
          
    
      intermediate_aerodynamic_rotation_matrix_type : 
          
    
      relative_body_aerodynamic_orientation_angle_type : 
          
    
      body_fixed_airspeed_based_velocity_type : 
          
    
      total_aerodynamic_g_load_type : 
          
    
      stagnation_point_heat_flux_type : No documentation found.
    
      local_temperature_type : 
          
    
      geodetic_latitude_type : 
          
    
      control_surface_deflection_type : 
          
    
      total_mass_rate_type : 
          
    
      tnw_to_inertial_frame_rotation_type : 
          
    
      rsw_to_inertial_frame_rotation_type : 
          
    
      periapsis_altitude_type : 
          
    
      apoapsis_altitude_type : 
          
    
      total_torque_norm_type : 
          
    
      single_torque_norm_type : 
          
    
      total_torque_type : 
          
    
      single_torque_type : 
          
    
      body_fixed_groundspeed_based_velocity_type : 
          
    
      keplerian_state_type : 
          
    
      modified_equinoctial_state_type : 
          
    
      spherical_harmonic_acceleration_terms_type : 
          
    
      spherical_harmonic_acceleration_norm_terms_type : 
          
    
      body_fixed_relative_cartesian_position_type : 
          
    
      body_fixed_relative_spherical_position_type : 
          
    
      total_gravity_field_variation_acceleration_type : 
          
    
      single_gravity_field_variation_acceleration_type : 
          
    
      single_gravity_field_variation_acceleration_terms_type : 
          
    
      acceleration_partial_wrt_body_translational_state_type : 
          
    
      local_dynamic_pressure_type : 
          
    
      euler_angles_to_body_fixed_type : 
          
    
      current_body_mass_type : 
          
    
      radiation_pressure_coefficient_type : 
          
    
      custom_type : No documentation found.
    
      gravity_field_potential_type : 
          
    
      gravity_field_laplacian_of_potential_type : 
          """
    __members__: typing.ClassVar[dict[str, PropagationDependentVariables]]
    acceleration_partial_wrt_body_translational_state_type: typing.ClassVar[PropagationDependentVariables]
    aerodynamic_force_coefficients_type: typing.ClassVar[PropagationDependentVariables]
    aerodynamic_moment_coefficients_type: typing.ClassVar[PropagationDependentVariables]
    airspeed_type: typing.ClassVar[PropagationDependentVariables]
    altitude_type: typing.ClassVar[PropagationDependentVariables]
    apoapsis_altitude_type: typing.ClassVar[PropagationDependentVariables]
    body_fixed_airspeed_based_velocity_type: typing.ClassVar[PropagationDependentVariables]
    body_fixed_groundspeed_based_velocity_type: typing.ClassVar[PropagationDependentVariables]
    body_fixed_relative_cartesian_position_type: typing.ClassVar[PropagationDependentVariables]
    body_fixed_relative_spherical_position_type: typing.ClassVar[PropagationDependentVariables]
    control_surface_deflection_type: typing.ClassVar[PropagationDependentVariables]
    current_body_mass_type: typing.ClassVar[PropagationDependentVariables]
    custom_type: typing.ClassVar[PropagationDependentVariables]
    euler_angles_to_body_fixed_type: typing.ClassVar[PropagationDependentVariables]
    geodetic_latitude_type: typing.ClassVar[PropagationDependentVariables]
    gravity_field_laplacian_of_potential_type: typing.ClassVar[PropagationDependentVariables]
    gravity_field_potential_type: typing.ClassVar[PropagationDependentVariables]
    intermediate_aerodynamic_rotation_matrix_type: typing.ClassVar[PropagationDependentVariables]
    keplerian_state_type: typing.ClassVar[PropagationDependentVariables]
    local_density_type: typing.ClassVar[PropagationDependentVariables]
    local_dynamic_pressure_type: typing.ClassVar[PropagationDependentVariables]
    local_temperature_type: typing.ClassVar[PropagationDependentVariables]
    mach_number_type: typing.ClassVar[PropagationDependentVariables]
    modified_equinoctial_state_type: typing.ClassVar[PropagationDependentVariables]
    periapsis_altitude_type: typing.ClassVar[PropagationDependentVariables]
    radiation_pressure_coefficient_type: typing.ClassVar[PropagationDependentVariables]
    radiation_pressure_type: typing.ClassVar[PropagationDependentVariables]
    relative_body_aerodynamic_orientation_angle_type: typing.ClassVar[PropagationDependentVariables]
    relative_distance_type: typing.ClassVar[PropagationDependentVariables]
    relative_position_type: typing.ClassVar[PropagationDependentVariables]
    relative_speed_type: typing.ClassVar[PropagationDependentVariables]
    relative_velocity_type: typing.ClassVar[PropagationDependentVariables]
    rotation_matrix_to_body_fixed_frame_type: typing.ClassVar[PropagationDependentVariables]
    rsw_to_inertial_frame_rotation_type: typing.ClassVar[PropagationDependentVariables]
    single_acceleration_norm_type: typing.ClassVar[PropagationDependentVariables]
    single_acceleration_type: typing.ClassVar[PropagationDependentVariables]
    single_gravity_field_variation_acceleration_terms_type: typing.ClassVar[PropagationDependentVariables]
    single_gravity_field_variation_acceleration_type: typing.ClassVar[PropagationDependentVariables]
    single_torque_norm_type: typing.ClassVar[PropagationDependentVariables]
    single_torque_type: typing.ClassVar[PropagationDependentVariables]
    spherical_harmonic_acceleration_norm_terms_type: typing.ClassVar[PropagationDependentVariables]
    spherical_harmonic_acceleration_terms_type: typing.ClassVar[PropagationDependentVariables]
    stagnation_point_heat_flux_type: typing.ClassVar[PropagationDependentVariables]
    tnw_to_inertial_frame_rotation_type: typing.ClassVar[PropagationDependentVariables]
    total_acceleration_norm_type: typing.ClassVar[PropagationDependentVariables]
    total_acceleration_type: typing.ClassVar[PropagationDependentVariables]
    total_aerodynamic_g_load_type: typing.ClassVar[PropagationDependentVariables]
    total_gravity_field_variation_acceleration_type: typing.ClassVar[PropagationDependentVariables]
    total_mass_rate_type: typing.ClassVar[PropagationDependentVariables]
    total_torque_norm_type: typing.ClassVar[PropagationDependentVariables]
    total_torque_type: typing.ClassVar[PropagationDependentVariables]

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

class SingleAccelerationDependentVariableSaveSettings(SingleDependentVariableSaveSettings):
    """`SingleDependentVariableSaveSettings`-derived class to save a single acceleration (norm or vector) during propagation.
    
    Class to define settings for saving a single acceleration (norm or vector) during propagation. Note: this acceleration is returned in the inertial frame!"""

    @property
    def acceleration_model_type(self) -> acceleration.AvailableAcceleration:
        """
        The type of the acceleration that is to be saved.
        
        :type: AvailableAcceleration
        """

class SingleDependentVariableSaveSettings(VariableSettings):
    """`VariableSettings`-derived class to define settings for dependent variables that are to be saved during propagation.
    
    Functional base class for defining settings for dependent variables that are to be computed and saved during propagation.
    Any dependent variable that requires additional information in addition to what can be provided here, should be
    defined by a dedicated derived class."""

    @property
    def associated_body(self) -> str:
        """
        Body associated with variable.
        
        :type: str
        """

    @property
    def component_index(self) -> int:
        """
        Index of the component to be saved.
        
        .. note::
        
            This parameter is only applicable for vectorial dependent variables. If negative, all components of the vector are saved.                    
                                            
        :type: int
        """

    @property
    def dependent_variable_type(self) -> PropagationDependentVariables:
        """
        Type of the dependent variable that is being saved.
        
        :type: PropagationDependentVariables
        """

    @property
    def secondary_body(self) -> str:
        """
        Secondary body w.r.t. which parameter is defined, e.g. relative position is defined for the associated body w.r.t. secondary body.
        
        .. note::
        
            This parameter is not necessarily defined. In this case an empty string is returned.                    
                                            
        :type: str
        """

class VariableSettings:
    """Functional base class to define settings for variables.
    
    This class is a functional base class for defining settings for variables.
    Any variable that requires additional information in addition to what can be provided here, should be defined by a
    dedicated derived class."""

def actual_cross_section(target_name: str, central_body_name: str, acceleration_type: str='radiation_pressure') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def aerodynamic_coefficients(target_name: str, central_body_name: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def aerodynamic_force_coefficients(body: str, central_body: str='') -> SingleDependentVariableSaveSettings:
    """Function to add the aerodynamic force coefficients to the dependent variables to save.
    
    Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic coefficient interface to be defined for the vehicle. The coefficients are returned in the following order: C_D, C_S, C_l (if coefficient interface defined in aerodynamic frame), or C_X, C_Y, C_Z (if coefficient interface defined in body frame).
    
    Parameters
    ----------
    body : str
        Body undergoing acceleration.
    central_body : str
        Body exerting acceleration (e.g. body with atmosphere).
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def aerodynamic_force_coefficients_control_surface_free(body: str, central_body: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def aerodynamic_force_coefficients_control_surface_increment(body: str, control_surface_name: str, central_body: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def aerodynamic_moment_coefficients(body: str, central_body: str='') -> SingleDependentVariableSaveSettings:
    """Function to add the aerodynamic moment coefficients to the dependent variables to save.
    
    Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic coefficient interface to be defined for the vehicle. The coefficients are returned in the following order: C_l, C_m, C_n , respectively about the X, Y, Z axes of the body-fixed frame, see (see :cite:t:`mooij1994`)
    
    Parameters
    ----------
    body : str
        Body undergoing acceleration.
    central_body : str
        Body exerting acceleration (e.g. body with atmosphere).
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def aerodynamic_moment_coefficients_control_surface_free(body: str, central_body: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def aerodynamic_moment_coefficients_control_surface_increment(body: str, control_surface_name: str, central_body: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def airspeed(body: str, body_with_atmosphere: str) -> SingleDependentVariableSaveSettings:
    """Function to add the airspeed to the dependent variables to save.
    
    Function to add the airspeed to the dependent variables to save. The calculation of the airspeed uses the rotation and wind models of the central body (to determine the motion of the atmosphere in inertial space), and the current state of the body for which the airspeed is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with atmosphere with respect to which the airspeed is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def altitude(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the altitude to the dependent variables to save.
    
    Function to add the altitude to the dependent variables to save. The calculation of the altitude uses the shape model of the central body and the current state of the body for which the altitude is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose altitude is to be saved.
    central_body : str
        Body with respect to which the altitude is computed (requires this body to have a shape model defined).
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def angle_of_attack(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the angle of attack to the dependent variables to save.
    
    Function to add the angle of attack angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the angle of attack is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def apoapsis_altitude(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the altitude of apoapsis to the dependent variables to save.
    
    Function to add the apoapsis altitude of the current osculating orbit to the dependent variables to save. The altitude depends on the shape of the central body. This function takes the current (osculating) orbit of the body w.r.t. the central body, and uses this Kepler orbit to extract the position/altitude of apoapsis.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the altitude of apoapsis is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def bank_angle(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the bank angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.
    
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the bank angle is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def body_fixed_airspeed_velocity(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the airspeed velocity vector to the dependent variables to save.
    
    Function to add the airspeed velocity vector to the dependent variables to save. The airspeed velocity vector is *not provided in an inertial frame*, but instead a frame centered on, and fixed to, the central body. It defines the velocity vector of a body w.r.t. the relative atmosphere It requires the central body to have an atmosphere.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the airspeed is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def body_fixed_groundspeed_velocity(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the groundspeed velocity vector to the dependent variables to save.
    
    Function to add the groundspeed velocity vector to the dependent variables to save. The groundspeed velocity vector is *not provided in an inertial frame*, but instead a frame centered on, and fixed to, the central body. It defines the velocity vector of a body w.r.t. 'the ground' or (alternatively and identically) the relative atmosphere in the case the atmosphere would be perfectly co-rotating with the central body.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the groundspeed is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def body_mass(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the current body mass to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body whose mass should be saved.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def center_of_mass(body: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def central_body_fixed_cartesian_position(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the relative Cartesian position, in the central body's fixed frame, to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body whose relative cartesian position is to be saved.
    central_body : str
        Body with respect to which the cartesian, body-fixed is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def central_body_fixed_spherical_position(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the spherical, body-fixed position to the dependent variables to save.
    
    Function to add the spherical position to the dependent variables to save. The spherical position is return as the radius, latitude, longitude, defined in the body-fixed frame of the central body
    
    Parameters
    ----------
    body : str
        Body whose spherical position is to be saved.
    central_body : str
        Body with respect to which the spherical, body-fixed is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def control_surface_deflection(body: str, control_surface: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def cross_section_change(target_name: str, source_name: str, acceleration_type: str='radiation_pressure') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def custom(custom_function: typing.Callable[[], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], variable_size: int) -> SingleDependentVariableSaveSettings:
    ...

def custom_dependent_variable(custom_function: typing.Callable[[], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], variable_size: int) -> SingleDependentVariableSaveSettings:
    """Function to compute a custom dependent variable.
    
    Function to compute a custom dependent variable, which can be implemented by the user as a Python function. The custom dependent variable is typically dependent on the current properties of the environment (e.g. bodies in the environment) or a user-defined guidance class (or similar)
    
    Parameters
    ----------
    custom_function : Callable[[], numpy.ndarray].
        Function taking no input, and returning the custom dependent variable (as a numpy Nx1 array).
    variable_size : int
        Size N of the array returned by the ``custom_function``
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    ``variable_size``"""

def density(body: str, body_with_atmosphere: str) -> SingleDependentVariableSaveSettings:
    """Function to add the local freestream density to the dependent variables to save.
    
    Function to add the freestream density (at a body's position) to the dependent variables to save. The calculation of the density uses the atmosphere model of the central body, and the current state of the body for which the density is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    body_with_atmosphere : str
        Body with atmosphere with respect to which the density is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def dynamic_pressure(body: str, body_with_atmosphere: str) -> SingleDependentVariableSaveSettings:
    """Function to add the local freestream dynamic pressure to the dependent variables to save.
    
    Function to add the freestream dynamic pressure (at a body's position) to the dependent variables to save. The calculation of the temperature uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    body_with_atmosphere : str
        Body with atmosphere with respect to which the dynamic pressure is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def flight_path_angle(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the flight path angle to the dependent variables to save.
    
    Function to add the flight path angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the flight path angle is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def full_body_paneled_geometry(target_name: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def geodetic_latitude(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the geodetic latitude to the dependent variables to save.
    
    Function to add the geodetic latitude, in the body-fixed frame of a central body, to the dependent variables to save. If the central body has a spherical shape model, this value is identical to the latitude. If the central body has an oblate spheroid shape model, the calculation of the geodetic latitude uses the flattening of the this shape model to determine the geodetic latitude
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the geodetic latitude is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def get_dependent_variable_id(dependent_variable_settings: SingleDependentVariableSaveSettings) -> str:
    '''Function to return a string representation of the dependent variable.
    
    Parameters
    ----------
    dependent_variable_settings : SingleDependentVariableSaveSettings
        Dependent variable settings for which the representation is to be returned.
    Returns
    -------
    str
        String representation of the dependent variable.
    Examples
    --------
    
    .. code-block:: python
    
        from tudatpy.dynamics import propagation_setup
    
        dep_var_settings = propagation_setup.dependent_variable.single_acceleration(propagation_setup.acceleration.point_mass_gravity_type, "SC", "Earth")
    
        print(propagation_setup.dependent_variable.get_dependent_variable_id(dep_var_settings))
        # prints "Single acceleration in inertial frame of type central gravity , acting on SC, exerted by Earth"'''

def get_dependent_variable_shape(dependent_variable_settings: SingleDependentVariableSaveSettings, bodies: environment.SystemOfBodies) -> tuple[int, int]:
    """Function to get the shape of a dependent variable.
    
    Parameters
    ----------
    dependent_variable_settings : SingleDependentVariableSaveSettings
        Dependent variable settings for which the shape is to be returned.
    bodies : SystemOfBodies
        Bodies that contain the bodies associated with the dependent variable.
    Returns
    -------
    tuple[int, int]
        Shape of the dependent variable.
    Examples
    --------
    
    .. code-block:: python
    
        from tudatpy.dynamics import propagation_setup
    
        bodies = ...
    
        dep_var_settings = propagation_setup.dependent_variable.single_acceleration(propagation_setup.acceleration.point_mass_gravity_type, "SC", "Earth")
    
        print(propagation_setup.dependent_variable.get_dependent_variable_shape(dep_var_settings, bodies))
        # prints (3, 1)"""

def get_dependent_variable_size(dependent_variable_settings: SingleDependentVariableSaveSettings, bodies: environment.SystemOfBodies) -> int:
    """Function to get the size of a dependent variable.
    
    Parameters
    ----------
    dependent_variable_settings : SingleDependentVariableSaveSettings
        Dependent variable settings for which the size is to be returned.
    bodies : SystemOfBodies
        Bodies that contain the bodies associated with the dependent variable.
    Returns
    -------
    int
        Size of the dependent variable.
    Examples
    --------
    
    .. code-block:: python
    
        from tudatpy.dynamics import propagation_setup
    
        bodies = ...
    
        dep_var_settings = propagation_setup.dependent_variable.single_acceleration(propagation_setup.acceleration.point_mass_gravity_type, "SC", "Earth")
    
        print(propagation_setup.dependent_variable.get_dependent_variable_size(dep_var_settings, bodies))
        # prints 3"""

def gravity_field_laplacian_of_potential(body_undergoing_acceleration: str, body_exerting_acceleration: str) -> SingleDependentVariableSaveSettings:
    """Function to add the laplacian of the gravitational potential to the dependent variables to save.
    
    Function to add the laplacian of the gravitational potential to the dependent variables to save. The laplacian is defined by the bodies undergoing and exerting the acceleration.
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body whose dependent variable should be saved.
    body_exerting_acceleration : str
        Body exerting acceleration.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def gravity_field_potential(body_undergoing_acceleration: str, body_exerting_acceleration: str) -> SingleDependentVariableSaveSettings:
    """Function to add the gravitational potential to the dependent variables to save.
    
    Function to add the gravitational potential to the dependent variables to save. The gravitational potential is defined by the bodies undergoing and exerting the acceleration.
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body whose dependent variable should be saved.
    body_exerting_acceleration : str
        Body exerting acceleration.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def heading_angle(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the heading angle to the dependent variables to save.
    
    Function to add the heading angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the heading angle is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def illuminated_panel_fraction(target_name: str, source_name: str, panel_type: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def inertia_tensor(body: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def inertial_to_body_fixed_313_euler_angles(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the 3-1-3 Euler angles for the rotation from inertial to body-fixed frame to the dependent variables to save.
    
    Function to add the 3-1-3 Euler angles for the rotation from inertial to body-fixed frame to the dependent variables to save. This requires the rotation of the body to be defined (either in the environment or the state vector).
    
    Parameters
    ----------
    body : str
        Body for which the rotation angles are to be saved.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def inertial_to_body_fixed_rotation_frame(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the rotation matrix from inertial to body-fixed frame to the dependent variables to save.
    
    Function to add the rotation matrix from inertial to body-fixed frame to the dependent variables to save. This requires the rotation of the body to be defined (either in the environment or the state vector). 
    
    .. note:: 
    
        A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).
    
        You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
        
        .. code-block:: python
    
           import numpy as np
    
           # Index where the rotation matrix starts in the dependent variable history
           start_index = ...
           rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]
    
           rotation_matrix_history = rotation_matrix_components.reshape(
               (rotation_matrix_components.shape[0], 3, 3)
           )
    
    Parameters
    ----------
    body : str
        Body for which the rotation matrix is to be saved.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    9"""

def intermediate_aerodynamic_rotation_matrix_variable(body: str, base_frame: environment_setup.aerodynamic_coefficients.AerodynamicsReferenceFrames, target_frame: environment_setup.aerodynamic_coefficients.AerodynamicsReferenceFrames, central_body: str='') -> SingleDependentVariableSaveSettings:
    """Function to add the rotation matrix between any two reference frames used in aerodynamic calculations.
    
    Function to add the rotation matrix between any two reference frames used in aerodynamic calculations. The list of available frames is defined by the :class:`AerodynamicsReferenceFrames` enum.
    
    .. note:: 
    
        A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).
    
        You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
        
        .. code-block:: python
    
           import numpy as np
    
           # Index where the rotation matrix starts in the dependent variable history
           start_index = ...
           rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]
    
           rotation_matrix_history = rotation_matrix_components.reshape(
               (rotation_matrix_components.shape[0], 3, 3)
           )
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    base_frame : AerodynamicsReferenceFrames
        Base reference frame for the rotation.
    target_frame : AerodynamicsReferenceFrames
        Target reference frame for the rotation.
    central_body : str
        Central body w.r.t. which the state of the body is considered.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    9"""

def keplerian_state(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the Keplerian state to the dependent variables to save.
    
    Function to add the Keplerian state to the dependent variables to save. The Keplerian state is returned in this order: 1: Semi-major Axis. 2: Eccentricity. 3: Inclination. 4: Argument of Periapsis. 5. Right Ascension of the Ascending Node. 6: True Anomaly.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the Keplerian state is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    6"""

def latitude(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the latitude to the dependent variables to save.
    
    Function to add the latitude of a body, in the body-fixed frame of a central body, to the dependent variables to save.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the latitude is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def local_aerodynamic_g_load(body: str, body_with_atmosphere: str) -> SingleDependentVariableSaveSettings:
    """Function to add the total aerodynamic G-load to the dependent variables to save.
    
    Function to add the total aerodynamic G-load of a body to the dependent variables to save. The calculation uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    body_with_atmosphere : str
        Body with atmosphere exerting the aerodynamic acceleration for which the g-load is to be computed
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def longitude(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the longitude to the dependent variables to save.
    
    Function to add the longitude of a body, in the body-fixed frame of a central body, to the dependent variables to save.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the longitude is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def mach_number(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the Mach number to the dependent variables to save.
    
    Function to add the Mach number to the dependent variables to save. The calculation of the altitude uses the atmosphere model of the central body and the current state of the body for which the Mach number is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose Mach number is to be saved.
    central_body : str
        Body with atmosphere with respect to which the Mach number is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1
    
    
    
    
    Examples
    --------
    
    To create settings for saving of a Mach number of a body name 'Spacecraft'
    w.r.t. the atmosphere of body 'Earth', use:
    
    .. code-block:: python
    
       # Define save settings for Mach number
       propagation_setup.dependent_variable.mach_number( "Spacecraft", "Earth" )"""

def minimum_body_distance(body_name: str, bodies_to_check: list[str]) -> SingleDependentVariableSaveSettings:
    """Function to compute the minimum distance between a given body, and a set of other bodies.
    
    Function to compute the minimum distance between a given body, and a set of other bodies. This function takes the instantaneous position of body ``body_name``, and each body in the list ``bodies_to_check``, and computes the body from this list closest to ``body_name``. In this calculation, the positions of the bodies are evaluated at the current propagation time, and therefore **light time is ignored**. In addition, this functions does not consider visibility requirements (e.g. is a planet between two bodies). The dependent variable is of size 2, and consists of: (0) The distance between the body, and the closest other body; (1) The index from ``bodies_to_check`` for which the distance (given by the first index) is closest to ``body`` Typically, this function is used to compute the closest body in a constellation of satellites.
    
    Parameters
    ----------
    body_name : str
        Body for which the distance to other bodies is to be computed.
    bodies_to_check : list[ str ]
        List of bodies for which it is to be checked which of these bodies is closest to ``body_name``.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    2"""

def minimum_visible_station_body_distances(body_name: str, station_name: str, bodies_to_check: list[str], minimum_elevation_angle: float) -> SingleDependentVariableSaveSettings:
    """Function to compute the minimum distance between a ground station, and a set of other bodies visible from that station.
    
    Function to compute the minimum distance between a ground station, and a set of other bodies visible from that station This function takes the instantaneous position of the ground station ``station_name`` on ``body_name``, and each body in the list ``bodies_to_check``, and computes the body from this list closest to this ground station, only taking into account those bodies from this list which are visible from teh ground station. For this function, visibility is defined by a single elevation angle cutoff (at the ground station) below which a body is deemed to not be visible. In this calculation, the positions of the bodies are evaluated at the current propagation time, and therefore **light time is ignored**. The dependent variable is of size 3, and consists of: (0) The distance between the ground station, and the closest visible body; (1) The index from ``bodies_to_check`` for which the distance (given by the first index) is closest to thee ground station, and the body is visible. (2) Elevation angle for closest body. In case, no body is visible from the station, this function returns [NaN, -1, NaN]. Typically, this function is used to compute the closest body between a ground station and a constellation of satellites.
    
    Parameters
    ----------
    body_name : str
        Body on which ground station is located, for which the distance to other bodies is to be computed.
    station_name : str
        Name of ground station, for which the distance to other bodies is to be computed.
    bodies_to_check : list[ str ]
        List of bodies for which it is to be checked which of these bodies is closest to ``station_name`` on ``body_name``.
    minimum_elevation_angle : float
        Minimum elevation angle (at ground station) below which the distance to the ``bodies_to_check`` is not considered.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def modified_equinoctial_state(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the modified equinoctial state to the dependent variables to save.
    
    Function to add the modified equinoctial state to the dependent variables to save. The value of the parameter I is automatically chosen as +1 or -1, depending on whether the inclination is smaller or larger than 90 degrees. The elements are returned in the order :math:`p`, :math:`f`, :math:`g`, :math:`h`, :math:`k`, :math:`L`
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the modified equinoctial state is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    6"""

def per_target_panel_radiation_pressure_force(target_name: str, source_name: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def periapsis_altitude(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the altitude of periapsis to the dependent variables to save.
    
    Function to add the periapsis altitude of the current osculating orbit to the dependent variables to save. The altitude depends on the shape of the central body. This function takes the current (osculating) orbit of the body w.r.t. the central body, and uses this Kepler orbit to extract the position/altitude of periapsis.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the altitude of periapsis is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def radiation_pressure(target_body: str, source_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the radiation pressure to the dependent variables to save.
    
    Function to add the local radiation pressure, in N/m^2, to the dependent variables to save. It requires a radiation source model to be defined for the radiating body.
    
    Parameters
    ----------
    target_body : str
        Name of body at the location of which the radiation pressure is to be returned
    source_body : str
        Name of body from which the radiation originates that causes the radiation pressure at the target body
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def radiation_pressure_coefficient(body: str, emitting_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the current radiation pressure coefficient to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    emitting_body : str
        Emitting body.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def radiation_pressure_source_panel_geometry(target_name: str, source_name: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def radiation_pressure_source_panel_irradiance(target_name: str, source_name: str) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def received_irradiance(target_body: str, source_body: str) -> SingleDependentVariableSaveSettings:
    """Function to save the received irradiance from a give source.
    
    Function to save the received irradiance (in W/m^{2}) from a given source, as used in the computation of radiation pressure
    
    Parameters
    ----------
    target_body : str
        Name of body at the location of which the irradiance is to be returned
    source_body : str
        Name of body from which the radiation originates that causes the irradiance at the target body
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def received_irradiance_shadow_function(target_body: str, source_body: str) -> SingleDependentVariableSaveSettings:
    """Function to save the shadow function that reduces the received irradiance from a given source.
    
    Function to save the shadow function that reduces the received irradiance from a given source, as a result of occultation of radiation between the source and the target. The shadow function is by definition between 0 (no irradiance is received, total occultation) and 1 (all irradiance is received, no occultation)
    
    Parameters
    ----------
    target_body : str
        Name of body at the location of which the irradiance is to be returned
    source_body : str
        Name of body from which the radiation originates that causes the irradiance at the target body
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def relative_distance(body: str, relative_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the relative distance to the dependent variables to save.
    
    Function to add a body's relative distance (norm of the position vector) with respect to a second body to the dependent variables to save. The relative distance is computed between the bodies' centers of mass.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    relative_body : str
        Body with respect to which the relative distance is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def relative_position(body: str, relative_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the relative position vector to the dependent variables to save.
    
    Function to add a body's relative position vector with respect to a second body to the dependent variables to save. The relative position is computed between the bodies' centers of mass.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    relative_body : str
        Body with respect to which the relative position is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def relative_speed(body: str, relative_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the relative speed to the dependent variables to save.
    
    Function to add a body's relative speed (norm of the relative velocity vector) with respect to a second body to the dependent variables to save. The relative speed is computed between the bodies' centers of mass.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    relative_body : str
        Body with respect to which the relative speed is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def relative_velocity(body: str, relative_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the relative velocity vector to the dependent variables to save.
    
    Function to add a body's relative velocity vector with respect to a second body to the dependent variables to save. The relative velocity is computed between the bodies' centers of mass.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    relative_body : str
        Body with respect to which the relative velocity is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def rsw_to_inertial_rotation_matrix(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the rotation matrix from the RSW to the inertial frame to the dependent variables to save.
    
    Function to add the rotation matrix from the RSW to the inertial frame to the dependent variables to save. It has the x-axis pointing along the position vector (away from the central body), the z-axis along the orbital angular momentum vector, and the y-axis completing the right-handed system.
    
    .. note:: 
    
        A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).
    
        You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
        
        .. code-block:: python
    
           import numpy as np
    
           # Index where the rotation matrix starts in the dependent variable history
           start_index = ...
           rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]
    
           rotation_matrix_history = rotation_matrix_components.reshape(
               (rotation_matrix_components.shape[0], 3, 3)
           )
    
    Parameters
    ----------
    body : str
        Body for which the rotation matrix is to be saved.
    central_body : str
        Body with respect to which the TNW frame is determined.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    9"""

def sideslip_angle(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the sideslip angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.
    
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    central_body : str
        Body with respect to which the sideslip angle is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def single_acceleration(acceleration_type: acceleration.AvailableAcceleration, body_undergoing_acceleration: str, body_exerting_acceleration: str) -> SingleAccelerationDependentVariableSaveSettings:
    """Function to add a single acceleration to the dependent variables to save.
    
    Function to add a single acceleration vector to the dependent variables to save. The requested acceleration is defined by its type, and the bodies undergoing and exerting the acceleration. This acceleration vector represents the acceleration in 3D in the inertial reference frame. NOTE: When requesting a third-body perturbation be saved, you may use either the direct acceleration type, or the third body type. For instance, for saving a point-mass third-body perturbation, you may specify either ``point_mass_gravity_type`` or ``third_body_point_mass_gravity_type`` as acceleration type.
    
    Parameters
    ----------
    acceleration_type : AvailableAcceleration
        Acceleration type to be saved.
    body_undergoing_acceleration : str
        Body undergoing acceleration.
    body_exerting_acceleration : str
        Body exerting acceleration.
    Returns
    -------
    SingleAccelerationDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3
    
    
    
    
    Examples
    --------
    
    To create settings for saving a point mass acceleration acting on body called 'Spacecraft',
    exerted by a body named 'Earth', use:
    
    .. code-block:: python
    
       # Define save settings for point-mass acceleration on Spacecraft by Earth
       propagation_setup.dependent_variable.single_acceleration(
               propagation_setup.acceleration.point_mass_gravity_type, 'Spacecraft', 'Earth' )"""

def single_acceleration_norm(acceleration_type: acceleration.AvailableAcceleration, body_undergoing_acceleration: str, body_exerting_acceleration: str) -> SingleAccelerationDependentVariableSaveSettings:
    """Function to add a single scalar acceleration to the dependent variables to save.
    
    Function to add a single scalar acceleration (norm of the acceleration vector) to the dependent variables to save. The requested acceleration is defined by its type, and the bodies undergoing and exerting the acceleration. NOTE: When requesting a third-body perturbation be saved, you may use either the direct acceleration type, or the third body type. For instance, for saving a point-mass third-body perturbation, you may specify either ``point_mass_gravity_type`` or ``third_body_point_mass_gravity_type`` as acceleration type.
    
    Parameters
    ----------
    acceleration_type : AvailableAcceleration
        Acceleration type to be saved
    body_undergoing_acceleration : str
        Body undergoing acceleration.
    body_exerting_acceleration : str
        Body exerting acceleration.
    Returns
    -------
    SingleAccelerationDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1
    
    
    
    
    Examples
    --------
    
    To create settings for saving norm of a point mass acceleration acting on body called 'Spacecraft',
    exerted by a body named 'Earth', use:
    
    .. code-block:: python
    
       # Define save settings for point-mass acceleration on Spacecraft by Earth
       propagation_setup.dependent_variable.single_acceleration_norm(
               propagation_setup.acceleration.point_mass_gravity_type, 'Spacecraft', 'Earth' )"""

def single_gravity_field_variation_acceleration(body_undergoing_acceleration: str, body_exerting_acceleration: str, deformation_type: environment_setup.gravity_field_variation.BodyDeformationTypes, identifier: str='') -> SingleDependentVariableSaveSettings:
    """Function to add the acceleration induced by a single time-variability of a gravity field to the dependent variables to save.
    
    Function to add the acceleration induced by a single time-variability of a gravity field to the dependent variables to save. The user specifies the type of variability for which the induced acceleration is to be saved.
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body whose dependent variable should be saved.
    body_exerting_acceleration : str
        Body exerting the acceleration.
    deformation_type : BodyDeformationTypes
        Type of gravity field variation for which the acceleration contribution is to be saved
    identifier : str, default=""
        Identifier for the deformation type. To be used in case multiple realizations of a single variation type are present in the given body. Otherwise, this entry can be left empty
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def single_per_term_gravity_field_variation_acceleration(body_undergoing_acceleration: str, body_exerting_acceleration: str, component_indices: list[tuple[int, int]], deformation_type: environment_setup.gravity_field_variation.BodyDeformationTypes, identifier: str='') -> SingleDependentVariableSaveSettings:
    """Function to add the acceleration induced by a single time-variability of a gravity field, at a given list of degrees/orders, to the dependent variables to save. This combines the functionality of the :func:`single_gravity_field_variation_acceleration` and :func:`spherical_harmonic_terms_acceleration` variables
    
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body whose dependent variable should be saved.
    body_exerting_acceleration : str
        Body exerting the acceleration.
    component_indices : list[tuple]
        Tuples of (degree, order) indicating the terms to save.
    deformation_type : BodyDeformationTypes
        Type of gravity field variation for which the acceleration contribution is to be saved
    identifier : str, default=""
        Identifier for the deformation type. To be used in case multiple realizations of a single variation type are present in the given body. Otherwise, this entry can be left empty
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3*``len(component_indices)``"""

def single_torque(torque_type: torque.AvailableTorque, body_undergoing_torque: str, body_exerting_torque: str) -> SingleDependentVariableSaveSettings:
    """Function to add a single torque vector to the dependent variables to save.
    
    
    Parameters
    ----------
    torque_type : AvailableTorque
        Torque type to be saved.
    body_undergoing_torque : str
        Body undergoing torque.
    body_exerting_torque : str
        Body exerting torque.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def single_torque_norm(torque_type: torque.AvailableTorque, body_undergoing_torque: str, body_exerting_torque: str) -> SingleDependentVariableSaveSettings:
    """Function to add a single torque (norm of the torque vector) to the dependent variables to save.
    
    
    Parameters
    ----------
    torque_type : AvailableTorque
        Torque type to be saved.
    body_undergoing_torque : str
        Body undergoing torque.
    body_exerting_torque : str
        Body exerting torque.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def spherical_harmonic_terms_acceleration(body_undergoing_acceleration: str, body_exerting_acceleration: str, component_indices: list[tuple[int, int]]) -> SingleDependentVariableSaveSettings:
    """Function to add single degree/order contributions of a spherical harmonic acceleration vector to the dependent variables to save.
    
    Function to add single degree/order contributions of a spherical harmonic acceleration vector to the dependent variables to save. The spherical harmonic acceleration consists of a (truncated) summation of contributions at degree :math:`l` and order :math:`m`. Using this function, you can save the contributions of separate :math:`l,m` entries to the total acceleration. For instance, when requesting dependent variables for :math:`l,m=2,2`, the contribution due to the combined influence of :math:`\\bar{C}_{22}` and `\\bar{S}_{22}` are provided
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body undergoing acceleration.
    body_exerting_acceleration : str
        Body exerting acceleration.
    component_indices : list[tuple]
        Tuples of (degree, order) indicating the terms to save.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3*``len(component_indices)``
    
    
    
    
    Examples
    --------
    
    To create settings for saving spherical harmonic acceleration contributions of degree/order
    2/0, 2/1 and 2/2, acting on a body names 'Spacecraft', exerted by a body named 'Earth',
    use the following for the acceleration. The resulting dependent variable will contain
    nine entries (three acceleration components for 2/0, 2/1 and 2/2, respectively).
    
    .. code-block:: python
    
       # Define degree/order combinations for which to save acceleration contributions
       spherical_harmonic_terms = [ (2,0), (2,1), (2,2) ]
    
       # Define save settings for separate spherical harmonic contributions
       propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration( "Spacecraft", "Earth", spherical_harmonic_terms )"""

def spherical_harmonic_terms_acceleration_norm(body_undergoing_acceleration: str, body_exerting_acceleration: str, component_indices: list[tuple[int, int]]) -> SingleDependentVariableSaveSettings:
    """Function to add a single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.
    
    Function to add single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body undergoing acceleration.
    body_exerting_acceleration : str
        Body exerting acceleration.
    component_indices : list[tuple]
        Tuples of (degree, order) indicating the terms to save.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    ``len(component_indices)``
    
    
    
    
    Examples
    --------
    
    To create settings for saving spherical harmonic acceleration contributions of degree/order
    2/0, 2/1 and 2/2, acting on a body names 'Spacecraft', exerted by a body named 'Earth',
    use the following for the acceleration. The resulting dependent variable will contain
    three entries (one acceleration norm for 2/0, 2/1 and 2/2, respectively).
    
    .. code-block:: python
    
    
       # Define degree/order combinations for which to save acceleration contributions
       spherical_harmonic_terms = [ (2,0), (2,1), (2,2) ]
    
       # Define save settings for separate spherical harmonic contributions
       propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration_norm( "Spacecraft", "Earth", spherical_harmonic_terms )"""

def temperature(body: str, body_with_atmosphere: str) -> SingleDependentVariableSaveSettings:
    """Function to add the local freestream temperature to the dependent variables to save.
    
    Function to add the freestream temperature (at a body's position) to the dependent variables to save. The calculation of the temperature uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    body_with_atmosphere : str
        Body with atmosphere with respect to which the temperature is computed.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def tnw_to_inertial_rotation_matrix(body: str, central_body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the rotation matrix from the TNW to the inertial frame to the dependent variables to save.
    
    Function to add the rotation matrix from the TNW to the inertial frame to the dependent variables to save. It has the x-axis pointing along the velocity vector, the z-axis along the orbital angular momentum vector, and the y-axis completing the right-handed system.
    
    .. note:: 
    
        A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).
    
        You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
        
        .. code-block:: python
    
           import numpy as np
    
           # Index where the rotation matrix starts in the dependent variable history
           start_index = ...
           rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]
    
           rotation_matrix_history = rotation_matrix_components.reshape(
               (rotation_matrix_components.shape[0], 3, 3)
           )
    
    Parameters
    ----------
    body : str
        Body for which the rotation matrix is to be saved.
    central_body : str
        Body with respect to which the TNW frame is determined.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    9"""

def total_acceleration(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the total acceleration vector acting on a body to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body undergoing acceleration.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def total_acceleration_norm(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the total scalar acceleration (norm of the vector) acting on a body to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body undergoing acceleration.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def total_gravity_field_variation_acceleration(body_undergoing_acceleration: str, body_exerting_acceleration: str) -> SingleDependentVariableSaveSettings:
    """Function to add the acceleration induced by the total time-variability of a gravity field to the dependent variables to save.
    
    Function to add the acceleration induced by the total time-variability of a gravity field to the dependent variables to save. This function does not distinguish between different sources of variations of the gravity field, and takes the full time-variation when computing the contribution to the acceleration. To select only one contribution, use the :func:`single_gravity_field_variation_acceleration` function.
    
    Parameters
    ----------
    body_undergoing_acceleration : str
        Body whose dependent variable should be saved.
    body_exerting_acceleration : str
        Body exerting the acceleration.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def total_mass_rate(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the total mass rate to the dependent variables to save.
    
    Function to add the total mass rate to the dependent variables to save. It requires the body mass to be numerically propagated.
    
    Parameters
    ----------
    body : str
        Body whose mass rate should be saved.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def total_spherical_harmonic_cosine_coefficien_variations(body: str, minimum_degree: int, maximum_degree: int, minimum_order: int, maximum_order: int) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def total_spherical_harmonic_cosine_coefficien_variations_from_indices(body: str, component_indices: list[tuple[int, int]]) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def total_spherical_harmonic_sine_coefficien_variations(body: str, minimum_degree: int, maximum_degree: int, minimum_order: int, maximum_order: int) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def total_spherical_harmonic_sine_coefficien_variations_from_indices(body: str, component_indices: list[tuple[int, int]]) -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def total_torque(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the total torque vector to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    3"""

def total_torque_norm(body: str) -> SingleDependentVariableSaveSettings:
    """Function to add the total torque (norm of the torque vector) to the dependent variables to save.
    
    
    Parameters
    ----------
    body : str
        Body whose dependent variable should be saved.
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""

def vehicle_panel_surface_normals_body_fixed_frame(body_name: str, part_name: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def vehicle_panel_surface_normals_inertial_frame(body_name: str, part_name: str='') -> SingleDependentVariableSaveSettings:
    """No documentation found."""

def visible_radiation_source_area(target_body: str, source_body: str) -> SingleDependentVariableSaveSettings:
    """Function to save the shadow function that reduces the received irradiance from a given source.
    
    Function to save the shadow function that reduces the received irradiance from a given source, as a result of occultation of radiation between the source and the target. The shadow function is by definition between 0 (no irradiance is received, total occultation) and 1 (all irradiance is received, no occultation)
    
    Parameters
    ----------
    target_body : str
        Name of body at the location of which the irradiance is to be returned
    source_body : str
        Name of body from which the radiation originates that causes the irradiance at the target body
    Returns
    -------
    SingleDependentVariableSaveSettings
        Dependent variable settings object.
    Variable Size
    -------------
    1"""
acceleration_partial_wrt_body_translational_state_type: PropagationDependentVariables
aerodynamic_force_coefficients_type: PropagationDependentVariables
aerodynamic_moment_coefficients_type: PropagationDependentVariables
airspeed_type: PropagationDependentVariables
altitude_type: PropagationDependentVariables
apoapsis_altitude_type: PropagationDependentVariables
body_fixed_airspeed_based_velocity_type: PropagationDependentVariables
body_fixed_groundspeed_based_velocity_type: PropagationDependentVariables
body_fixed_relative_cartesian_position_type: PropagationDependentVariables
body_fixed_relative_spherical_position_type: PropagationDependentVariables
control_surface_deflection_type: PropagationDependentVariables
current_body_mass_type: PropagationDependentVariables
custom_type: PropagationDependentVariables
euler_angles_to_body_fixed_type: PropagationDependentVariables
geodetic_latitude_type: PropagationDependentVariables
gravity_field_laplacian_of_potential_type: PropagationDependentVariables
gravity_field_potential_type: PropagationDependentVariables
intermediate_aerodynamic_rotation_matrix_type: PropagationDependentVariables
keplerian_state_type: PropagationDependentVariables
local_density_type: PropagationDependentVariables
local_dynamic_pressure_type: PropagationDependentVariables
local_temperature_type: PropagationDependentVariables
mach_number_type: PropagationDependentVariables
modified_equinoctial_state_type: PropagationDependentVariables
periapsis_altitude_type: PropagationDependentVariables
radiation_pressure_coefficient_type: PropagationDependentVariables
radiation_pressure_type: PropagationDependentVariables
relative_body_aerodynamic_orientation_angle_type: PropagationDependentVariables
relative_distance_type: PropagationDependentVariables
relative_position_type: PropagationDependentVariables
relative_speed_type: PropagationDependentVariables
relative_velocity_type: PropagationDependentVariables
rotation_matrix_to_body_fixed_frame_type: PropagationDependentVariables
rsw_to_inertial_frame_rotation_type: PropagationDependentVariables
single_acceleration_norm_type: PropagationDependentVariables
single_acceleration_type: PropagationDependentVariables
single_gravity_field_variation_acceleration_terms_type: PropagationDependentVariables
single_gravity_field_variation_acceleration_type: PropagationDependentVariables
single_torque_norm_type: PropagationDependentVariables
single_torque_type: PropagationDependentVariables
spherical_harmonic_acceleration_norm_terms_type: PropagationDependentVariables
spherical_harmonic_acceleration_terms_type: PropagationDependentVariables
stagnation_point_heat_flux_type: PropagationDependentVariables
tnw_to_inertial_frame_rotation_type: PropagationDependentVariables
total_acceleration_norm_type: PropagationDependentVariables
total_acceleration_type: PropagationDependentVariables
total_aerodynamic_g_load_type: PropagationDependentVariables
total_gravity_field_variation_acceleration_type: PropagationDependentVariables
total_mass_rate_type: PropagationDependentVariables
total_torque_norm_type: PropagationDependentVariables
total_torque_type: PropagationDependentVariables