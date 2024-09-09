import numpy
import typing
__all__ = ['AerodynamicAngleCalculator', 'AerodynamicAngleRotationalEphemeris', 'AerodynamicCoefficientFrames', 'AerodynamicCoefficientGenerator36', 'AerodynamicCoefficientInterface', 'AerodynamicCoefficientsIndependentVariables', 'AerodynamicsReferenceFrameAngles', 'AerodynamicsReferenceFrames', 'AtmosphericFlightConditions', 'Body', 'ConstantEphemeris', 'ConstantFrequencyInterpolator', 'ControlSurfaceIncrementAerodynamicInterface', 'CustomBodyFixedDirectionCalculator', 'CustomControlSurfaceIncrementAerodynamicInterface', 'CustomInertialDirectionBasedRotationalEphemeris', 'DirectLongitudeLibrationCalculator', 'EngineModel', 'Ephemeris', 'FlightConditions', 'GcrsToItrsRotationModel', 'GravityFieldModel', 'GroundStation', 'GroundStationState', 'HypersonicLocalInclinationAnalysis', 'InertialBodyFixedDirectionCalculator', 'KeplerEphemeris', 'LongitudeLibrationCalculator', 'MultiArcEphemeris', 'PointingAnglesCalculator', 'PolyhedronGravityField', 'RotationalEphemeris', 'ShapeModel', 'SphericalHarmonicsGravityField', 'StationFrequencyInterpolator', 'SynchronousRotationalEphemeris', 'SystemOfBodies', 'TabulatedEphemeris', 'Tle', 'TleEphemeris', 'VehicleSystems', 'aerodynamic_frame', 'altitude_dependent', 'angle_of_attack', 'angle_of_attack_dependent', 'angle_of_sideslip', 'anomalous_o_number_density_dependent', 'ar_number_density_dependent', 'bank_angle', 'body_frame', 'control_surface_deflection_dependent', 'corotating_frame', 'flight_path_angle', 'get_default_local_inclination_angle_of_attack_points', 'get_default_local_inclination_mach_points', 'get_default_local_inclination_sideslip_angle_points', 'get_local_inclination_mesh', 'get_local_inclination_total_vehicle_area', 'h_number_density_dependent', 'he_number_density_dependent', 'heading_angle', 'inertial_frame', 'latitude_angle', 'longitude_angle', 'mach_number_dependent', 'n2_number_density_dependent', 'n_number_density_dependent', 'negative_aerodynamic_frame_coefficients', 'negative_body_fixed_frame_coefficients', 'o2_number_density_dependent', 'o_number_density_dependent', 'positive_aerodynamic_frame_coefficients', 'positive_body_fixed_frame_coefficients', 'save_vehicle_mesh_to_file', 'sideslip_angle_dependent', 'temperature_dependent', 'time_dependent', 'trajectory_frame', 'transform_to_inertial_orientation', 'undefined_independent_variable', 'velocity_dependent', 'vertical_frame']

class AerodynamicAngleCalculator:
    """Object to calculate (aerodynamic) orientation angles, and frame transformations,
	from current vehicle state.
	
	
		Object to calculate (aerodynamic) orientation angles (list given by the :class:`~AerodynamicsReferenceFrameAngles` enum)
		and transformations between frames (list given by the :class:`~AerodynamicsReferenceFrames` enum) from current vehicle state.
	"""

    def get_angle(self, angle_type: AerodynamicsReferenceFrameAngles) -> float:
        """
        Function to get a single orientation angle
        
        
        	Function to get a single orientation angle. This function
        	is meant to be used only *during* a numerical propagation, in particular
        	for the definition of a custom (e.g. guidance) model.
        
        
        	:param original_frame:
        		The identifier for the angle that is to be returnd
        
        	:return:
        		Value of requested angle
        """

    def get_rotation_matrix_between_frames(self, original_frame: AerodynamicsReferenceFrames, target_frame: AerodynamicsReferenceFrames) -> numpy.ndarray:
        """
        Function to get the rotation matrix between two frames.
        
        
        	Function to get the rotation matrix between two frames. This function
        	is meant to be used only *during* a numerical propagation, in particular
        	for the definition of a custom (e.g. guidance) model.
        
        
        	:param original_frame:
        		The frame :math:`A` from which the rotation matrix is to be calculated
        
        	:param target_frame:
        		The frame :math:`B` to which the rotation matrix is to be calculated
        
        	:return:
        		Rotation matrix :math:`\\mathbf{R}^{B/A}` from frame :math:`A` to frame `B`
        """

    def set_body_orientation_angle_functions(self, angle_of_attack_function: typing.Callable[[], float]=None, angle_of_sideslip_function: typing.Callable[[], float]=None, bank_angle_function: typing.Callable[[], float]=None, angle_update_function: typing.Callable[[float], None]=None, silence_warnings: bool=False) -> None:
        ...

    def set_body_orientation_angles(self, angle_of_attack: float=..., angle_of_sideslip: float=..., bank_angle: float=..., silence_warnings: bool=False) -> None:
        ...

class AerodynamicAngleRotationalEphemeris(RotationalEphemeris):

    def reset_aerodynamic_angle_function(self, arg0: typing.Callable[[float], numpy.ndarray]) -> None:
        ...

class AerodynamicCoefficientFrames:
    """Members:
	
	positive_body_fixed_frame_coefficients :
	
	negative_body_fixed_frame_coefficients :
	
	positive_aerodynamic_frame_coefficients :
	
	negative_aerodynamic_frame_coefficients :
	"""
    __members__: typing.ClassVar[dict[str, AerodynamicCoefficientFrames]]
    negative_aerodynamic_frame_coefficients: typing.ClassVar[AerodynamicCoefficientFrames]
    negative_body_fixed_frame_coefficients: typing.ClassVar[AerodynamicCoefficientFrames]
    positive_aerodynamic_frame_coefficients: typing.ClassVar[AerodynamicCoefficientFrames]
    positive_body_fixed_frame_coefficients: typing.ClassVar[AerodynamicCoefficientFrames]

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

class AerodynamicCoefficientGenerator36(AerodynamicCoefficientInterface):
    """<no_doc, only_dec>
	"""

class AerodynamicCoefficientInterface:
    """
		"""

    def current_control_surface_force_coefficient_increment(self, control_surface_name: str) -> numpy.ndarray:
        ...

    def current_control_surface_moment_coefficient_increment(self, control_surface_name: str) -> numpy.ndarray:
        ...

    def set_control_surface_increments(self, control_surface_list: dict[str, ...]) -> None:
        ...

    def update_coefficients(self, independent_variables: list[float], time: float) -> None:
        ...

    def update_full_coefficients(self, independent_variables: list[float], control_surface_independent_variables: dict[str, list[float]], time: float, check_force_contribution: bool=True) -> None:
        ...

    @property
    def control_surface_independent_variable_names(self) -> dict[str, list[AerodynamicCoefficientsIndependentVariables]]:
        ...

    @property
    def current_coefficients(self) -> numpy.ndarray:
        ...

    @property
    def current_control_surface_free_force_coefficients(self) -> numpy.ndarray:
        ...

    @property
    def current_control_surface_free_moment_coefficients(self) -> numpy.ndarray:
        ...

    @property
    def current_force_coefficients(self) -> numpy.ndarray:
        ...

    @property
    def current_moment_coefficients(self) -> numpy.ndarray:
        ...

    @property
    def force_coefficient_frame(self) -> AerodynamicCoefficientFrames:
        ...

    @property
    def independent_variable_names(self) -> list[AerodynamicCoefficientsIndependentVariables]:
        ...

    @property
    def moment_coefficient_frame(self) -> AerodynamicCoefficientFrames:
        ...

    @property
    def reference_area(self) -> float:
        ...

class AerodynamicCoefficientsIndependentVariables:
    """Enumeration of the independent variables that can be used to compute aerodynamic coefficients.
	
	
	
	:member mach_number_dependent:
	:member angle_of_attack_dependent:
	:member sideslip_angle_dependent:
	:member altitude_dependent:
	:member time_dependent:
	:member control_surface_deflection_dependent:
	:member undefined_independent_variable:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, AerodynamicCoefficientsIndependentVariables]]
    altitude_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    angle_of_attack_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    anomalous_o_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    ar_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    control_surface_deflection_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    h_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    he_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    mach_number_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    n2_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    n_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    o2_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    o_number_density_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    sideslip_angle_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    temperature_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    time_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    undefined_independent_variable: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]
    velocity_dependent: typing.ClassVar[AerodynamicCoefficientsIndependentVariables]

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

class AerodynamicsReferenceFrameAngles:
    """Members:
	
	latitude_angle
	
	longitude_angle
	
	heading_angle
	
	flight_path_angle
	
	angle_of_attack
	
	angle_of_sideslip
	
	bank_angle
	"""
    __members__: typing.ClassVar[dict[str, AerodynamicsReferenceFrameAngles]]
    angle_of_attack: typing.ClassVar[AerodynamicsReferenceFrameAngles]
    angle_of_sideslip: typing.ClassVar[AerodynamicsReferenceFrameAngles]
    bank_angle: typing.ClassVar[AerodynamicsReferenceFrameAngles]
    flight_path_angle: typing.ClassVar[AerodynamicsReferenceFrameAngles]
    heading_angle: typing.ClassVar[AerodynamicsReferenceFrameAngles]
    latitude_angle: typing.ClassVar[AerodynamicsReferenceFrameAngles]
    longitude_angle: typing.ClassVar[AerodynamicsReferenceFrameAngles]

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

class AerodynamicsReferenceFrames:
    """Enumeration of reference frame identifiers typical for aerodynamic calculations.
	
	Enumeration of reference frame identifiers typical for aerodynamic calculations. Note that the frames are also defined
	in the absence of any aerodynamic forces and/or atmosphere. They define frames of a body w.r.t. a central body, with
	the details given by Mooij (1994). The chain of frames starts from the inertial frame, to the frame fixed to the
	central body (corotating), to the vertical frame (defined by the body's relative position), the trajectory and aerodynamic frames
	(defined by the body's relative velocity) and finally the body's own body-fixed frame.
	
	
	:member inertial_frame:
	:member corotating_frame:
	:member vertical_frame:
	:member trajectory_frame:
	:member aerodynamic_frame:
	:member body_frame:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, AerodynamicsReferenceFrames]]
    aerodynamic_frame: typing.ClassVar[AerodynamicsReferenceFrames]
    body_frame: typing.ClassVar[AerodynamicsReferenceFrames]
    corotating_frame: typing.ClassVar[AerodynamicsReferenceFrames]
    inertial_frame: typing.ClassVar[AerodynamicsReferenceFrames]
    trajectory_frame: typing.ClassVar[AerodynamicsReferenceFrames]
    vertical_frame: typing.ClassVar[AerodynamicsReferenceFrames]

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

class AtmosphericFlightConditions(FlightConditions):
    """Object that calculates various state-derived quantities typically
	relevant for flight dynamics, for flight in an atmosphere.
	
	
		Object that calculates various state-derived quantities typically
		relevant for flight dynamics, for flight in an atmosphere, such
		as latitude,  longitude, altitude, density, Mach number etc. It
		also contains an ``AerodynamicAngleCalculator`` that computes
		derived angles (flight path, heading angle, etc.). This object is
		derived from ``FlightConditions``, which performs computations for
		non-atmospheric flight only. This object is stored inside a Body
		object, and represents the flight conditions of a single body
		w.r.t. a single central body.
	"""

    @property
    def aero_coefficient_independent_variables(self) -> list[float]:
        """
        List of current values of independent variables of aerodynamic
        coefficients. This list is only defined if the body has an
        :py:class:`~AerodynamicCoefficientInterface` that has
        dependencies on environmental variables (e.g. Mach number,
        angle of attack, etc.).
        """

    @property
    def aerodynamic_coefficient_interface(self) -> AerodynamicCoefficientInterface:
        """
        Object extracted from the same Body object as this
        :py:class:`~AtmosphericFlightConditions` object, which defines
        the aerodynamic coefficients.
        """

    @property
    def airspeed(self) -> float:
        """
        The airspeed of the body w.r.t. the atmosphere.
        """

    @property
    def airspeed_velocity(self) -> numpy.ndarray:
        """
        The velocity vector of the body w.r.t. the freestream
        atmosphere (e.g. vectorial counterpart of airspeed).
        """

    @property
    def control_surface_aero_coefficient_independent_variables(self) -> dict[str, list[float]]:
        """
        List of lists current values of independent variables of
        aerodynamic coefficients for control surfaces. The outer list
        defines the control surface, the inner list the values of the
        independent variables. This list is only defined if the body
        has an :py:class:`~AerodynamicCoefficientInterface` with
        control surfaces that have dependencies on environmental
        variables (e.g. Mach number, angle of attack, etc.).
        """

    @property
    def density(self) -> float:
        """
        The freestream atmospheric density at the body's current
        location.
        """

    @property
    def dynamic_pressure(self) -> float:
        """
        The freestream atmospheric dynamic pressure at the body's
        current location.
        """

    @property
    def mach_number(self) -> float:
        """
        The freestream Mach number of the body.
        """

    @property
    def pressure(self) -> float:
        """
        The freestream atmospheric static pressure at the body's
        current location.
        """

    @property
    def speed_of_sound(self) -> float:
        """
        The freestream atmospheric speed of sound at the body's current
        location.
        """

    @property
    def temperature(self) -> float:
        """
        The freestream atmospheric temperature at the body's current
        location.
        """

class Body:
    """Object that stores the environment properties and current state of
	a single body.
	
	
		Object that stores the environment properties and current state
		of a single celestial body (natural or artificial). Each separate
		environment model (gravity field, ephemeris, etc.) is stored as a
		member object in this class. During each time step, the Body gets
		updated to the current time/propagated state, and the current
		properties, in as much as they are time-dependent, can be
		extracted from this object
	"""
    ephemeris_frame_to_base_frame: ...
    inertia_tensor: numpy.ndarray
    rigid_body_properties: ...
    system_models: VehicleSystems

    def get_ground_station(self, station_name: str) -> GroundStation:
        ...

    def set_constant_mass(self, mass: float) -> None:
        ...

    def state_in_base_frame_from_ephemeris(self, time: float) -> numpy.ndarray:
        ...

    @property
    def aerodynamic_coefficient_interface(self) -> AerodynamicCoefficientInterface:
        """
        Object defining the aerodynamic coefficients of a vehicle (force-only, or force and moment)
        as a function of any number of independent variables. Depending on the selected type of model, the type of this attribute
        is of type AerodynamicCoefficientInterface, or a derived class thereof.
        """

    @aerodynamic_coefficient_interface.setter
    def aerodynamic_coefficient_interface(self, arg1: AerodynamicCoefficientInterface) -> None:
        ...

    @property
    def atmosphere_model(self) -> ...:
        """
        Atmosphere model of this body, used to calculate density, temperature, etc. at a given
        state/time. Depending on the selected type of model, the type of this attribute
        is of type AtmosphereModel, or a derived class thereof.
        """

    @atmosphere_model.setter
    def atmosphere_model(self, arg1: ...) -> None:
        ...

    @property
    def body_fixed_angular_velocity(self) -> numpy.ndarray:
        """
        Angular velocity vector of the body, expressed in body-fixed
        frame (see :py:attr:`~inertial_to_body_fixed_frame`).
        """

    @property
    def body_fixed_to_inertial_frame(self) -> numpy.ndarray:
        """
        The rotation from this Body's body-fixed frame to inertial
        frame (see :py:attr:`~inertial_to_body_fixed_frame`).
        """

    @property
    def body_fixed_to_inertial_frame_derivative(self) -> numpy.ndarray:
        """
        Time derivative of rotation matrix from this Body's body-fixed
        frame to inertial frame
        (see :py:attr:`~inertial_to_body_fixed_frame`).
        """

    @property
    def ephemeris(self) -> Ephemeris:
        """
        Ephemeris model of this body, used to calculate its current state as a function of time.
        Depending on the selected type of model, the type of this attribute
        is of type Ephemeris, or a derived class thereof.
        """

    @property
    def flight_conditions(self) -> FlightConditions:
        """
        Current flight conditions of the body, which can be accessed during the propagation to get the current altitude, aerodynamic angle calculator, longitude, etc.
        """

    @flight_conditions.setter
    def flight_conditions(self, arg1: FlightConditions) -> None:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
        Attribute of convenience, equivalent to ``.gravity_field_model.gravitational_parameter``
        """

    @property
    def gravity_field_model(self) -> GravityFieldModel:
        """
        Gravity field model of this body, used to define the exterior gravitational potential, and
        its gradient(s). Depending on the selected type of model, the type of this attribute
        is of type GravityFieldModel, or a derived class thereof.
        """

    @gravity_field_model.setter
    def gravity_field_model(self, arg1: GravityFieldModel) -> None:
        ...

    @property
    def ground_station_list(self) -> dict[str, GroundStation]:
        ...

    @property
    def inertial_angular_velocity(self) -> numpy.ndarray:
        """
        Angular velocity vector of the body, expressed in inertial
        frame (see :py:attr:`~inertial_to_body_fixed_frame`).
        """

    @property
    def inertial_to_body_fixed_frame(self) -> numpy.ndarray:
        """
        The rotation from inertial frame (with global frame
        orientation) to this Body's body-fixed frame. The rotation is
        always returned here as a rotation matrix.  If the body's
        rotational state is numerically propagated, this property gets
        extracted from the propagated state vector. If it is not
        propagated, the state is extracted from this body's rotational
        ephemeris.
        
        .. note:: This function is **only** valid during the
                  numerical propagation if any aspects of the dynamics
                  or dependent variables require the body's rotational
                  state.
        """

    @property
    def inertial_to_body_fixed_frame_derivative(self) -> numpy.ndarray:
        """
        Time derivative of rotation matrix from inertial frame to this
        Body's body-fixed frame
        (see :py:attr:`~inertial_to_body_fixed_frame`).
        """

    @property
    def mass(self) -> float:
        """
        Denotes the current mass of the vehicle, as used in the calculation of
        non-conservative acceleration. Note that defining a mass for a vehicle
        does *not* define a gravity field (this is done through a gravity field model).
        However, defining a gravity field model in a body automatically assigns it
        a mass here.
        Unlike the attributes containing the state, orientation, angular velocity
        of the Body, this attribute may be used to retrieve the state during the
        propagation *and* to define the mass of a vehicle
        """

    @mass.setter
    def mass(self, arg1: float) -> None:
        ...

    @property
    def position(self) -> numpy.ndarray:
        """
        The translational position of the Body, as set during the
        current step of the numerical propagation
        (see :py:attr:`~state`).
        """

    @property
    def rotation_model(self) -> RotationalEphemeris:
        """
        Object defining the orientation of a body, used to calculate the rotation to/from a body-fixed
        frame (and its derivate). Depending on the selected type of model, the type of this attribute
        is of type RotationalEphemeris, or a derived class thereof.
        """

    @rotation_model.setter
    def rotation_model(self, arg1: RotationalEphemeris) -> None:
        ...

    @property
    def shape_model(self) -> ShapeModel:
        """
        Shape model of this body, used to define the exterior shape of the body, for instance for
        the calculation of vehicle's altitude. Depending on the selected type of model, the type of this attribute
        is of type BodyShapeModel, or a derived class thereof.
        """

    @shape_model.setter
    def shape_model(self, arg1: ShapeModel) -> None:
        ...

    @property
    def state(self) -> numpy.ndarray:
        """
        The translational state of the Body, as set during the current
        step of the numerical propagation. The translational state
        stored here is always in Cartesian elements, w.r.t. the global
        frame origin, with axes along the global frame orientation. If
        the body's translational state is numerically propagated, this
        property gets extracted from the propagated state vector. If it
        is not propagated, the state is extracted from this body's
        ephemeris. In both cases, any required state transformations
        are automatically applied. Note that this function  is *only*
        valid during the numerical propagation if any aspects of the
        dynamics or dependent variables require the body's state.
        """

    @property
    def velocity(self) -> numpy.ndarray:
        """
        The translational velocity of the Body, as set during the
        current step of the numerical propagation
        (see :py:attr:`~state`).
        """

class ConstantEphemeris(Ephemeris):
    """
		"""

    @typing.overload
    def __init__(self, constant_state_function: typing.Callable[[], numpy.ndarray], reference_frame_origin: str='SSB', reference_frame_orientation: str='ECLIPJ2000') -> None:
        ...

    @typing.overload
    def __init__(self, constant_state: numpy.ndarray, reference_frame_origin: str='SSB', reference_frame_orientation: str='ECLIPJ2000') -> None:
        ...

    def update_constant_state(self, new_state: numpy.ndarray) -> None:
        ...

class ConstantFrequencyInterpolator(StationFrequencyInterpolator):

    def __init__(self, frequency: float) -> None:
        ...

class ControlSurfaceIncrementAerodynamicInterface:
    """<no_doc, only_dec>
	"""

class CustomBodyFixedDirectionCalculator(InertialBodyFixedDirectionCalculator):
    inertial_body_axis_direction_function: typing.Callable[[float], numpy.ndarray]

class CustomControlSurfaceIncrementAerodynamicInterface(ControlSurfaceIncrementAerodynamicInterface):
    """<no_doc, only_dec>
	"""

    def __init__(self, coefficient_function: typing.Callable[[list[float]], numpy.ndarray], independent_variable_names: list[AerodynamicCoefficientsIndependentVariables]) -> None:
        ...

class CustomInertialDirectionBasedRotationalEphemeris(RotationalEphemeris):

    @property
    def inertial_body_axis_calculator(self) -> ...:
        ...

class DirectLongitudeLibrationCalculator(LongitudeLibrationCalculator):

    def __init__(self, scaled_libration_amplitude: float) -> None:
        ...

class EngineModel:

    @property
    def thrust_magnitude_calculator(self) -> ...:
        ...

class Ephemeris:
    """
		"""

    def cartesian_position(self, current_time: float) -> numpy.ndarray:
        ...

    def cartesian_state(self, current_time: float) -> numpy.ndarray:
        ...

    def cartesian_velocity(self, current_time: float) -> numpy.ndarray:
        ...

    @property
    def frame_orientation(self) -> str:
        ...

    @property
    def frame_origin(self) -> str:
        ...

class FlightConditions:
    """Object that calculates various state-derived quantities typically
	relevant for flight dynamics.
	
	
		Object that calculates various state-derived quantities typically
		relevant for flight dynamics, such as latitude, longitude,
		altitude, etc. It also contains an
		:py:class:`~AerodynamicAngleCalculator` that computes derived
		angles (flight path, heading angle, etc.). This object is limited
		to non-atmospheric flight. For flight through Body objects
		endowed with an atmosphere model, the derived class
		:py:class:`~AtmosphericFlightConditions` is used. This object is
		stored inside a Body object, and represents the flight conditions
		of a single body w.r.t. a single central body.
	"""

    def update_conditions(self, current_time: float) -> None:
        ...

    @property
    def aerodynamic_angle_calculator(self) -> AerodynamicAngleCalculator:
        """
        The object that is responsible for computing various relevant
        flight dynamics angles and frame rotations.
        """

    @property
    def altitude(self) -> float:
        """
        The current time, at which this object was last updated
        """

    @property
    def body_centered_body_fixed_state(self) -> numpy.ndarray:
        """
        Cartesian translational state, expressed in a frame centered
        at, and fixed to, the central body. Note that, due to the
        rotation of the central body, the norm of the body-fixed,
        body-centered, velocity differs from the norm of the inertial
        body-centered velocity.
        """

    @property
    def geodetic_latitude(self) -> float:
        """
        The body-fixed geographic latitude of the body w.r.t. its
        central body.
        """

    @property
    def latitude(self) -> float:
        """
        The body-fixed geographic latitude of the body w.r.t. its
        central body.
        """

    @property
    def longitude(self) -> float:
        """
        The body-fixed longitude of the body w.r.t. its central body.
        """

    @property
    def time(self) -> float:
        """
        The current time, at which this object was last updated
        """

class GcrsToItrsRotationModel(RotationalEphemeris):
    pass

class GravityFieldModel:
    gravitational_parameter: float

    def __init__(self, gravitational_parameter: float, update_inertia_tensor: typing.Callable[[], None]=None) -> None:
        ...

    def get_gravitational_parameter(self) -> float:
        ...

class GroundStation:

    def set_pressure_function(self, pressure_function: typing.Callable[[float], float]) -> None:
        ...

    def set_relative_humidity_function(self, relative_humidity_function: typing.Callable[[float], float]) -> None:
        ...

    def set_temperature_function(self, temperature_function: typing.Callable[[float], float]) -> None:
        ...

    def set_transmitting_frequency_calculator(self, transmitting_frequency_calculator: ...) -> None:
        ...

    def set_water_vapor_partial_pressure_function(self, water_vapor_partial_pressure_function: typing.Callable[[float], float]) -> None:
        ...

    @property
    def pointing_angles_calculator(self) -> ...:
        ...

    @property
    def pressure_function(self) -> typing.Callable[[float], float]:
        ...

    @property
    def relative_humidity_function(self) -> typing.Callable[[float], float]:
        ...

    @property
    def station_state(self) -> GroundStationState:
        ...

    @property
    def temperature_function(self) -> typing.Callable[[float], float]:
        ...

class GroundStationState:

    def get_cartesian_position(self, seconds_since_epoch: float, target_frame_origin: str) -> numpy.ndarray:
        ...

    def get_cartesian_state(self, seconds_since_epoch: float, target_frame_origin: str) -> numpy.ndarray:
        ...

    @property
    def cartesian_positon_at_reference_epoch(self) -> numpy.ndarray:
        ...

    @property
    def geodetic_positon_at_reference_epoch(self) -> numpy.ndarray:
        ...

    @property
    def rotation_matrix_body_fixed_to_topocentric(self, arg1: float) -> numpy.ndarray:
        ...

    @property
    def spherical_positon_at_reference_epoch(self) -> numpy.ndarray:
        ...

class HypersonicLocalInclinationAnalysis(AerodynamicCoefficientGenerator36):

    def __init__(self, independent_variable_points: list[list[float]], body_shape: ..., number_of_lines: list[int], number_of_points: list[int], invert_orders: list[bool], selected_methods: list[list[int]], reference_area: float, reference_length: float, moment_reference_point: numpy.ndarray, save_pressure_coefficients: bool=False) -> None:
        """
        Class constructor, taking the shape of the vehicle, and various analysis options as input.
        
        	:param independent_variable_points:
        		List containing three lists, with each sublist containing the data points of each of the
        		independent variables for the coefficient generation. The physical meaning of each of the
        		three independent variables is: 0 = mach number, 1 = angle of attack, 2 = angle of sideslip.
        		Each of the subvectors must be sorted in ascending order.
        
        	:param body_shape:
        		Class that defines the shape of the vehicle as a continuous surface. The local inclination analysis
        		discretizes the surface of the vehicle into quadrilateral panels, defined by the other inputs to
        		this constructor. In case the :class:`tudat.geometry.SurfaceGeometry` object is made up of multiple
        		sub-shapes, different settings may be used for each
        
        	:param number_of_lines:
        		Number of discretization points in the first independent surface variable of each of the subparts of body_shape.
        		The size of this list should match the number of parts of which the body_shape is composed. The first independent
        		variable of a subpart typically runs along the longitudinal vehicle direction
        
        	:param number_of_points:
        		Number of discretization points in the second independent surface variable of each of the subparts of body_shape.
        		The size of this list should match the number of parts of which the body_shape is composed. The first independent
        		variable of a subpart typically runs along the lateral vehicle direction
        
        	:param invert_orders:
        		Booleans to denote whether the surface normals of the panels of each discretized body_shape subpart are to be inverted
        		(i.e. inward-facing->outward facing or vice versa). The size of this list should match the number of parts of which the body_shape is composed.
        
        	:param selected_methods:
        		Double list of selected local inclination methods, the first index (outer list) represents compression or expansion (0 and 1),
        		the second index (inner list) denotes the vehicle part index. The size of this inner list should match the number of parts of which the body_shape is composed.
        		The int defining the method type is interpreted as follows.
        		For the compression methods, the following are available:
        		*  0: Newtonian Method.
        		*  1: Modified Newtonian.
        		*  2 and 3: not available at this moment.
        		*  4: Tangent-wedge method.
        		*  5: Tangent-cone method.
        		*  6: Modified Dahlem-Buck method.
        		*  7: VanDyke unified pressure method.
        		*  8: Smyth Delta Wing method.
        		*  9: Hankey flat surface method
        		The expansion method has the following options:
        		*  0: Vacuum Pressure coefficient method.
        		*  1: Zero Pressure function.
        		*  4: High Mach base pressure method.
        		*  3 or 5: Prandtl-Meyer method.
        		*  6: ACM empirical pressure coefficient.
        
        	:param reference_area:
        		Reference area used to non-dimensionalize aerodynamic forces and moments.
        
        	:param moment_reference_point:
        		Reference point wrt which aerodynamic moments are calculated.
        
        	:param save_pressure_coefficients:
        		Boolean denoting whether to save the pressure coefficients that are computed to files
        """

    def clear_data(self) -> None:
        ...

class InertialBodyFixedDirectionCalculator:
    pass

class KeplerEphemeris(Ephemeris):
    pass

class LongitudeLibrationCalculator:
    pass

class MultiArcEphemeris(Ephemeris):

    def __init__(self, single_arc_ephemerides: dict[float, Ephemeris], reference_frame_origin: str='SSB', reference_frame_orientation: str='ECLIPJ2000') -> None:
        ...

class PointingAnglesCalculator:

    def calculate_azimuth_angle(self, inertial_vector_to_target: numpy.ndarray, time: float) -> float:
        ...

    def calculate_elevation_angle(self, inertial_vector_to_target: numpy.ndarray, time: float) -> float:
        ...

    def convert_inertial_vector_to_topocentric(self, inertial_vector: numpy.ndarray, time: float) -> numpy.ndarray:
        ...

class PolyhedronGravityField(GravityFieldModel):

    @property
    def vertices_coordinates(self) -> numpy.ndarray:
        ...

    @property
    def vertices_defining_each_facet(self) -> numpy.ndarray:
        ...

    @property
    def volume(self) -> float:
        ...

class RotationalEphemeris:
    """
		"""

    def angular_velocity_in_body_fixed_frame(self, time: float) -> numpy.ndarray:
        ...

    def angular_velocity_in_inertial_frame(self, time: float) -> numpy.ndarray:
        ...

    def body_fixed_to_inertial_rotation(self, time: float) -> numpy.ndarray:
        ...

    def inertial_to_body_fixed_rotation(self, time: float) -> numpy.ndarray:
        ...

    def time_derivative_body_fixed_to_inertial_rotation(self, time: float) -> numpy.ndarray:
        ...

    def time_derivative_inertial_to_body_fixed_rotation(self, time: float) -> numpy.ndarray:
        ...

    @property
    def body_fixed_frame_name(self) -> str:
        ...

    @property
    def inertial_frame_name(self) -> str:
        ...

class ShapeModel:

    def get_average_radius(self) -> float:
        ...

    @property
    def average_radius(self) -> float:
        ...

class SphericalHarmonicsGravityField(GravityFieldModel):
    cosine_coefficients: numpy.ndarray
    sine_coefficients: numpy.ndarray

    @property
    def maximum_degree(self) -> float:
        ...

    @property
    def maximum_order(self) -> float:
        ...

    @property
    def reference_radius(self) -> float:
        ...

class StationFrequencyInterpolator:
    """
		"""

class SynchronousRotationalEphemeris(RotationalEphemeris):
    libration_calculator: LongitudeLibrationCalculator

class SystemOfBodies:
    """Object that contains a set of Body objects and associated frame
	information.
	
	
		Object that contains a set of Body objects and associated frame
		information. This object stored the entire environment for a
		typical Tudat numerical simulation, and is fundamental for the
		overall Tudat architecture.
	"""

    def add_body(self, body_to_add: Body, body_name: str, process_body: bool=1) -> None:
        """
        This function adds an existing body, which the user has
        separately created, to the :py:class:`~SystemOfBodies`.
        
        
        	:param body_to_add:
        		Body object that is to be added.
        
        	:param body_name:
        		Name of the Body that is to be added.
        
        	:param process_body:
        		Variable that defines whether this new Body will have its
        		global frame origin/orientation set to conform to rest of
        		the environment.
        
        		.. warning:: Only in very rare cases should this variable be
        		             anything other than ``True``. Users are
        		             recommended to keep this default value intact.
        """

    def create_empty_body(self, body_name: str, process_body: bool=1) -> None:
        """
        This function creates a new empty body.
        
        	This function creates a new empty body, and adds it to the
        	:py:class:`~SystemOfBodies`. Since the body is empty, it will
        	not have any environment models defined. These must all be
        	added manually by a user.
        
        
        	:param body_name:
        		Name of the Body that is to be added
        
        	:param process_body:
        		Variable that defines whether this new Body will have its
        		global frame origin/orientation set to conform to rest of
        		the environment.
        
        		.. warning:: Only in very rare cases should
        		             this variable be anything other than ``True``.
        		             Users are recommended to keep this default value
        		             intact.
        """

    def does_body_exist(self, body_name: str) -> bool:
        ...

    def get(self, body_name: str) -> Body:
        """
        This function extracts a single Body object from the SystemOfBodies.
        
        	:param body_name:
        		Name of the Body that is to be retrieved.
        
        	:return:
        		Body object of the requested name
        """

    def get_body(self, body_name: str) -> Body:
        """
        Deprecated version of :py:func:`~get`
        """

    def global_frame_orientation(self) -> str:
        ...

    def global_frame_origin(self) -> str:
        ...

    def list_of_bodies(self) -> list[str]:
        ...

    def remove_body(self, body_name: str) -> None:
        """
        This function removes an existing body from the
        :py:class:`~SystemOfBodies`.
        
        
        
        	.. warning:: This function does *not* necessarily delete the
        	             Body object, it only removes it from this object.
        	             If any existing models in the simulation refer to
        	             this Body, it will persist in memory.
        
        
        	:param body_name:
        		Name of the Body that is to be removed.
        """

class TabulatedEphemeris(Ephemeris):

    @property
    def interpolator(self) -> ...:
        ...

    @interpolator.setter
    def interpolator(*args, **kwargs):
        """
        (arg0: tudatpy.numerical_simulation.environment.expose_environment.TabulatedEphemeris, arg1: tudat::interpolators::OneDimensionalInterpolator<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>>) -> None
        """

class Tle:

    @typing.overload
    def __init__(self, lines: str) -> None:
        ...

    @typing.overload
    def __init__(self, line_1: str, line_2: str) -> None:
        ...

    def get_arg_of_perigee(self) -> float:
        ...

    def get_b_star(self) -> float:
        ...

    def get_eccentricity(self) -> float:
        ...

    @typing.overload
    def get_epoch(self) -> float:
        ...

    @typing.overload
    def get_epoch(self) -> float:
        ...

    def get_inclination(self) -> float:
        ...

    def get_mean_anomaly(self) -> float:
        ...

    def get_mean_motion(self) -> float:
        ...

    def get_right_ascension(self) -> float:
        ...

class TleEphemeris(Ephemeris):

    def __init__(self, frame_origin: str='Earth', frame_orientation: str='J2000', tle: Tle=None, use_sdp: bool=False) -> None:
        ...

class VehicleSystems:
    """
		"""

    def __init__(self) -> None:
        ...

    def get_control_surface_deflection(self, control_surface_id: str) -> float:
        ...

    def get_engine_model(self, engine_name: str) -> ...:
        ...

    def set_control_surface_deflection(self, control_surface_id: str, deflection_angle: float) -> None:
        ...

    def set_transponder_turnaround_ratio(self, transponder_ratio_per_uplink_and_downlink_frequency_band: dict[tuple[..., ...], float]) -> None:
        ...

def get_default_local_inclination_angle_of_attack_points() -> list[float]:
    ...

def get_default_local_inclination_mach_points(mach_regime: str='Full') -> list[float]:
    ...

def get_default_local_inclination_sideslip_angle_points() -> list[float]:
    ...

def get_local_inclination_mesh(local_inclination_analysis_object: HypersonicLocalInclinationAnalysis) -> tuple[list[numpy.ndarray], list[numpy.ndarray]]:
    ...

def get_local_inclination_total_vehicle_area(local_inclination_analysis_object: HypersonicLocalInclinationAnalysis) -> float:
    ...

def save_vehicle_mesh_to_file(local_inclination_analysis_object: HypersonicLocalInclinationAnalysis, output_directory: str, output_file_prefix: str='') -> None:
    """Function to save the mesh used for a hypersonic local inclination analysis to a file.
	
	Function to save the mesh used for a hypersonic local inclination analysis to a file. This function saves
	two files to the specified directory, with filenames: "ShapeFile.dat" and "SurfaceNormalFile.dat", where these
	files names may be prefixed by an optional string (see below). The first of these files contains four columns defining
	the surface points that define mesh, with Column 0: point index; Column 1: x-position of point; Column 1: y-position of point;
	Column 2: z-position of point. The second file contains four columns with Column 0: point index; Column 1: x-component of surface normal;
	Column 1: y-position of surface normal; Column 2: z-position of surface normal.
	
	
	:param local_inclination_analysis_object:
			Object used to calculate the aerodynamics of the vehicle
	
	:param output_directory:
			Directory to which the files are to be saved
	
	:param output_file_prefix:
			Optional prefix of output file names
	"""

def transform_to_inertial_orientation(state_in_body_fixed_frame: numpy.ndarray, current_time: float, rotational_ephemeris: RotationalEphemeris) -> numpy.ndarray:
    ...
aerodynamic_frame: AerodynamicsReferenceFrames
altitude_dependent: AerodynamicCoefficientsIndependentVariables
angle_of_attack: AerodynamicsReferenceFrameAngles
angle_of_attack_dependent: AerodynamicCoefficientsIndependentVariables
angle_of_sideslip: AerodynamicsReferenceFrameAngles
anomalous_o_number_density_dependent: AerodynamicCoefficientsIndependentVariables
ar_number_density_dependent: AerodynamicCoefficientsIndependentVariables
bank_angle: AerodynamicsReferenceFrameAngles
body_frame: AerodynamicsReferenceFrames
control_surface_deflection_dependent: AerodynamicCoefficientsIndependentVariables
corotating_frame: AerodynamicsReferenceFrames
flight_path_angle: AerodynamicsReferenceFrameAngles
h_number_density_dependent: AerodynamicCoefficientsIndependentVariables
he_number_density_dependent: AerodynamicCoefficientsIndependentVariables
heading_angle: AerodynamicsReferenceFrameAngles
inertial_frame: AerodynamicsReferenceFrames
latitude_angle: AerodynamicsReferenceFrameAngles
longitude_angle: AerodynamicsReferenceFrameAngles
mach_number_dependent: AerodynamicCoefficientsIndependentVariables
n2_number_density_dependent: AerodynamicCoefficientsIndependentVariables
n_number_density_dependent: AerodynamicCoefficientsIndependentVariables
negative_aerodynamic_frame_coefficients: AerodynamicCoefficientFrames
negative_body_fixed_frame_coefficients: AerodynamicCoefficientFrames
o2_number_density_dependent: AerodynamicCoefficientsIndependentVariables
o_number_density_dependent: AerodynamicCoefficientsIndependentVariables
positive_aerodynamic_frame_coefficients: AerodynamicCoefficientFrames
positive_body_fixed_frame_coefficients: AerodynamicCoefficientFrames
sideslip_angle_dependent: AerodynamicCoefficientsIndependentVariables
temperature_dependent: AerodynamicCoefficientsIndependentVariables
time_dependent: AerodynamicCoefficientsIndependentVariables
trajectory_frame: AerodynamicsReferenceFrames
undefined_independent_variable: AerodynamicCoefficientsIndependentVariables
velocity_dependent: AerodynamicCoefficientsIndependentVariables
vertical_frame: AerodynamicsReferenceFrames