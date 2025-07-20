import numpy
import pybind11_stubgen.typing_ext
from ...astro import time_representation
from ...dynamics.environment_setup import aerodynamic_coefficients
from ...math import geometry
from ...math import interpolators
import typing
__all__ = ['AerodynamicAngleCalculator', 'AerodynamicAngleRotationalEphemeris', 'AerodynamicCoefficientGenerator36', 'AerodynamicCoefficientInterface', 'AtmosphereModel', 'AtmosphericFlightConditions', 'Body', 'BodyShapeModel', 'CannonballRadiationPressureTargetModel', 'ConstantEphemeris', 'ConstantTransmittingFrequencyCalculator', 'ContinuousInterpolatedMeteoData', 'ControlSurfaceIncrementAerodynamicInterface', 'CustomBodyFixedDirectionCalculator', 'CustomControlSurfaceIncrementAerodynamicInterface', 'CustomInertialDirectionBasedRotationalEphemeris', 'DirectLongitudeLibrationCalculator', 'EarthOrientationAnglesCalculator', 'EngineModel', 'Ephemeris', 'FlightConditions', 'GcrsToItrsRotationModel', 'GravityFieldModel', 'GravityFieldVariationModel', 'GroundStation', 'GroundStationState', 'HypersonicLocalInclinationAnalysis', 'InertialBodyFixedDirectionCalculator', 'IonosphereModel', 'KeplerEphemeris', 'LongitudeLibrationCalculator', 'MeteoDataEntries', 'MultiArcEphemeris', 'PiecewiseLinearFrequencyInterpolator', 'PointingAnglesCalculator', 'PolyhedronGravityField', 'RadiationPressureTargetModel', 'RadiationSourceModel', 'RigidBodyProperties', 'RotationalEphemeris', 'SphericalHarmonicsGravityField', 'StationMeteoData', 'SynchronousRotationalEphemeris', 'SystemOfBodies', 'TabulatedEphemeris', 'TimeDependentSphericalHarmonicsGravityField', 'TimingSystem', 'Tle', 'TleEphemeris', 'TransmittingFrequencyCalculator', 'VehicleSystems', 'dew_point_meteo_data', 'get_default_local_inclination_angle_of_attack_points', 'get_default_local_inclination_mach_points', 'get_default_local_inclination_sideslip_angle_points', 'get_local_inclination_mesh', 'get_local_inclination_total_vehicle_area', 'pressure_meteo_data', 'relative_humidity_meteo_data', 'save_vehicle_mesh_to_file', 'temperature_meteo_data', 'transform_to_inertial_orientation', 'water_vapor_pressure_meteo_data']

class AerodynamicAngleCalculator:
    """Object to calculate (aerodynamic) orientation angles, and frame transformations,
    from current vehicle state.
    
    
    Object to calculate (aerodynamic) orientation angles (list given by the :class:`~AerodynamicsReferenceFrameAngles` enum)
    and transformations between frames (list given by the :class:`~AerodynamicsReferenceFrames` enum) from current vehicle state."""

    def get_angle(self, angle_type: aerodynamic_coefficients.AerodynamicsReferenceFrameAngles) -> float:
        """
                 Function to get a single orientation angle
        
        
                 Function to get a single orientation angle. This function
                 is meant to be used only *during* a numerical propagation, in particular
                 for the definition of a custom (e.g. guidance) model.
        
        
                 Parameters
                 ----------
                 original_frame : AerodynamicsReferenceFrameAngles
                     The identifier for the angle that is to be returned
        
                 Returns
                 -------
                 double
                     Value of requested angle
        """

    def get_rotation_matrix_between_frames(self, original_frame: aerodynamic_coefficients.AerodynamicsReferenceFrames, target_frame: aerodynamic_coefficients.AerodynamicsReferenceFrames) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 Function to get the rotation matrix between two frames.
        
        
                 Function to get the rotation matrix between two frames. This function
                 is meant to be used only *during* a numerical propagation, in particular
                 for the definition of a custom (e.g. guidance) model.
        
        
                 Parameters
                 ----------
                 original_frame : AerodynamicsReferenceFrames
                     The frame :math:`A` from which the rotation matrix is to be calculated
        
                 target_frame : AerodynamicsReferenceFrames
                     The frame :math:`B` to which the rotation matrix is to be calculated
        
                 Returns
                 -------
                 numpy.ndarray
                     Rotation matrix :math:`\\mathbf{R}^{B/A}` from frame :math:`A` to frame `B`
        """

    def set_body_orientation_angle_functions(self, angle_of_attack_function: typing.Callable[[], float]=None, angle_of_sideslip_function: typing.Callable[[], float]=None, bank_angle_function: typing.Callable[[], float]=None, angle_update_function: typing.Callable[[float], None]=None, silence_warnings: bool=False) -> None:
        ...

    def set_body_orientation_angles(self, angle_of_attack: float=..., angle_of_sideslip: float=..., bank_angle: float=..., silence_warnings: bool=False) -> None:
        ...

class AerodynamicAngleRotationalEphemeris(RotationalEphemeris):

    def reset_aerodynamic_angle_function(self, arg0: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]) -> None:
        ...

class AerodynamicCoefficientGenerator36(AerodynamicCoefficientInterface):
    """<no_doc, only_dec>"""

class AerodynamicCoefficientInterface:
    """Base class for computing the current aerodynamic coefficients of the body
    
    
    Base class for computing the current aerodynamic coefficients of the body. The implementation of the computation
    depends on the choice of aerodynamic coefficient model (see :ref:`aerodynamic_coefficients` for available options).
    During the propagation, this object is automatically updated to the current state by the :class:`~AtmosphericFlightConditions` object.
    The user may override the current aerodynamic coefficients when using, for instance, a custom aerodynamic guidance model
    (see `here <https://docs.tudat.space/en/latest/_src_getting_started/_src_examples/notebooks/propagation/reentry_trajectory.html>`_ for an example).
    using the member functions of this class."""

    def current_control_surface_force_coefficient_increment(self, control_surface_name: str) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 Function to get the contribution from a single control surface to the aerodynamic force coefficient, as compute by last call to :meth:`~update_full_coefficients`
        
        
        
                 Parameters
                 ----------
                 control_surface_name : str
                     The name of the control surface for which the contribution is to be retrieved
        
                 Returns
                 -------
                 numpy.ndarray
                     Contribution from the requested control surface to the aerodynamic force coefficient
        """

    def current_control_surface_moment_coefficient_increment(self, control_surface_name: str) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 Function to get the contribution from a single control surface to the aerodynamic moment coefficients, as compute by last call to :meth:`~update_full_coefficients`
        
        
        
                 Parameters
                 ----------
                 control_surface_name : str
                     The name of the control surface for which the contribution is to be retrieved
        
                 Returns
                 -------
                 numpy.ndarray
                     Contribution from the requested control surface to the aerodynamic moment coefficients
        """

    def set_control_surface_increments(self, control_surface_list: dict[str, ...]) -> None:
        """
        No documentation found.
        """

    def update_coefficients(self, independent_variables: list[float], time: float) -> None:
        """
                 Function to update the aerodynamic coefficients of the body only
        
        
                 Function to update the aerodynamic coefficients of the body only (without the control surface contribution),
                 based on the current state. This function may be called by the user, but will set *only* the
                 :attr:`~current_force_coefficients` and :attr:`~current_moment_coefficients` (while leaving the
                 :attr:`~current_control_surface_free_force_coefficients` and :attr:`~current_control_surface_free_moment_coefficients` unchanged)
        
        
                 Parameters
                 ----------
                 independent_variables : list[float]
                     List of inputs from which the aerodynamic coefficients are to be computed, with each entry corresponding to the
                     value of the physical variable defined by the :attr:`independent_variable_names` attribute.
        
                 time : float
                     Current time (in seconds since J2000)
        
                 Returns
                 -------
                 numpy.ndarray
                     Contribution from the requested control surface to the aerodynamic moment coefficients
        """

    def update_full_coefficients(self, independent_variables: list[float], control_surface_independent_variables: dict[str, list[float]], time: float, check_force_contribution: bool=True) -> None:
        """
                 Function to update the aerodynamic coefficients, from both the body and its control surfaces
        
        
                 Function to update the aerodynamic coefficients of both the body and its control surfaces,
                 based on the current state. This function will call the :meth:`~update_coefficients` function to update the body coefficients.
                 This function may be called by the user, and will set the following attributes:
                 :attr:`~current_force_coefficients`, :attr:`~current_moment_coefficients` ,
                 :attr:`~current_control_surface_free_force_coefficients` and :attr:`~current_control_surface_free_moment_coefficients`.
                 In addition, it will modify the coefficients returned by the :meth:`~current_control_surface_force_coefficient_increment` and
                 :meth:`~current_control_surface_moment_coefficient_increment` functions
        
        
                 Parameters
                 ----------
                 independent_variables : list[float]
                     List of inputs from which the aerodynamic coefficients of the body are to be computed, with each entry corresponding to the
                     value of the physical variable defined by the :attr:`independent_variable_names` attribute.
        
                 control_surface_independent_variables : dict[str,list[float]]
                     List of inputs from which the control surface aerodynamic coefficients are to be computed (with dictionary key the control surface name),
                     with each entry corresponding to the
                     value of the physical variable defined by the :attr:`control_surface_independent_variable_names` attribute.
        
                 time : float
                     Current time (in seconds since J2000)
        
                 check_force_contribution : bool, default = True
                     Boolean that determines if the force contribution to the aerodynamic moments should be added. Note that this input is
                     only used if the :attr:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings.add_force_contribution_to_moments` attribute is set to True.
        """

    @property
    def control_surface_independent_variable_names(self) -> dict[str, list[aerodynamic_coefficients.AerodynamicCoefficientsIndependentVariables]]:
        """
                 **read-only**
        
                 List of independent variables from which the aerodynamic coefficients of each control surface are computed, with dictionary key being the control surface name (e.g. required input to :meth:`~update_full_coefficients` function).
        
        
                 :type: dict[str,list[AerodynamicCoefficientsIndependentVariables]]
        """

    @property
    def current_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
                 **read-only**
        
                 Concatenation of :attr:`~current_force_coefficients` and :attr:`~current_moment_coefficients`
        
        
                 :type: numpy.ndarray
        """

    @property
    def current_control_surface_free_force_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Same as :attr:`current_force_coefficients`, but without contribution (if any) from control surfaces
        
        
                 :type: numpy.ndarray
        """

    @property
    def current_control_surface_free_moment_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Same as :attr:`current_moment_coefficients`, but without contribution (if any) from control surfaces
        
        
                 :type: numpy.ndarray
        """

    @property
    def current_force_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 The current aerodynamic force coefficients, in the frame defined by the :attr:`~force_coefficient_frame` attribute,
                 as computed by the last call to the :meth:`~update_coefficients` function.
        
        
                 :type: numpy.ndarray
        """

    @property
    def current_moment_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 The current aerodynamic moment coefficients, in the frame defined by the :attr:`~moment_coefficient_frame` attribute,
                 as computed by the last call to the :meth:`~update_coefficients` function.
        
        
                 :type: numpy.ndarray
        """

    @property
    def force_coefficient_frame(self) -> aerodynamic_coefficients.AerodynamicCoefficientFrames:
        """
                 **read-only**
        
                 Reference frame in which the  :attr:`~current_force_coefficients` are defined
        
        
                 :type: AerodynamicCoefficientFrames
        """

    @property
    def independent_variable_names(self) -> list[aerodynamic_coefficients.AerodynamicCoefficientsIndependentVariables]:
        """
                 **read-only**
        
                 List of independent variables from which the aerodynamic coefficients are computed (e.g. required input to :meth:`~update_coefficients` function).
        
        
                 :type: list[AerodynamicCoefficientsIndependentVariables]
        """

    @property
    def moment_coefficient_frame(self) -> aerodynamic_coefficients.AerodynamicCoefficientFrames:
        """
                 **read-only**
        
                 Reference frame in which the  :attr:`~current_moment_coefficients` are defined
        
        
                 :type: AerodynamicCoefficientFrames
        """

    @property
    def reference_area(self) -> float:
        """
                 **read-only**
        
                 The aerodynamic reference area :math:`A` of the coefficients
        
        
                 :type: float
        """

class AtmosphereModel:
    """Object that provides the atmospheric properties of the body.
    
    Object that provides the atmospheric properties of the body, as a function of altitude, latitude, longitude and time. Depending on the implementation, the
    this dependence may be limited to altitude-only (e.g. standard atmosphere models). This object can be accessed directly by the user to compute atmospheric properties
    outside the loop of the propagation by calling one of its member functions. During the propagation, each body undergoing aerodynamic forces has a :class:`AtmosphericFlightConditions`
    object associated with it (accessed from a boyd through :attr:`~Body.flight_condition`) that links the atmosphere model to the aerodynamic model."""

    def get_density(self, altitude: float, longitude: float, latitude: float, time: float) -> float:
        """
                 Function to compute the atmospheric freestream density at a given location.
        
                 Parameters
                 ----------
                 altitude : float
                     Local altitude above the body surface at which the property is to be computed
                 latitude : float
                     Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 longitude : float
                     Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 time : float
                     Time (in seconds since J2000 TDB) at which the property is to be computed.
        
                 Returns
                 -------
                 float
                     Freestream density at the given time and location
        """

    def get_number_density(self, species: aerodynamic_coefficients.AtmosphericCompositionSpecies, altitude: float, longitude: float, latitude: float, time: float) -> float:
        """
                 Function to compute the atmospheric freestream number density of a given specie at a given location.
        
                 Parameters
                 ----------
                 species : AtmosphericCompositionSpecies
                     Atmospheric species for which the number density is to be computed
                 altitude : float
                     Local altitude above the body surface at which the property is to be computed
                 latitude : float
                     Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 longitude : float
                     Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 time : float
                     Time (in seconds since J2000 TDB) at which the property is to be computed.
        
                 Returns
                 -------
                 float
                     Freestream number density of the requested specie at the given time and location
        """

    def get_pressure(self, altitude: float, longitude: float, latitude: float, time: float) -> float:
        """
                 Function to compute the atmospheric freestream static pressure at a given location.
        
                 Parameters
                 ----------
                 altitude : float
                     Local altitude above the body surface at which the property is to be computed
                 latitude : float
                     Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 longitude : float
                     Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 time : float
                     Time (in seconds since J2000 TDB) at which the property is to be computed.
        
                 Returns
                 -------
                 float
                     Freestream static pressure at the given time and location
        """

    def get_speed_of_sound(self, altitude: float, longitude: float, latitude: float, time: float) -> float:
        """
                 Function to compute the atmospheric freestream speed of sound at a given location.
        
                 Parameters
                 ----------
                 altitude : float
                     Local altitude above the body surface at which the property is to be computed
                 latitude : float
                     Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 longitude : float
                     Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 time : float
                     Time (in seconds since J2000 TDB) at which the property is to be computed.
        
                 Returns
                 -------
                 float
                     Freestream speed of sound at the given time and location
        """

    def get_temperature(self, altitude: float, longitude: float, latitude: float, time: float) -> float:
        """
                 Function to compute the atmospheric freestream temperature at a given location.
        
                 Parameters
                 ----------
                 altitude : float
                     Local altitude above the body surface at which the property is to be computed
                 latitude : float
                     Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 longitude : float
                     Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
                 time : float
                     Time (in seconds since J2000 TDB) at which the property is to be computed.
        
                 Returns
                 -------
                 float
                     Freestream temperature at the given time and location
        """

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
    w.r.t. a single central body."""

    @property
    def aero_coefficient_independent_variables(self) -> list[float]:
        """
                 **read-only**
        
                 List of current values of independent variables of aerodynamic
                 coefficients. This list is only defined if the body has an
                 :py:class:`~AerodynamicCoefficientInterface` that has
                 dependencies on environmental variables (e.g. Mach number,
                 angle of attack, etc.).
        
        
                 :type: numpy.ndarray
        """

    @property
    def aerodynamic_coefficient_interface(self) -> AerodynamicCoefficientInterface:
        """
                 **read-only**
        
                 Object extracted from the same Body object as this
                 :py:class:`~AtmosphericFlightConditions` object, which defines
                 the aerodynamic coefficients.
        
        
                 :type: AerodynamicCoefficientInterface
        """

    @property
    def airspeed(self) -> float:
        """
                 **read-only**
        
                 The airspeed of the body w.r.t. the atmosphere.
        
        
                 :type: float
        """

    @property
    def airspeed_velocity(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 The velocity vector of the body w.r.t. the freestream
                 atmosphere (e.g. vectorial counterpart of airspeed).
        
        
                 :type: numpy.ndarray
        """

    @property
    def control_surface_aero_coefficient_independent_variables(self) -> dict[str, list[float]]:
        """
                 **read-only**
        
                 List of lists current values of independent variables of
                 aerodynamic coefficients for control surfaces. The outer list
                 defines the control surface, the inner list the values of the
                 independent variables. This list is only defined if the body
                 has an :py:class:`~AerodynamicCoefficientInterface` with
                 control surfaces that have dependencies on environmental
                 variables (e.g. Mach number, angle of attack, etc.).
        
        
                 :type: numpy.ndarray
        """

    @property
    def density(self) -> float:
        """
                 **read-only**
        
                 The freestream atmospheric density at the body's current
                 location.
        
        
                 :type: float
        """

    @property
    def dynamic_pressure(self) -> float:
        """
                 **read-only**
        
                 The freestream atmospheric dynamic pressure at the body's
                 current location.
        
        
                 :type: float
        """

    @property
    def mach_number(self) -> float:
        """
                 **read-only**
        
                 The freestream Mach number of the body.
        
        
                 :type: float
        """

    @property
    def pressure(self) -> float:
        """
                 **read-only**
        
                 The freestream atmospheric static pressure at the body's
                 current location.
        
        
                 :type: float
        """

    @property
    def speed_of_sound(self) -> float:
        """
                 **read-only**
        
                 The freestream atmospheric speed of sound at the body's current
                 location.
        
        
                 :type: float
        """

    @property
    def temperature(self) -> float:
        """
                 **read-only**
        
                 The freestream atmospheric temperature at the body's current
                 location.
        
        
                 :type: float
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
    extracted from this object"""

    def get_ground_station(self, station_name: str) -> GroundStation:
        """
                 This function extracts a ground station object from the body.
        
                 This function extracts a ground station object, for a station of a given name, from the body.
                 If no station of this name exists, an exception is thrown
        
        
                 Parameters
                 ----------
                 station_name : str
                     Name of the ground station that is to be retrieved.
        
                 Returns
                 -------
                 GroundStation
                     Ground station object of the station of requested name
        """

    def get_ionosphere_model(self) -> IonosphereModel:
        ...

    def state_in_base_frame_from_ephemeris(self, time: time_representation.Time) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
                 This function returns the body's state, as computed from its ephemeris model (extracted from :attr:`~Body.ephemeris`) at the current time, and (if needed)
                 translates this state to the global frame origin. For the case where the origin of the body's ephemeris (extracted from :attr:`~Ephemeris.frame_origin`) is equal to the
                 global frame origin of the system of bodies it is in (extracted from :attr:`SystemOfBodies.global_frame_origin`), this function is equal to ``Body.ephemeris.cartesian_state( time )``.
                 Where the global frame origin and ephemeris origin is not equal, other bodies' ephemerides are queried as needed to provide this body's state w.r.t. the global frame origin
        
        
                 Parameters
                 ----------
                 time : float
                     Time (in TDB seconds since J2000) at which the state is to be computed
                 Returns
                 -------
                 numpy.ndarray
                     Cartesian state (position and velocity) of the body w.r.t. the global frame origin at the requested time.
        """

    @property
    def aerodynamic_coefficient_interface(self) -> AerodynamicCoefficientInterface:
        """
                 Object defining the aerodynamic coefficients of the body (force-only, or force and moment)
                 as a function of any number of independent variables. Depending on the selected type of model, the type of this attribute
                 is of type AerodynamicCoefficientInterface, or a derived class thereof.
        
        
                 :type: AerodynamicCoefficientInterface
        """

    @aerodynamic_coefficient_interface.setter
    def aerodynamic_coefficient_interface(self, arg1: AerodynamicCoefficientInterface) -> None:
        ...

    @property
    def atmosphere_model(self) -> AtmosphereModel:
        """
                 Object defining the atmosphere model of this body, used to calculate density, temperature, etc. at a given
                 state/time. Depending on the selected type of model, the type of this attribute
                 is of type :class:`~AtmosphereModel`, or a derived class thereof.
        
        
                 :type: AtmosphereModel
        """

    @atmosphere_model.setter
    def atmosphere_model(self, arg1: AtmosphereModel) -> None:
        ...

    @property
    def body_fixed_angular_velocity(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Angular velocity vector of the body, expressed in body-fixed
                 frame (see :py:attr:`~inertial_to_body_fixed_frame`).
        
        
                 :type: numpy.ndarray
        """

    @property
    def body_fixed_to_inertial_frame(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 **read-only**
        
                 The rotation from this Body's body-fixed frame to inertial
                 frame (see :py:attr:`~inertial_to_body_fixed_frame`).
        
        
                 :type: numpy.ndarray
        """

    @property
    def body_fixed_to_inertial_frame_derivative(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 **read-only**
        
                 Time derivative of rotation matrix from this Body's body-fixed
                 frame to inertial frame
                 (see :py:attr:`~inertial_to_body_fixed_frame`).
        
        
                 :type: numpy.ndarray
        """

    @property
    def ephemeris(self) -> Ephemeris:
        """
                 Object defining the ephemeris model of this body, used to calculate its current state as a function of time.
                 Depending on the selected type of model, the type of this attribute
                 is of type :class:`~Ephemeris`, or a derived class thereof.
        
        
                 :type: Ephemeris
        """

    @ephemeris.setter
    def ephemeris(self, arg1: Ephemeris) -> None:
        ...

    @property
    def flight_conditions(self) -> FlightConditions:
        """
                 Object used to calculated and store the current flight conditions of a vehicle (altitude, latitude, longitude,
                 flight-path angle, etc.) w.r.t. a central body. In case the central body contains an atmosphere, this object
                 also stores current local density, Mach number, etc. This object is typically used for aerodynamic accelerations,
                 guidance models or other central-body-related custom models.
        
        
                 :type: FlightConditions
        """

    @flight_conditions.setter
    def flight_conditions(self, arg1: FlightConditions) -> None:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
                 **read-only**
        
                 Attribute of convenience, equivalent to ``.gravity_field_model.gravitational_parameter``
        
        
                 :type: float
        """

    @property
    def gravity_field_model(self) -> GravityFieldModel:
        """
                 Object defining the a gravity field model of this body, used to define the exterior gravitational potential, and
                 its gradient(s). Depending on the selected type of model, the type of this attribute
                 is of type GravityFieldModel, or a derived class thereof.
        
        
                 :type: GravityFieldModel
        """

    @gravity_field_model.setter
    def gravity_field_model(self, arg1: GravityFieldModel) -> None:
        ...

    @property
    def ground_station_list(self) -> dict[str, GroundStation]:
        """
                 Dictionary of all ground stations that exist in the body, with dictionary key being the name of the station,
                 and the ground station object the key of the dictionary.
        
        
                 :type: dict[str,GroundStation]
        """

    @property
    def inertia_tensor(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 The current inertia tensor :math:`\\mathbf{I}` of the vehicle, as used in the calculation of
                 (for instance) the response to torques. This attribute is a shorthand for accessing the
                 inertia tensor as computed/stored in the :attr:`~Body.rigid_body_properties` attribute. For certain
                 types of rigid-body properties, this attribute cannot be used to (re)set the current
                 mass.
        
                 Unlike the attributes containing the state, orientation, angular velocity
                 of the Body, this attribute may be used to retrieve the state during the
                 propagation *and* to define the mass of a vehicle.
        
        
                 :type: numpy.ndarray
        """

    @inertia_tensor.setter
    def inertia_tensor(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]) -> None:
        ...

    @property
    def inertial_angular_velocity(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Angular velocity vector of the body, expressed in inertial
                 frame (see :py:attr:`~inertial_to_body_fixed_frame`).
        
        
                 :type: numpy.ndarray
        """

    @property
    def inertial_to_body_fixed_frame(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 **read-only**
        
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
        
        
                 :type: numpy.ndarray
        """

    @property
    def inertial_to_body_fixed_frame_derivative(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 **read-only**
        
                 Time derivative of rotation matrix from inertial frame to this
                 Body's body-fixed frame
                 (see :py:attr:`~inertial_to_body_fixed_frame`).
        
        
                 :type: numpy.ndarray
        """

    @property
    def mass(self) -> float:
        """
                 The current mass :math:`m` of the vehicle, as used in the calculation of
                 non-conservative acceleration. This attribute is a shorthand for accessing the
                 mass as computed/stored in the :attr:`Body.rigid_body_properties` attribute. For certain
                 types of rigid-body properties, this attribute cannot be used to (re)set the current
                 mass. If the body has no    :attr:`Body.rigid_body_properties`, and this function is used to
                 set a mass, a new object is automatically created, with settings analogous to the
                 the :func:`~tudatpy.dynamics.environment_setup.rigid_body.constant_rigid_body_properties` setting.
        
                 Unlike the attributes containing the state, orientation, angular velocity
                 of the Body, this attribute may be used to retrieve the state during the
                 propagation *and* to define the mass of a vehicle.
        
        
                 :type: float
        """

    @mass.setter
    def mass(self, arg1: float) -> None:
        ...

    @property
    def position(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 The translational position of the Body, as set during the
                 current step of the numerical propagation
                 (see :py:attr:`~state`).
        
        
                 :type: numpy.ndarray
        """

    @property
    def radiation_pressure_source_model(self) -> RadiationSourceModel:
        """
                Object that defines the radiation that a body emits, primarily for the calculation of radiation pressure acceleration.
                It computes the irradiance at a given target location.
        
                This attribute is a list of :class:`~RadiationSourceModel`, or a derived class thereof.
        
        
                :type: RadiationSourceModel
        """

    @radiation_pressure_source_model.setter
    def radiation_pressure_source_model(self, arg1: RadiationSourceModel) -> None:
        ...

    @property
    def radiation_pressure_target_models(self) -> list[RadiationPressureTargetModel]:
        """
                List of radiation pressure target models that exist in the body. These objects define how incoming radiation interacts with the body to produce
                a force/torque. A single body may be endowed with multiple target models, which may be selected
                for an acceleration depending on the application. For instance, a body may have a cannonball target model and a panelled target model available,
                and use one for solar radiation pressure acceleration, and the other for planetary radiation pressure acceleration (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure`).
                This attribute is a list of :class:`~RadiationPressureTargetModel`, or a derived class thereof.
        
        
                :type: list[RadiationPressureTargetModel]
        """

    @radiation_pressure_target_models.setter
    def radiation_pressure_target_models(self, arg1: list[RadiationPressureTargetModel]) -> None:
        ...

    @property
    def rigid_body_properties(self) -> RigidBodyProperties:
        """
                Object defining the mass, center of mass and inertia tensor of the body. This object is distinct from
                the gravity field of a body (defined by the :attr:`Body.gravity_field` object). A body endowed with this property does *not*
                automatically have a gravity field created for it. However, the whenever a body is endowed with a gravity field,
                a rigid body properties attribute is created to be consistent with this gravity field (e.g. for a spherical harmonic gravity field
                the mass, center of mass and inertia tensor are created from the gravitational parameter, degree-1 coefficients, and degree-2 coefficients plus mean moment of inertia, respectively).
        
                :type: RigidBodyProperties
        """

    @rigid_body_properties.setter
    def rigid_body_properties(self, arg1: RigidBodyProperties) -> None:
        ...

    @property
    def rotation_model(self) -> RotationalEphemeris:
        """
                 Object defining the orientation of the body, used to calculate the rotation to/from a body-fixed
                 frame (and its derivate). Depending on the selected type of model, the type of this attribute
                 is of type RotationalEphemeris, or a derived class thereof.
        
        
                 :type: RotationalEphemeris
        """

    @rotation_model.setter
    def rotation_model(self, arg1: RotationalEphemeris) -> None:
        ...

    @property
    def shape_model(self) -> BodyShapeModel:
        """
                 Object defining the a shape model of this body, used to define the exterior shape of the body, for instance for
                 the calculation of vehicle's altitude. Depending on the selected type of model, the type of this attribute
                 is of type BodyShapeModel, or a derived class thereof.
        
        
                 :type: BodyShapeModel
        """

    @shape_model.setter
    def shape_model(self, arg1: BodyShapeModel) -> None:
        ...

    @property
    def state(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
                 **read-only**
        
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
        
        
                 :type: numpy.ndarray
        """

    @property
    def system_models(self) -> VehicleSystems:
        """
                 Object used to store physical (hardware) properties of a vehicle, such as engines, control surfaces, etc. This
                 object is typically created automatically whenever such a hardware model needs to be assigned to a vehicle.
        
        
                 :type: VehicleSystems
        """

    @system_models.setter
    def system_models(self, arg1: VehicleSystems) -> None:
        ...

    @property
    def velocity(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 The translational velocity of the Body, as set during the
                 current step of the numerical propagation
                 (see :py:attr:`~state`).
        
        
                 :type: numpy.ndarray
        """

class BodyShapeModel:
    """Object that provides a shape model for a natural body.
    
    Object (typically stored inside a :class:`~Body` object) that provides a shape model for a body, for instance to compute the altitude from a body-centered state, or w.r.t. which
    to place ground stations. This shape model is typically only associated with natural bodies. Shape models for spacecraft (for non-conservative force models) use properties stored inside the
    :class:`~VehicleSystems` object."""

    @property
    def average_radius(self) -> float:
        """
                 **read-only**
        
                 Average radius of the body, for use in computations that assume a spherical body shape.
        
                 :type: float
        """

class CannonballRadiationPressureTargetModel(RadiationPressureTargetModel):
    radiation_pressure_coefficient: float

class ConstantEphemeris(Ephemeris):
    """No documentation found."""

    @typing.overload
    def __init__(self, constant_state_function: typing.Callable[[], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]], reference_frame_origin: str='SSB', reference_frame_orientation: str='ECLIPJ2000') -> None:
        ...

    @typing.overload
    def __init__(self, constant_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], reference_frame_origin: str='SSB', reference_frame_orientation: str='ECLIPJ2000') -> None:
        ...

    def update_constant_state(self, new_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> None:
        """
        No documentation found.
        """

class ConstantTransmittingFrequencyCalculator(TransmittingFrequencyCalculator):

    def __init__(self, frequency: float) -> None:
        ...

class ContinuousInterpolatedMeteoData(StationMeteoData):

    def __init__(self, interpolator: interpolators.OneDimensionalInterpolatorVector, vector_entries: dict[MeteoDataEntries, int]) -> None:
        ...

class ControlSurfaceIncrementAerodynamicInterface:
    """<no_doc, only_dec>"""

class CustomBodyFixedDirectionCalculator(InertialBodyFixedDirectionCalculator):
    inertial_body_axis_direction_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]

class CustomControlSurfaceIncrementAerodynamicInterface(ControlSurfaceIncrementAerodynamicInterface):
    """<no_doc, only_dec>"""

    def __init__(self, coefficient_function: typing.Callable[[list[float]], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]], independent_variable_names: list[aerodynamic_coefficients.AerodynamicCoefficientsIndependentVariables]) -> None:
        ...

class CustomInertialDirectionBasedRotationalEphemeris(RotationalEphemeris):

    @property
    def inertial_body_axis_calculator(self) -> ...:
        ...

class DirectLongitudeLibrationCalculator(LongitudeLibrationCalculator):

    def __init__(self, scaled_libration_amplitude: float) -> None:
        ...

class EarthOrientationAnglesCalculator:
    """Object for computing high-accuracy Earth orientation angles"""

    def get_gcrs_to_itrs_rotation_angles(self, epoch: time_representation.Time, time_scale: time_representation.TimeScales=...) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(5, 1)], time_representation.Time]:
        """
                 Function to compute high-accuracy Earth orientation angles
        
                 Function to compute high-accuracy Earth orientation angle quantities :math:`X,Y,s,x_{p},y_{p}` and UT1 (from which :math:`\\theta_{E}` is computed)
                 as described in :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`
        
                 Parameters
                 ----------
                 epoch : float
                     Epoch at which the Earth orientation angles are to be compute
                 time_scale : TimeScales
                     Time scale in which the input epoch is given
        
                 Returns
                 -------
                 tuple[list[float],float]
                     Pair (tuple of size two) with the first entry a list of orientation angles :math:`X,Y,s,x_{p},y_{p}` (in that order) and the second entry the current UT1.
        """

class EngineModel:

    @property
    def thrust_magnitude_calculator(self) -> ...:
        ...

class Ephemeris:
    """Object that computes the state of a body as a function of time
    
    
    Object (typically stored inside a :class:`~Body` object) that computes the state of a body as a function of time,
    both outside of a propagation, and during a propagation if the given body's translational state is not propagated.
    Note that this object computes the state w.r.t. its own origin (defined by ``frame_origin``), which need not be the same as the global frame origin
    of the environment."""

    def cartesian_position(self, current_time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 As ``cartesian_state``, but only the three position components
        
        
                 Parameters
                 ----------
                 current_time : float
                     Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.
        
                 Returns
                 -------
                 numpy.ndarray
                     Requested Cartesian position
        """

    def cartesian_state(self, current_time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
                 This function returns the Cartesian state (position and velocity) at the given time, w.r.t. the ``frame_origin``.
        
        
                 Parameters
                 ----------
                 current_time : float
                     Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.
        
                 Returns
                 -------
                 numpy.ndarray
                     Requested Cartesian state
        """

    def cartesian_velocity(self, current_time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 As ``cartesian_state``, but only the three velocity components
        
        
                 Parameters
                 ----------
                 current_time : float
                     Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.
        
                 Returns
                 -------
                 numpy.ndarray
                     Requested Cartesian velocity
        """

    @property
    def frame_orientation(self) -> str:
        """
                 **read-only**
        
                 Name of the frame orientation w.r.t which this object provides its states
        
        
        
                 :type: str
        """

    @property
    def frame_origin(self) -> str:
        """
                 **read-only**
        
                 Name of the reference body/point w.r.t. which this object provides its states
        
        
                 :type: str
        """

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
    of a single body w.r.t. a single central body."""

    def update_conditions(self, current_time: float) -> None:
        ...

    @property
    def aerodynamic_angle_calculator(self) -> AerodynamicAngleCalculator:
        """
                 **read-only**
        
                 The object that is responsible for computing various relevant
                 flight dynamics angles and frame rotations.
        
        
                 :type: AerodynamicAngleCalculator
        """

    @property
    def altitude(self) -> float:
        """
                 **read-only**
        
                 The current time, at which this object was last updated
        
        
                 :type: float
        """

    @property
    def body_centered_body_fixed_state(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
                 **read-only**
        
                 Cartesian translational state, expressed in a frame centered
                 at, and fixed to, the central body. Note that, due to the
                 rotation of the central body, the norm of the body-fixed,
                 body-centered, velocity differs from the norm of the inertial
                 body-centered velocity.
        
        
                 :type: numpy.ndarray
        """

    @property
    def geodetic_latitude(self) -> float:
        """
                 **read-only**
        
                 The body-fixed geographic latitude of the body w.r.t. its
                 central body.
        
        
                 :type: float
        """

    @property
    def latitude(self) -> float:
        """
                 **read-only**
        
                 The body-fixed geographic latitude of the body w.r.t. its
                 central body.
        
        
                 :type: float
        """

    @property
    def longitude(self) -> float:
        """
                 **read-only**
        
                 The body-fixed longitude of the body w.r.t. its central body.
        
        
                 :type: float
        """

    @property
    def time(self) -> float:
        """
                 **read-only**
        
                 The current time, at which this object was last updated
        
        
                 :type: float
        """

class GcrsToItrsRotationModel(RotationalEphemeris):
    """Object for high-accuracy GCRS<->ITRS rotation.
    
    Object derived from :class:`~RotationalEphemeris` that implements the high-accuracy GCRS<->ITRS rotation as per the IERS 2010 Conventions. The details of the model are described in
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`
    With the exception of :math:`s'`, the list of angles used to compute the full rotation are computed by an object of type :class:`~EarthOrientationAnglesCalculator` (which can be retrieved from this rotation model
    through :attr:`~GcrsToItrsRotationModel.angles_calculator`."""

    @property
    def angles_calculator(self) -> EarthOrientationAnglesCalculator:
        """
                 **read-only**
        
                 Object that computes the Earth rotation angles :math:`X,Y,s,\\theta_{E},x_{p},y_{p}`
        
        
                 :type: EarthOrientationAnglesCalculator
        """

class GravityFieldModel:
    """Object that provides the gravity field of a body
    
    
    Object (typically stored inside a :class:`~Body` object) that provides the gravity field of a body, typically (but not exclusively) for
    use in gravitational acceleration and torque models. This base class allows access to the gravitational parameter of the body.
    Specific derived classes are implemented to provide models for more detailed gravity field models (e.g. spherical harmonics, polyhedron)."""

    def __init__(self, gravitational_parameter: float, update_inertia_tensor: typing.Callable[[], None]=None) -> None:
        ...

    def get_gravitational_parameter(self) -> float:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
                 Value of the gravity field's gravitational parameters :math:`\\mu`
        
        
                 :type: float
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

class GravityFieldVariationModel:
    """Object that computes a single type of gravity field variation.
    
    Object that computes a single type of gravity field variation. This object is typically not used directly, but internally by the :class:`~TimeDependentSphericalHarmonicsGravityField` class."""

class GroundStation:
    """Object used to define and store properties of a ground station.
    
    Object (typically stored inside a :class:`~Body` object) used to define and store properties of a ground station, typically used in modelling tracking observations to/from a ground station."""

    def set_station_meteo_data(self, meteo_data: ...) -> None:
        ...

    def set_timing_system(self, timing_system: TimingSystem) -> None:
        ...

    def set_transmitting_frequency_calculator(self, transmitting_frequency_calculator: ...) -> None:
        ...

    @property
    def pointing_angles_calculator(self) -> ...:
        """
                 **read-only**
        
                 Object that performs computations of the azimuth and elevation of an arbitrary target as observed by the ground station
        
                 :type: PointingAnglesCalculator
        """

    @property
    def pressure_function(self) -> typing.Callable[[float], float]:
        """
                 Function that provides the local pressure at the ground station (typically use for media corrections) as a function of time
        
                 :type: :type: Callable[[float], float]
        """

    @property
    def relative_humidity_function(self) -> typing.Callable[[float], float]:
        """
                 Function that provides the local relative humidity at the ground station (typically use for media corrections) as a function of time
        
                 :type: :type: Callable[[float], float]
        """

    @property
    def station_state(self) -> GroundStationState:
        """
                 **read-only**
        
                 Object that performs computations of the current (body-fixed) position and frame conversions of the ground station.
        
                 :type: GroundStationState
        """

    @property
    def temperature_function(self) -> typing.Callable[[float], float]:
        """
                 Function that provides the local temperature at the ground station (typically use for media corrections) as a function of time
        
                 :type: :type: Callable[[float], float]
        """

    @property
    def transmitting_frequency_calculator(self) -> ...:
        """
                 Object that provides the transmission frequency as a function of time for (radio) tracking stations. This attribute is typically set automatically when loading tracking data files (e.g. ODF, IFMS, TNF, etc.)
        
                 :type: TransmittingFrequencyCalculator
        """

    @transmitting_frequency_calculator.setter
    def transmitting_frequency_calculator(self, arg1: ...) -> None:
        ...

class GroundStationState:
    """Object that performs computations of the current (body-fixed) position and frame conversions of the ground station.
    
    Object that performs computations of the current (body-fixed) position and frame conversions of the ground station. In the simplest situation,
    only a Cartesian position is provided, which is then assumed constant. If time variations (for instance due to tides or plate motion) are present,
    their impact on station position is computed in this object."""

    def get_cartesian_position(self, current_time: float, target_frame_origin: str='') -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 This function computes the position of the station as a function of time.
        
                 This function computes the position of the station as a function of time, in a frame with body-fixed orientation.
                 Some time-variations of the station position depend on the *origin* of the frame in which the computation is to be
                 used. For instance, relativistic correction to the Earth-fixed position is different in a geocentric or barycentric frame.
                 However, the output of this function is always given in the body-fixed, body-centered frame.
        
                 Parameters
                 ----------
                 current_time : float
                     Time (in seconds since J2000 TDB) at which the position is to be computed.
        
                 target_frame_origin: str, default = ""
                     Identifier for the frame origin w.r.t. which the computed position is to be used.
        
                 Returns
                 -------
                 numpy.ndarray
                     Cartesian position of the station at the current epoch, in a body-centered, body-fixed frame
        """

    def get_cartesian_state(self, current_time: float, target_frame_origin: str='') -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        ...

    @property
    def cartesian_positon_at_reference_epoch(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Cartesian position of the ground station, at the reference epoch, in a body-fixed, body-centered frame.
        
                 :type: numpy.ndarray[numpy.float64[3, 1]]
        """

    @property
    def geodetic_position_at_reference_epoch(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Geodetic position of the ground station (altitude w.r.t. body shape model, geodetic latitude, longitude), at the reference epoch, in a body-fixed, body-centered frame.
        
                 :type: numpy.ndarray[numpy.float64[3, 1]]
        """

    @property
    def rotation_matrix_body_fixed_to_topocentric(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 **read-only**
        
                 Rotation matrix from the body-fixed frame (of the station's central body) to the topocentric frame
                 of the ground station. The body-fixed frame is defined by the rotation model of the body object (:attr:`~tudatpy.dynamics.environment.Body.rotation_model`).
                 The axes of the topocentric frame are defined such that the x-axis is in East direction, the z-direction is upwards, perpendicular to the body's surface sphere
                 (with properties defined by the central body's shape model :attr:`~tudatpy.dynamics.environment.Body.shape_model`). The y-axis completes the frame, and is in northern direction.
                 For time-varying ground station positions, this function uses the station position at reference epoch for the computation of the axes.
        
                 :type: numpy.ndarray[numpy.float64[3, 3]]
        """

    @property
    def spherical_positon_at_reference_epoch(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 **read-only**
        
                 Spherical position of the ground station (distance w.r.t. body center, latitude, longitude), at the reference epoch, in a body-fixed, body-centered frame.
        
                 :type: numpy.ndarray[numpy.float64[3, 1]]
        """

class HypersonicLocalInclinationAnalysis(AerodynamicCoefficientGenerator36):

    def __init__(self, independent_variable_points: list[list[float]], body_shape: geometry.SurfaceGeometry, number_of_lines: list[int], number_of_points: list[int], invert_orders: list[bool], selected_methods: list[list[int]], reference_area: float, reference_length: float, moment_reference_point: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], save_pressure_coefficients: bool=False) -> None:
        """
                 Class constructor, taking the shape of the vehicle, and various analysis options as input.
        
        
                 Parameters
                 ----------
                 independent_variable_points : list[list[float]]
                     List containing three lists, with each sublist containing the data points of each of the
                     independent variables for the coefficient generation. The physical meaning of each of the
                     three independent variables is: 0 = mach number, 1 = angle of attack, 2 = angle of sideslip.
                     Each of the subvectors must be sorted in ascending order.
        
                 body_shape : SurfaceGeometry
                     Class that defines the shape of the vehicle as a continuous surface. The local inclination analysis
                     discretizes the surface of the vehicle into quadrilateral panels, defined by the other inputs to
                     this constructor. In case the :class:`tudat.geometry.SurfaceGeometry` object is made up of multiple
                     sub-shapes, different settings may be used for each
        
                 number_of_lines : List[ float ]
                     Number of discretization points in the first independent surface variable of each of the subparts of body_shape.
                     The size of this list should match the number of parts of which the body_shape is composed. The first independent
                     variable of a subpart typically runs along the longitudinal vehicle direction
        
                 number_of_points : List[ float ]
                     Number of discretization points in the second independent surface variable of each of the subparts of body_shape.
                     The size of this list should match the number of parts of which the body_shape is composed. The first independent
                     variable of a subpart typically runs along the lateral vehicle direction
        
                 invert_orders : List[ bool ]
                     Booleans to denote whether the surface normals of the panels of each discretized body_shape subpart are to be inverted
                     (i.e. inward-facing->outward facing or vice versa). The size of this list should match the number of parts of which the body_shape is composed.
        
                 selected_methods : List[ List[ int ] ]
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
        
                 reference_area : float
                     Reference area used to non-dimensionalize aerodynamic forces and moments.
        
                 moment_reference_point : numpy.ndarray
                     Reference point wrt which aerodynamic moments are calculated.
        
                 save_pressure_coefficients : bool
                     Boolean denoting whether to save the pressure coefficients that are computed to files
        """

    def clear_data(self) -> None:
        ...

class InertialBodyFixedDirectionCalculator:
    pass

class IonosphereModel:
    """Base class for ionospheric models.
    
    Provides the vertical total electron content (VTEC) in TECU (1 TECU = 1e16 e-/m²) based on geodetic position (latitude, longitude) and time.
    
    This is the base class from which models like TabulatedIonosphereModel or GlobalIonosphereModelVtecCalculator retrieve electron content data. The model is typically stored
    inside a `Body` instance and used in observation corrections or environmental queries."""

class KeplerEphemeris(Ephemeris):
    pass

class LongitudeLibrationCalculator:
    pass

class MeteoDataEntries:
    """Members:
    
      temperature_meteo_data
    
      pressure_meteo_data
    
      water_vapor_pressure_meteo_data
    
      relative_humidity_meteo_data
    
      dew_point_meteo_data"""
    __members__: typing.ClassVar[dict[str, MeteoDataEntries]]
    dew_point_meteo_data: typing.ClassVar[MeteoDataEntries]
    pressure_meteo_data: typing.ClassVar[MeteoDataEntries]
    relative_humidity_meteo_data: typing.ClassVar[MeteoDataEntries]
    temperature_meteo_data: typing.ClassVar[MeteoDataEntries]
    water_vapor_pressure_meteo_data: typing.ClassVar[MeteoDataEntries]

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

class MultiArcEphemeris(Ephemeris):

    def __init__(self, single_arc_ephemerides: dict[float, Ephemeris], reference_frame_origin: str='SSB', reference_frame_orientation: str='ECLIPJ2000') -> None:
        ...

class PiecewiseLinearFrequencyInterpolator(TransmittingFrequencyCalculator):

    def __init__(self, start_times: list[time_representation.Time], end_times: list[time_representation.Time], ramp_rates: list[float], start_frequency: list[float]) -> None:
        ...

class PointingAnglesCalculator:

    def calculate_azimuth_angle(self, inertial_vector_to_target: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], time: float) -> float:
        ...

    def calculate_elevation_angle(self, inertial_vector_to_target: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], time: float) -> float:
        ...

    def convert_inertial_vector_to_topocentric(self, inertial_vector: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        ...

class PolyhedronGravityField(GravityFieldModel):

    @property
    def vertices_coordinates(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        ...

    @property
    def vertices_defining_each_facet(self) -> typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        ...

    @property
    def volume(self) -> float:
        ...

class RadiationPressureTargetModel:
    pass

class RadiationSourceModel:
    pass

class RigidBodyProperties:
    """Object that defines the mass, center of mass, and inertia tensor as a function of time.
    
    Object that defines the mass, center of mass, and inertia tensor as a function of time, typically used for evaluation of torques and non-conservative forces
    in numerical state propagation. Note that this object does *not* define properties of a gravity field (it defines the inertial mass rather than the gravitational mass)"""

    def update(self, time: float) -> None:
        """
                 Function to update the body properties to the current time. This function is called automatically during a propagation loop.
                 In case these properties are not time-dependent (e.g. when using the :func:`~tudatpy.dynamics.environment_setup.rigid_body.constant_rigid_body_properties` settings)
                 this function does nothing (since no update is needed).
        
                 Parameters
                 ----------
                 current_time : float
                     Time (in seconds since J2000 TDB) to which this object is to be updated
        
                 Returns
                 -------
        """

    @property
    def current_center_of_mass(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 Position of the center of mass of the object (in the body-centered, body-fixed frame), as set by the latest call to the ``update`` function of this object.
        """

    @property
    def current_inertia_tensor(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                Inertia tensor of the object (with axes along those of the body-fixed frame), as set by the latest call to the ``update`` function of this object.
        """

    @property
    def current_mass(self) -> float:
        """
                 Mass of the object, as set by the latest call to the ``update`` function of this object.
        """

class RotationalEphemeris:
    """Object that stores the rotational state of the bodies.
    
    
    Object that stores the rotational state of the bodies. This object can be used to calculate rotation matrices,
    which are used to transform coordinates between reference frames."""

    def angular_velocity_in_body_fixed_frame(self, time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 Function to get the body's angular velocity vector, expressed in the body-fixed frame.
        
        
                 Function to get the body's angular velocity vector :math:`\\boldsymbol{\\omega}^{(B)}`, expressed in the body-fixed frame :math:`B`.
                 The calculation of the angular velocity depends on the specific rotation model that has been defined,
                 either from an a priori definition (see :ref:`rotation_model` submodule) or from processing
                 the results of propagation of the rotational equations of motion.
                 Note that when numerically propagating rotational dynamics, this angular velocity vector is typically directly defined
                 in the last three entries of the state vector.
        
        
                 Parameters
                 ----------
                 current_time : float
                     The time (in seconds since epoch J2000, TDB time scale) at which the angular velocity vector is evaluated
        
                 Returns
                 -------
                 numpy.ndarray
                     Angular velocity vector of the body  :math:`\\boldsymbol{\\omega}^{(B)}` expressed in the body-fixed frame :math:`B`
        """

    def angular_velocity_in_inertial_frame(self, time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
        """
                 Function to get the body's angular velocity vector, expressed in the inertial frame.
        
        
                 Function to get the body's angular velocity vector :math:`\\boldsymbol{\\omega}^{(I)}`, expressed in the body-fixed frame :math:`I`.
                 This quantity is computed from :math:`\\mathbf{R}^{I/B}\\boldsymbol{\\omega}^{(B)}`, see the ``angular_velocity_in_body_fixed_frame`` and
                 ``body_fixed_to_inertial_rotation`` functions.
        
        
                 Parameters
                 ----------
                 current_time : float
                     The time (in seconds since epoch J2000, TDB time scale) at which the angular velocity vector is evaluated
        
                 Returns
                 -------
                 numpy.ndarray
                     Angular velocity vector of the body  :math:`\\boldsymbol{\\omega}^{(B)}` expressed in the body-fixed frame :math:`B`
        """

    def body_fixed_to_inertial_rotation(self, time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 Function to get rotation matrix from body-fixed frame to inertial frame over time.
        
        
                 Function to get rotation matrix from body-fixed (target) frame to inertial (base) frame over time.
                 The calculation of this rotation matrix depends on the specific rotation model that has been defined,
                 either from an a priori definition (see :ref:`rotation_model` submodule) or from processing
                 the results of propagation of the rotational equations of motion.
        
        
                 Parameters
                 ----------
                 current_time : float
                     The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix is evaluated
        
                 Returns
                 -------
                 numpy.ndarray
                     Rotation matrix :math:`\\mathbf{R}^{I/B}` from body-fixed frame :math:`B` to inertial frame `I`
        """

    def inertial_to_body_fixed_rotation(self, time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 Function to get rotation matrix from inertial frame to body-fixed frame over time.
        
        
                 Function computes the inverse (equal to transpose) rotation of the ``body_fixed_to_inertial_rotation`` function.
        
        
                 Parameters
                 ----------
                 current_time : float
                     The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix is evaluated
        """

    def time_derivative_body_fixed_to_inertial_rotation(self, time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 Function to get time derivative of rotation matrix from body-fixed frame to inertial frame over time.
        
        
                 Function to get time derivative of rotation matrix from body-fixed frame to inertial frame over time (see ``body_fixed_to_inertial_rotation``),
                 denoted :math:`\\dot{\\mathbf{R}}^{(I/B)}`,
        
        
                 Parameters
                 ----------
                 current_time : float
                     The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix derivative is evaluated
        
                 Returns
                 -------
                 numpy.ndarray
                     Rotation matrix :math:`\\dot{\\mathbf{R}}^{I/B}` from body-fixed frame :math:`B` to inertial frame `I`
        """

    def time_derivative_inertial_to_body_fixed_rotation(self, time: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
        """
                 Function to get time derivative of rotation matrix from inertial frame to body-fixed frame over time.
        
        
                 Function to get time derivative of rotation matrix from inertial frame to body-fixed frame over time (see ``inertial_to_body_fixed_rotation``),
                 denoted :math:`\\dot{\\mathbf{R}}^{(B/I)}`,
        
        
                 Parameters
                 ----------
                 current_time : float
                     The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix derivative is evaluated
        
                 Returns
                 -------
                 numpy.ndarray
                     Rotation matrix :math:`\\dot{\\mathbf{R}}^{B/I}` from inertial frame `I` to body-fixed frame :math:`B`
        """

    @property
    def body_fixed_frame_name(self) -> str:
        """
                 **read-only**
        
                 The identifier of the body-fixed frame, used in other parts of the simulation to identify it.
        
        
                 :type: str
        """

    @property
    def inertial_frame_name(self) -> str:
        """
                 **read-only**
        
                 The identifier of the inertial frame, used in other parts of the simulation to identify it.
        
        
                 :type: str
        """

class SphericalHarmonicsGravityField(GravityFieldModel):
    """Object that provides a spherical harmonic gravity field of a body.
    
    Object (typically stored inside a :class:`~Body` object) that provides a spherical harmonic gravity field of a body, typically (but not exclusively) for
    use in gravitational acceleration and torque models. This class is derived from :class:`~GravityFieldModel`.  This object is typically created using the :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic`
    settings function. If any time variations of the gravity field are provided, an object of the derived class :class:`~TimeVariableSphericalHarmonicsGravityField` is created."""

    @property
    def cosine_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Matrix with cosine spherical harmonic coefficients :math:`\\bar{C}_{lm}` (geodesy normalized). Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.
        
                 :type: numpy.ndarray[numpy.float64[l, m]]
        """

    @cosine_coefficients.setter
    def cosine_coefficients(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

    @property
    def maximum_degree(self) -> float:
        """
                 **read-only**
        
                 Maximum spherical harmonic degree :math:`l_{max}` for which coefficients are defined
        
                 :type: int
        """

    @property
    def maximum_order(self) -> float:
        """
                 **read-only**
        
                 Maximum spherical harmonic order :math:`m_{max}` for which coefficients are defined
        
                 :type: int
        """

    @property
    def reference_radius(self) -> float:
        """
                 **read-only**
        
                 Reference radius :math:`R` of the gravity field
        
                 :type: float
        """

    @property
    def sine_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Matrix with sine spherical harmonic coefficients :math:`\\bar{S}_{lm}` (geodesy normalized). Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.
        
                 :type: numpy.ndarray[numpy.float64[l, m]]
        """

    @sine_coefficients.setter
    def sine_coefficients(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

class StationMeteoData:
    pass

class SynchronousRotationalEphemeris(RotationalEphemeris):
    libration_calculator: LongitudeLibrationCalculator

class SystemOfBodies:
    """Object that contains a set of Body objects and associated frame
    information.
    
    
    Object that contains a set of Body objects and associated frame
    information. This object stored the entire environment for a
    typical Tudat numerical simulation, and is fundamental for the
    overall Tudat architecture."""

    def add_body(self, body_to_add: Body, body_name: str, process_body: bool=1) -> None:
        """
                 This function adds an existing body, which the user has
                 separately created, to the :py:class:`~SystemOfBodies`.
        
        
        
                 Parameters
                 ----------
                 body_to_add : Body
                     Body object that is to be added.
        
                 body_name : numpy.ndarray
                     Name of the Body that is to be added.
        
                 process_body : bool, default=True
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
        
        
                 Parameters
                 ----------
                 body_name : string
                     Name of the Body that is to be added
        
                 process_body : bool, default=True
                     Variable that defines whether this new Body will have its
                     global frame origin/orientation set to conform to rest of
                     the environment.
        
                     .. warning:: Only in very rare cases should
                                  this variable be anything other than ``True``.
                                  Users are recommended to keep this default value
                                  intact.
        
        
        
        
        
                 Examples
                 --------
        
                 This function is often used early on in the environment
                 creation segment of a simulation, following the creation of
                 a :py:class:`~SystemOfBodies` from the default settings
                 for celestial bodies.
        
                 .. code-block:: python
                    :emphasize-lines: 18
        
                    # Define string names for bodies to be created from default.
                    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]
        
                    # Use "Earth"/"J2000" as global frame origin and orientation.
                    global_frame_origin = "Earth"
                    global_frame_orientation = "J2000"
        
                    # Create default body settings, usually from `spice`.
                    body_settings = environment_setup.get_default_body_settings(
                        bodies_to_create,
                        global_frame_origin,
                        global_frame_orientation)
        
                    # Create system of selected celestial bodies
                    bodies = environment_setup.create_system_of_bodies(body_settings)
        
                    # Create vehicle objects.
                    bodies.create_empty_body("Delfi-C3")
        """

    def does_body_exist(self, body_name: str) -> bool:
        """
                 Function to check if a body with a given name exists in the SystemOfBodies
        
                 Parameters
                 ----------
                 body_name : string
                     Name of the Body whose existence is to be checked
        
                 Returns
                 -------
                 bool
                     True if the body exists in this object, false if not
        """

    def get(self, body_name: str) -> Body:
        """
                 This function extracts a single Body object from the SystemOfBodies.
        
        
                 Parameters
                 ----------
                 body_name : str
                     Name of the Body that is to be retrieved.
        
                 Returns
                 -------
                 Body
                     Body object of the requested name
        """

    def get_body(self, body_name: str) -> Body:
        """
                 Deprecated version of :py:func:`~get`
        """

    def global_frame_orientation(self) -> str:
        """
                 Common global frame orientation for all bodies in this SystemOfBodies, described in more detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#frame-orientation>`__.
        """

    def global_frame_origin(self) -> str:
        """
                 Common global frame origin for all bodies in this SystemOfBodies, described in more detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#global-origin>`__.
        """

    def list_of_bodies(self) -> list[str]:
        """
                 List of names of bodies that are stored in this SystemOfBodies
        """

    def remove_body(self, body_name: str) -> None:
        """
                 This function removes an existing body from the
                 :py:class:`~SystemOfBodies`.
        
        
        
                 .. warning:: This function does *not* necessarily delete the
                              Body object, it only removes it from this object.
                              If any existing models in the simulation refer to
                              this Body, it will persist in memory.
        
        
                 Parameters
                 ----------
                 body_name : numpy.ndarray
                     Name of the Body that is to be removed.
        """

class TabulatedEphemeris(Ephemeris):
    interpolator: interpolators.OneDimensionalInterpolatorVector

class TimeDependentSphericalHarmonicsGravityField(SphericalHarmonicsGravityField):
    """Derived class of :class:`~SphericalHarmonicsGravityField` that is created when any gravity field variations are detected.
    
    Derived class of :class:`~SphericalHarmonicsGravityField` that is created instead when any gravity field variations are detected during object creation
    (typically during a call of :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`)
    This object computes the time-variability of spherical harmonic coefficients from a list of :class:`~GravityFieldVariationModel` objects.
    The ``cosine_coefficients`` and ``sine_coefficients`` attributes provide the instantaneous coefficients (including the time-variability)
    The ``nominal_cosine_coefficients`` and ``nominal_sine_coefficients`` provide the static (e.g. without time-variations) coefficients."""

    @property
    def gravity_field_variation_models(self) -> list[...]:
        """
                 **read-only**
        
                 List of gravity field variation models that the object uses to update the spherical harmonic coefficients at every time step
        
                 :type: list[GravityFieldVariationModel]
        """

    @property
    def nominal_cosine_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Matrix with cosine spherical harmonic coefficients :math:`\\bar{C}_{lm}` (geodesy normalized) *excluding* time-variations. Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.
        
                 :type: numpy.ndarray[numpy.float64[l, m]]
        """

    @nominal_cosine_coefficients.setter
    def nominal_cosine_coefficients(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

    @property
    def nominal_sine_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Matrix with sine spherical harmonic coefficients :math:`\\bar{S}_{lm}` (geodesy normalized) *excluding* time-variations. Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.
        
                 :type: numpy.ndarray[numpy.float64[l, m]]
        """

    @nominal_sine_coefficients.setter
    def nominal_sine_coefficients(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

class TimingSystem:
    """No documentation found."""

    @typing.overload
    def __init__(self, arc_times: list[time_representation.Time], all_arcs_polynomial_drift_coefficients: list[float]=[], clock_noise_generation_function: typing.Callable[[float, float, float], typing.Callable[[float], float]]=None, clock_noise_time_step: float=0.001) -> None:
        ...

    @typing.overload
    def __init__(self, arc_times: list[time_representation.Time], polynomial_drift_coefficients: list[list[float]], clock_noise_generation_function: typing.Callable[[float, float, float], typing.Callable[[float], float]]=None, clock_noise_time_step: float=0.001) -> None:
        ...

    @typing.overload
    def __init__(self, polynomial_drift_coefficients: list[list[float]], stochastic_clock_noise_functions: list[typing.Callable[[float], float]], arc_times: list[time_representation.Time]) -> None:
        ...

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

class TransmittingFrequencyCalculator:
    """No documentation found."""

class VehicleSystems:
    """Object used to store physical (hardware) properties of a vehicle."""

    def __init__(self) -> None:
        ...

    def get_control_surface_deflection(self, control_surface_id: str) -> float:
        """
                 Function to retrieve the current deflection of an aerodynamic control surface.
        
        
                 Function to retrieve the current deflection of an aerodynamic control surface,
                 identified by its name. To extract the control surface deflection, the control
                 surface has to exist. A control surface is created whenever control surfaces are
                 defined in a body's aerodynamic coefficient interface.
        
        
                 Parameters
                 ----------
                 control_surface_id : str
                     The identified (name) of the given control surface
        
                 Returns
                 -------
                 float
                     Current deflection (in radians) that the control surface
        """

    def get_engine_model(self, engine_name: str) -> ...:
        """
                 Function to retrieve an engine model from the vehicle
        
        
        
                 Parameters
                 ----------
                 engine_name : str
                     The identifier for the engine model that is to be retrieved
        
                 Returns
                 -------
                 EngineModel
                     Model for the engine that is requested
        """

    def set_control_surface_deflection(self, control_surface_id: str, deflection_angle: float) -> None:
        """
                 Function to set the current deflection of an aerodynamic control surface.
        
        
                 Function to set the current deflection of an aerodynamic control surface,
                 identified by its name. To set the control surface deflection, the control
                 surface has to exist. A control surface is created whenever control surfaces are
                 defined in a body's aerodynamic coefficient interface.
        
        
                 Parameters
                 ----------
                 control_surface_id : str
                     The identified (name) of the given control surface
        
                 deflection_angle : float
                     The deflection (in radians) that the control surface is to be set to. This will
                     typically influence the aerodynamic coefficients of the vehicle
        """

    def set_default_transponder_turnaround_ratio_function(self) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_reference_point(self, reference_point: str, location: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], frame_origin: str='', frame_orientation: str='') -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_reference_point(self, reference_point: str, ephemeris: ...) -> None:
        """
        No documentation found.
        """

    def set_timing_system(self, timing_system: ...) -> None:
        ...

    def set_transponder_turnaround_ratio(self, transponder_ratio_per_uplink_and_downlink_frequency_band: dict[tuple[..., ...], float]) -> None:
        """
        No documentation found.
        """

def get_default_local_inclination_angle_of_attack_points() -> list[float]:
    ...

def get_default_local_inclination_mach_points(mach_regime: str='Full') -> list[float]:
    ...

def get_default_local_inclination_sideslip_angle_points() -> list[float]:
    ...

def get_local_inclination_mesh(local_inclination_analysis_object: HypersonicLocalInclinationAnalysis) -> tuple[list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]]:
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
    
    
    Parameters
    ----------
    local_inclination_analysis_object : HypersonicLocalInclinationAnalysis
        Object used to calculate the aerodynamics of the vehicle
    
    output_directory : str
        Directory to which the files are to be saved
    
    output_file_prefix : str, default=''
        Optional prefix of output file names"""

def transform_to_inertial_orientation(state_in_body_fixed_frame: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], current_time: float, rotational_ephemeris: RotationalEphemeris) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Function to convert a Cartesian state vector from a body-fixed to an inertial frame
    
    Function to convert a Cartesian state vector from a body-fixed to an inertial frame, using a :class:`~tudatpy.dynamics.environment.RotationalEphemeris`
    object as a model for the rotation. The body-fixed frame from which the conversion takes place is the :attr:`~tudatpy.dynamics.environment.RotationalEphemeris.body_fixed_frame_name` frame,
    the (assumedly) inertial frame to which the conversion is done is :attr:`~tudatpy.dynamics.environment.RotationalEphemeris.inertial_frame_name`.
    
    This function calls :func:`~tudatpy.astro.element_conversion.rotate_state_to_frame` (with frame :math:`A` the inertial frame, and frame :math:`B` the body-fixed frame). The present function
    computes the required rotation matrix and its time derivative from the ``rotational_ephemeris`` input given here.
    
    Parameters
    ----------
    state_in_body_fixed_frame : numpy.ndarray[numpy.float64[6, 1]]
        Cartesian state (position and velocity) in the body-fixed frame
    
    current_time : float
        Time at which the transformation is to be computed
    
    rotational_ephemeris : RotationalEphemeris
        Boy rotation model that is to be used to convert the body-fixed state to inertial state
    
    Returns
    -------
    numpy.ndarray[numpy.float64[6, 1]]
        Cartesian state transformed to inertial frame, using ``rotational_ephemeris`` model, from body-fixed ``state_in_body_fixed_frame``"""
dew_point_meteo_data: MeteoDataEntries
pressure_meteo_data: MeteoDataEntries
relative_humidity_meteo_data: MeteoDataEntries
temperature_meteo_data: MeteoDataEntries
water_vapor_pressure_meteo_data: MeteoDataEntries