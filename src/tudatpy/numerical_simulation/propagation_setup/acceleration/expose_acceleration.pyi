import numpy
import tudatpy.numerical_simulation.environment_setup.radiation_pressure.expose_radiation_pressure
import tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust
import typing
__all__ = ['AccelerationSettings', 'AvailableAcceleration', 'CustomAccelerationSettings', 'DirectTidalDissipationAccelerationSettings', 'EmpiricalAccelerationSettings', 'MomentumWheelDesaturationAccelerationSettings', 'MutualSphericalHarmonicAccelerationSettings', 'RelativisticAccelerationCorrectionSettings', 'SphericalHarmonicAccelerationSettings', 'ThrustAccelerationSettings', 'aerodynamic', 'aerodynamic_type', 'cannonball_radiation_pressure', 'cannonball_radiation_pressure_type', 'custom', 'custom_acceleration', 'custom_acceleration_type', 'direct_tidal_dissipation_acceleration', 'direct_tidal_dissipation_in_central_body_acceleration_type', 'direct_tidal_dissipation_in_orbiting_body_acceleration_type', 'empirical', 'empirical_acceleration_type', 'mutual_spherical_harmonic_gravity', 'mutual_spherical_harmonic_gravity_type', 'point_mass_gravity', 'point_mass_gravity_type', 'polyhedron_gravity', 'polyhedron_gravity_type', 'quasi_impulsive_shots_acceleration', 'quasi_impulsive_shots_acceleration_type', 'radiation_pressure', 'radiation_pressure_type', 'relativistic_correction', 'relativistic_correction_acceleration_type', 'ring_gravity', 'ring_gravity_type', 'spherical_harmonic_gravity', 'spherical_harmonic_gravity_type', 'thrust_acceleration_type', 'thrust_and_isp_from_custom_function', 'thrust_from_all_engines', 'thrust_from_custom_function', 'thrust_from_direction_and_magnitude', 'thrust_from_engine', 'thrust_from_engines', 'undefined_acceleration_type', 'yarkovsky']

class AccelerationSettings:
    """Functional base class to define settings for accelerations.
	
	Class for providing settings for acceleration model. This class is a functional (base) class for
	settings of acceleration models that  require no information in addition to their type.
	Classes defining settings for acceleration models requiring additional information must be derived from this class.
	Bodies exerting and undergoing acceleration are set externally from this class.
	This class can be used for the easy setup of acceleration models
	(see createAccelerationModels.h), but users may also chose to do so manually.
	(Derived) Class members are all public, for ease of access and modification.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class AvailableAcceleration:
    """Members:
	
	undefined_acceleration_type :
	
	point_mass_gravity_type :
	
	aerodynamic_type :
	
	cannonball_radiation_pressure_type :
	
	spherical_harmonic_gravity_type :
	
	mutual_spherical_harmonic_gravity_type :
	
	polyhedron_gravity_type :
	
	ring_gravity_type :
	
	thrust_acceleration_type :
	
	relativistic_correction_acceleration_type :
	
	empirical_acceleration_type :
	
	direct_tidal_dissipation_in_central_body_acceleration_type :
	
	direct_tidal_dissipation_in_orbiting_body_acceleration_type :
	
	quasi_impulsive_shots_acceleration_type :
	
	custom_acceleration_type :
	
	radiation_pressure_type :
	"""
    __members__: typing.ClassVar[dict[str, AvailableAcceleration]]
    aerodynamic_type: typing.ClassVar[AvailableAcceleration]
    cannonball_radiation_pressure_type: typing.ClassVar[AvailableAcceleration]
    custom_acceleration_type: typing.ClassVar[AvailableAcceleration]
    direct_tidal_dissipation_in_central_body_acceleration_type: typing.ClassVar[AvailableAcceleration]
    direct_tidal_dissipation_in_orbiting_body_acceleration_type: typing.ClassVar[AvailableAcceleration]
    empirical_acceleration_type: typing.ClassVar[AvailableAcceleration]
    mutual_spherical_harmonic_gravity_type: typing.ClassVar[AvailableAcceleration]
    point_mass_gravity_type: typing.ClassVar[AvailableAcceleration]
    polyhedron_gravity_type: typing.ClassVar[AvailableAcceleration]
    quasi_impulsive_shots_acceleration_type: typing.ClassVar[AvailableAcceleration]
    radiation_pressure_type: typing.ClassVar[AvailableAcceleration]
    relativistic_correction_acceleration_type: typing.ClassVar[AvailableAcceleration]
    ring_gravity_type: typing.ClassVar[AvailableAcceleration]
    spherical_harmonic_gravity_type: typing.ClassVar[AvailableAcceleration]
    thrust_acceleration_type: typing.ClassVar[AvailableAcceleration]
    undefined_acceleration_type: typing.ClassVar[AvailableAcceleration]

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

class CustomAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for custom acceleration.
	
	Class to provide settings for custom accelerations. This is done by means of a function and, if necessary,
	an associated scaling function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class DirectTidalDissipationAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for direct tidal dissipation acceleration.
	
	Class to provide settings for direct tidal dissipation accelerations. Creates settings for tidal accelerations.
	The direct of tidal effects in a satellite system is applied directly as an acceleration
	(as opposed to a modification of spherical harmonic coefficients).
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class EmpiricalAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for the empirical acceleration.
	
	Class to provide settings for empirical accelerations. These are expressed in the
	RSW frame, for which the magnitude is determined empirically (typically during an orbit determination process).
	The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
	a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
	the RSW frame.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class MomentumWheelDesaturationAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for momentum wheel desaturation acceleration.
	
	Class to provide settings for momentum wheel desaturation acceleration. Settings for the direction and magnitude
	of the thrust are included.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class MutualSphericalHarmonicAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for the mutual spherical harmonic acceleration.
	
	Class for providing settings for the mutual spherical harmonics acceleration model,
	including the maximum degree and order up to which the fields of the bodies are to be expanded. Note that
	the minimum degree and order are currently always set to zero.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RelativisticAccelerationCorrectionSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for the relativistic acceleration correction.
	
	Class to provide settings for typical relativistic corrections to the dynamics of an orbiter: the
	Schwarzschild, Lense-Thirring and de Sitter terms (see 'General relativity and Space Geodesy' by L. Combrinck,
	2012).
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class SphericalHarmonicAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for the spherical harmonic acceleration.
	
	Class for providing settings for spherical harmonics acceleration model,
	including the maximum degree and order up to which the field is to be expanded. Note that
	the minimum degree and order are currently always set to zero.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ThrustAccelerationSettings(AccelerationSettings):
    """`AccelerationSettings`-derived class to define settings for thrust acceleration, listing the engine models that are to be used
	
	Class to provide settings for thrust acceleration, listing the engine models that are to be used
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def direction_settings(self) -> tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust.ThrustDirectionSettings:
        ...

    @property
    def magnitude_settings(self) -> tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust.ThrustMagnitudeSettings:
        ...

def aerodynamic() -> AccelerationSettings:
    """Creates settings for the aerodynamic acceleration.
	
	Creates settings for the aerodynamic acceleration. The acceleration is computed from:
	
	.. math::
	   \\mathbf{a}=-
		ext{Aero})}\\left(
	
	
	
	
	with :math:`\\mathbf{R}^{(I/	 ext{Aero})}` the rotation matrix from the aerodynamic frame of the body undergoing acceleration to the inertial frame (computed from the body's current state, and the rotation of the body exerting the acceleration), :math:`
	
	The body exerting the acceleration needs to have an
	atmosphere (:ref:`\\`\\`gravity_field\\`\\`` module), shape (:ref:`\\`\\`shape\\`\\`` module) and rotation model (:ref:`\\`\\`rotation_model\\`\\`` module) defined. The body undergoing the acceleration needs to have aerodynamic coefficients (:ref:`\\`\\`aerodynamic_coefficients\\`\\`` module) defined.
	
	:return:
			Acceleration settings object.
	"""

def cannonball_radiation_pressure() -> AccelerationSettings:
    """Creates settings for the cannonball radiation pressure acceleration.
	
	Creates settings for the radiation pressure acceleration, for which a cannonball model is used. The acceleration is computed from:
	
	.. math::
	
	   \\mathbf{a}=\\left(
	
	
	
	
	
	
	with :math:`P` the total emitted radiation power for the body exerting the acceleration, :math:`C_{r}` the radiation pressure coefficient with reference area :math:`S_{ref}`, :math:`\\mathbf{r}` the vector from the body exerting the acceleration to the body undergoing the acceleration, and :math:`m` the mass of the body undergoing acceleration
	
	In this model,
	the effective acceleration is colinear with the vector connecting the source of radiation and the target.
	The body undergoing the acceleration needs to have a radiation pressure model defined, while the body emitting
	radiation needs to have radiative properties defined.
	
	:return:
			Acceleration settings object.
	"""

def custom(acceleration_function: typing.Callable[[float], numpy.ndarray]) -> AccelerationSettings:
    ...

def custom_acceleration(acceleration_function: typing.Callable[[float], numpy.ndarray]) -> AccelerationSettings:
    """Creates settings for custom acceleration.
	
	Creates settings for a custom accelerations, this acceleration must be parameterized as a function of time,
	and expressed with an inertial orientation.
	
	
	:param acceleration_function:
			Custom acceleration function with time as an independent variable, returning the acceleration in an inertial frame (*e.g.* with global frame orientation) as a function of time.
	:return:
			Custom acceleration settings object.
	"""

@typing.overload
def direct_tidal_dissipation_acceleration(k2_love_number: float, time_lag: float, include_direct_radial_component: bool=True, use_tide_raised_on_planet: bool=True, explicit_libraional_tide_on_satellite: bool=False) -> AccelerationSettings:
    """Creates settings for custom acceleration.
	
	Creates settings for tidal accelerations. The direct of tidal effects in a satellite system is applied directly as
	an acceleration (as opposed to a modification of spherical harmonic coefficients).
	The model is based on Lainey et al. (2007, 2012). It can compute the acceleration due to tides, and in
	particular tidal dissipation, on a planetary satellite. The acceleration computed can account for either the
	effect of tide raised on the satellite by the planet or on the planet by the satellite. The satellite is assumed
	to be tidally locked to the planet.
	
	
	:param k2_love_number:
			Value of the k2 Love number.
	:param time_lag:
			Value of the tidal time lag.
	:param include_direct_radial_component:
			It denotes whether the term independent of the time lag is to be computed.
	:param use_tide_raised_on_planet:
			It denotes whether the tide raised on the planet is to be modelled (if true) or the tide raised on the satellite (if false).
	:return:
			Direct tidal dissipation acceleration settings object.
	"""

@typing.overload
def direct_tidal_dissipation_acceleration(k2_love_number: float, inverse_tidal_quality_factor: float, tidal_period: float, include_direct_radial_component: bool=True, use_tide_raised_on_planet: bool=True, explicit_libraional_tide_on_satellite: bool=False) -> AccelerationSettings:
    """Creates settings for custom acceleration.
	
	Creates settings for tidal accelerations. The direct of tidal effects in a satellite system is applied directly as
	an acceleration (as opposed to a modification of spherical harmonic coefficients).
	The model is based on Lainey et al. (2007, 2012). It can compute the acceleration due to tides, and in
	particular tidal dissipation, on a planetary satellite. The acceleration computed can account for either the
	effect of tide raised on the satellite by the planet or on the planet by the satellite. The satellite is assumed
	to be tidally locked to the planet.
	
	
	:param k2_love_number:
			Value of the k2 Love number.
	:param time_lag:
			Value of the tidal time lag.
	:param include_direct_radial_component:
			It denotes whether the term independent of the time lag is to be computed.
	:param use_tide_raised_on_planet:
			It denotes whether the tide raised on the planet is to be modelled (if true) or the tide raised on the satellite (if false).
	:return:
			Direct tidal dissipation acceleration settings object.
	"""

def empirical(constant_acceleration: numpy.ndarray=..., sine_acceleration: numpy.ndarray=..., cosine_acceleration: numpy.ndarray=...) -> AccelerationSettings:
    """Creates settings for empirical acceleration.
	
	Creates settings for empirical accelerations. These are expressed in the
	RSW frame, for which the magnitude is determined empirically (typically during an orbit determination process).
	The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
	a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
	the RSW frame. The empirical acceleration is calculated as:
	
	 .. math::
	
		\\mathbf{a}=R^{I/RSW}\\left(\\mathbf{a}_{	  ext{const.}}+\\mathbf{a}_{\\sin}\\sin	  heta+\\mathbf{a}_{\\cos}\\cos	  heta
	
	
	Here, :math:`R^{I/RSW}` is the rotation matrix from the RSW frame (of the body undergoing the acceleration w.r.t. the
	body exerting the acceleration), :math:`theta` is the true anomaly, and the three constituent acceleration vectors are
	the inputs provided in the above code block. The body 'exerting' the acceleration is considered to be the
	central body, w.r.t. which the true anomaly is calculated.
	
	
	:param constant_acceleration:
			Constant term, defined in the RSW frame.
	:param sine_acceleration:
			Sine term (function of the true anomaly), defined in the RSW frame..
	:param cosine_acceleration:
			Cosine term (function of the true anomaly), defined in the RSW frame..
	:return:
			Empirical acceleration settings object.
	"""

def mutual_spherical_harmonic_gravity(maximum_degree_body_exerting: int, maximum_order_body_exerting: int, maximum_degree_body_undergoing: int, maximum_order_body_undergoing: int, maximum_degree_central_body: int=0, maximum_order_central_body: int=0) -> AccelerationSettings:
    """Creates settings for the mutual spherical harmonic gravity acceleration.
			
				Creates settings for the mutual spherical harmonic gravity acceleration. This model computes the total spherical harmonic acceleration exerted by a body :math:`B` on a body :math:`A`, where the influence of the gravity field coefficients of body :math:`A` itself has been included. The model includes couplings between the mass of each body, and the gravity field coefficients of the other body. It does not include the 'figure-figure' interactions (coupling between the two-bodies' gravity field coefficients). It corresponds to the model presented by Lainey et al. (2004); Dirkx et al. (2016).
				The model combines the spherical harmonic accelerations of the two bodies (see :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic`) on each other. The direct acceleration (acceleration w.r.t. an inertial origin) is computed from:
			
				.. math::
			
				   \\mathbf{a}={-
		rac{\\mu_{_{B}}}{{r}^{2}}\\hat{\\mathbf{r}}}+{\\mathbf{R}^{(I/B)}
		abla^{(B)}U_{\\hat{B}}(\\mathbf{r})}-{
		rac{\\mu_{_{B}}}{\\mu_{_{A}}}\\mathbf{R}^{(I/A)}
		abla^{(A)}U_{\\hat{A}}(-\\mathbf{r})}
			
				where :math:`U_{\\hat{B}}` and :math:`U_{\\hat{A}}` denote the spherical harmonic gravity fields a degree :math:`>=1` of bodies :math:`B` and :math:`A`, respectively.
				Both the body exerting the acceleration and the body undergoing it need to
				have spherical harmonic gravity field and rotation models defined.
			
				Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct spherical harmonic attraction (as above), a central spherical harmonic attraction (:math:`\\mu_{B}
		ightarrow\\mu_{B}+\\mu_{A}`, in the above equation and in :math:`U_{\\hat{B}}`) or a third-body spherical harmonic attraction (see `here <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/propagation_setup/acceleration_models/third_body_acceleration.html>`_ for more details).
			
				For the case where a third-body mutual spherical harmonic acceleration,
				additional parameters have to be provided that denote the expansion degree/order of the central body (``maximum_degree_central_body`` and ``maximum_order_central_body``)
			
			
				:param maximum_degree_body_exerting:
						Maximum degree of the spherical harmonic expansion for the body exerting the acceleration.
				:param maximum_order_body_exerting:
						Maximum order of the spherical harmonic expansion for the body exerting the acceleration.
				:param maximum_degree_body_undergoing:
						Maximum degree of the spherical harmonic expansion for the body undergoing the acceleration.
				:param maximum_order_body_undergoing:
						Maximum order of the spherical harmonic expansion for the body undergoing the acceleration.
				:param maximum_degree_central_body:
						Maximum degree of the spherical harmonic expansion for the central body, if needed.
				:param maximum_order_central_body:
						Maximum order of the spherical harmonic expansion for the central body, if needed.
				:return:
						Spherical harmonic acceleration settings object.
			
		"""

def point_mass_gravity() -> AccelerationSettings:
    """Creates settings for the point-mass gravity acceleration.
	
	Creates settings for the point-mass gravity acceleration. The direct acceleration (acceleration w.r.t. an inertial frame) is computed from:
	
	.. math::
	   \\mathbf{a}=
	
	
	with :math:`\\mathbf{r}` the position vector measured from the center of mass of the body exerting the acceleration.
	
	The body exerting the acceleration needs to have a gravity field model (:ref:`\\`\\`gravity_field\\`\\`` module) defined to use this acceleration.
	
	Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation \\math:`C`, choosing this option may create a direct point-mass attraction (:math:`\\mu=\\mu_{B}`), a central point-mass attraction (:math:`\\mu=\\mu_{B}+\\mu_{A}`) or a third-body point-mass attraction (see `here <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/propagation_setup/acceleration_models/third_body_acceleration.html>`_ for more details).
	
	:return:
			Acceleration settings object.
	"""

def polyhedron_gravity() -> AccelerationSettings:
    ...

def quasi_impulsive_shots_acceleration(thrust_mid_times: list[float], delta_v_values: list[numpy.ndarray], total_maneuver_time: float, maneuver_rise_time: float) -> AccelerationSettings:
    """Creates settings for incorporating quasi-impulsive shots into the acceleration.
	
	The acceleration model is purpose-built to represent short bursts of thrust, such as a momentum wheel desaturation.
	A typical use case is precise orbit determination, but the functionality can be used just as well in propagation
	(for instance to model an impulsive manuever in a continuous manner when going from preliminary modelling to
	full modelling). The thrust is modelled similarly to Fig. 3 of Alessi et al. (2012), with the main difference
	being that a third-order polynomial to go from zero acceleration to the maximum acceleration level is employed.
	By using a 3rd-order polynomial and imposing continuity in the value and first derivative of the acceleration,
	defining the rise time (time it takes acceleration to go from 0 to its maximum level), the total time where
	there is non-zero thrust (total maneuver time), and the total Delta V exerted by a single maneuver,
	the acceleration profile is fully defined.
	
	
	:param thrust_mid_times:
			Set of middle point in times in the maneuver denoting the epoch of each maneuver.
	:param delta_v_values:
			Set of delta V, one for each maneuver.
	:param total_maneuver_time:
			Total duration of every maneuver.
	:param maneuver_rise_time:
			Time taken by the acceleration to go from zero to its maximum level.
	:return:
			Momentum wheel desaturation acceleration settings object.
	"""

def radiation_pressure(target_type: tudatpy.numerical_simulation.environment_setup.radiation_pressure.expose_radiation_pressure.RadiationPressureTargetModelType=...) -> AccelerationSettings:
    ...

def relativistic_correction(use_schwarzschild: bool=False, use_lense_thirring: bool=False, use_de_sitter: bool=False, de_sitter_central_body: str='', lense_thirring_angular_momentum: numpy.ndarray=...) -> AccelerationSettings:
    """Creates settings for the relativistic acceleration correction.
	
	Creates settings for typical relativistic acceleration corrections: the Schwarzschild, Lense-Thirring and de
	Sitter terms, where each of the three terms can be toggled on or of (see 'General relativity and Space Geodesy' by L. Combrinck, 2012). It implements the model of
	2010 Conventions (chapter 10, section 3). Here, the ‘primary body’ for a planetary orbiter should always be set
	as the Sun (only relevant for de Sitter correction). The angular momentum vector of the orbited body is only
	relevant for Lense-Thirring correction.
	
	
	:param use_schwarzschild:
			Boolean defining whether or not to use the Schwarzschild contribution to the acceleration correction
	:param use_lense_thirring:
			Boolean defining whether or not to use the Lense-Thirring contribution to the acceleration correction
	:param use_de_sitter:
			Boolean defining whether or not to use the de Sitter contribution to the acceleration correction
	:param de_sitter_central_body:
			Body used as 'perturbed' in the calculation of the de Sitter acceleration. For the case of an Earth-orbiting satellite, this would be the Sun
	:param lense_thirring_angular_momentum:
			Angular momentum vector (in global frame) that is to be used for the calculation of the Lense-Thirring acceleration
	:return:
			Relativistic acceleration correction settings object.
	"""

def ring_gravity() -> AccelerationSettings:
    ...

def spherical_harmonic_gravity(maximum_degree: int, maximum_order: int) -> AccelerationSettings:
    """Creates settings for the spherical harmonic gravity acceleration.
			
				Creates settings for the spherical harmonic gravity acceleration, accounting for a finite (given as input) number
				of degrees and orders. The direct acceleration (acceleration w.r.t. an inertial origin) is computed from:
			
				.. math::
				   \\mathbf{a}=\\mathbf{R}^{(I/B)}
		abla^{(B)}U(\\mathbf{r})
			
				with :math:`\\mathbf{r}` the position vector measured from the center of mass of the body exerting the acceleration, :math:`\\mathbf{R}^{(I/B)}` the rotation matrix from body-fixed to inertial frame, and :math:`
		abla^{(B)}` the gradient operator in a body-fixed frame, and :math:`U` the spherical harmonic gravitational potential, expanded up to the provided ``maximum_degree`` and ``maximum_order``.
			
				The body exerting the acceleration needs to have a spherical harmonic gravity field model (see :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`) and a rotation model (:ref:`\\`\\`rotation_model\\`\\`` module) defined.
			
				Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct spherical harmonic attraction (:math:`\\mu=\\mu_{B}`), a central spherical harmonic attraction (:math:`\\mu=\\mu_{B}+\\mu_{A}`) or a third-body spherical harmonic attraction (see `here <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/propagation_setup/acceleration_models/third_body_acceleration.html>`_ for more details).
			
			
				:param maximum_degree:
						Maximum degree of the spherical harmonic expansion.
				:param maximum_order:
						Maximum order of the spherical harmonic expansion.
				:return:
						Spherical harmonic acceleration settings object.
			
		"""

def thrust_and_isp_from_custom_function(thrust_force_function: typing.Callable[[float], numpy.ndarray], constant_specific_impulse: float, thrust_frame: tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust.ThrustFrames=..., central_body: str='') -> AccelerationSettings:
    ...

def thrust_from_all_engines() -> AccelerationSettings:
    """Creates settings for thrust acceleration using a single engine models.
	
	Creates settings for thrust acceleration by combining thryst from all engines defined in the body.. See the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/propagation_setup/acceleration_models/thrust.html>`_
	for more details on the definition of a thrust model in Tudat.
	
	:return:
			Thrust acceleration settings object.
	"""

def thrust_from_custom_function(thrust_force_function: typing.Callable[[float], numpy.ndarray], specific_impulse_function: typing.Callable[[float], float], thrust_frame: tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust.ThrustFrames=..., central_body: str='') -> AccelerationSettings:
    ...

def thrust_from_direction_and_magnitude(thrust_direction_settings: tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust.ThrustDirectionSettings, thrust_magnitude_settings: tudatpy.numerical_simulation.propagation_setup.thrust.expose_thrust.ThrustMagnitudeSettings) -> AccelerationSettings:
    ...

def thrust_from_engine(engine_name: str) -> AccelerationSettings:
    """Creates settings for thrust acceleration using a single engine models.
	
	Creates settings for thrust acceleration using a single engine models. See the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/propagation_setup/acceleration_models/thrust.html>`_
	for more details on the definition of a thrust model in Tudat.
	
	
	:param engine_name:
			Name of engine to use when computing thrust.
	:return:
			Thrust acceleration settings object.
	"""

def thrust_from_engines(engine_names: list[str]) -> AccelerationSettings:
    """Creates settings for thrust acceleration using a list of engine models.
	
	Creates settings for thrust acceleration using a list of engine models. See the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/propagation_setup/acceleration_models/thrust.html>`_
	for more details on the definition of a thrust model in Tudat.
	
	
	:param engine_names:
			List of engine names to use when computing thrust.
	:return:
			Thrust acceleration settings object.
	"""

def yarkovsky(yarkovsky_parameter: float) -> AccelerationSettings:
    ...
aerodynamic_type: AvailableAcceleration
cannonball_radiation_pressure_type: AvailableAcceleration
custom_acceleration_type: AvailableAcceleration
direct_tidal_dissipation_in_central_body_acceleration_type: AvailableAcceleration
direct_tidal_dissipation_in_orbiting_body_acceleration_type: AvailableAcceleration
empirical_acceleration_type: AvailableAcceleration
mutual_spherical_harmonic_gravity_type: AvailableAcceleration
point_mass_gravity_type: AvailableAcceleration
polyhedron_gravity_type: AvailableAcceleration
quasi_impulsive_shots_acceleration_type: AvailableAcceleration
radiation_pressure_type: AvailableAcceleration
relativistic_correction_acceleration_type: AvailableAcceleration
ring_gravity_type: AvailableAcceleration
spherical_harmonic_gravity_type: AvailableAcceleration
thrust_acceleration_type: AvailableAcceleration
undefined_acceleration_type: AvailableAcceleration