import numpy
import tudatpy.math.interpolators.expose_interpolators
import typing
__all__ = ['IAUConventions', 'PlanetaryRotationModelSettings', 'RotationModelSettings', 'RotationModelType', 'SimpleRotationModelSettings', 'aerodynamic_angle_based', 'constant_rotation_model', 'custom_inertial_direction_based', 'custom_rotation_model', 'gcrs_to_itrs', 'gcrs_to_itrs_rotation_model', 'iau_2000_a', 'iau_2000_b', 'iau_2006', 'mars_high_accuracy', 'orbital_state_direction_based', 'planetary_rotation_model', 'simple', 'simple_from_spice', 'simple_rotational_model', 'spice', 'spice_rotation_model', 'synchronous', 'synchronous_rotation_model', 'zero_pitch_moment_aerodynamic_angle_based']

class IAUConventions:
    """Enumeration of IAU conventions for Earth rotation.
	
	Enumeration of IAU conventions for Earth rotation supported by tudat.
	
	
	:member iau_2000_a:
	:member iau_2000_b:
	:member iau_2006:
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, IAUConventions]]
    iau_2000_a: typing.ClassVar[IAUConventions]
    iau_2000_b: typing.ClassVar[IAUConventions]
    iau_2006: typing.ClassVar[IAUConventions]

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

class PlanetaryRotationModelSettings(RotationModelSettings):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RotationModelSettings:
    """Base class for providing settings for automatic rotation model creation.
	
	This class is a functional base class for settings of rotation models that require no information in addition to their type.
	Basic rotation model has constant orientation of the rotation axis (body-fixed z-axis) and constant rotation rate about this axis.
	Rotation models requiring additional information must be created using the factory functions which create the specific object derived from this base class.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def base_frame(self) -> str:
        """
        Name of the base frame of rotation model.
        	
        """

    @base_frame.setter
    def base_frame(self, arg1: str) -> None:
        ...

    @property
    def rotation_type(self) -> RotationModelType:
        """
        Type of rotation model that is to be created.
        	
        """

    @property
    def target_frame(self) -> str:
        """
        Name of the target frame of rotation model.
        	
        """

class RotationModelType:
    """Enumeration of rotation model types.
	
	Enumeration of rotation model types supported by tudat.
	
	
	:member simple_rotation_model:
	:member spice_rotation_model:
	:member gcrs_to_itrs_rotation_model:
	:member synchronous_rotation_model:
	:member planetary_rotation_model:
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, RotationModelType]]
    gcrs_to_itrs_rotation_model: typing.ClassVar[RotationModelType]
    planetary_rotation_model: typing.ClassVar[RotationModelType]
    simple_rotational_model: typing.ClassVar[RotationModelType]
    spice_rotation_model: typing.ClassVar[RotationModelType]
    synchronous_rotation_model: typing.ClassVar[RotationModelType]

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

class SimpleRotationModelSettings(RotationModelSettings):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def aerodynamic_angle_based(central_body: str, base_frame: str, target_frame: str, angle_funcion: typing.Callable[[float], numpy.ndarray]=None) -> RotationModelSettings:
    """Factory function for creating rotation model settings based on custom aerodynamic angles (attack, sideslip, bank).
	
	Factory function for creating rotation model settings based on custom aerodynamic angles:
	angle of attack :math:`\x07lpha`, sideslip angle :math:`\x08eta` and bank angle :math:`\\sigma`. The use of this function is typical for
	simulating the dynamics of a (guided) re-entry vehicle. It calculates the rotation matrix from inertial frame to the body-fixed frame
	of the current body B (typically a vehicle) w.r.t. the body-fixed frame of a central body C (e.g., the body at which the re-entry is taking place.
	The full algorithm for :math:`R^{(I/B)}` is described by Mooij (1994), and is composed of:
	
	*  The rotation from inertial frame to the body fixed frame of body C, using the existing rotation model of body C
	*  The rotation from body-fixed frame of body C to the vehicle's vertical frame V. This rotation uses the current latitude and longitude angles.
	*  The rotation of the vehicle's vertical frame V to its trajectory frame T. This rotation uses the current heading and flight path angles.
	*  The rotation of the vehicle's trajectory frame T to its aerodynamic frame A. This rotation uses the current bank angle
	*  The rotation of the vehicle's aerodynamic frame A to its body-fixed frame. This rotation uses the current angle of attack and sideslip angles
	
	In the above algorithm, the latitude, longitude, heading and flight-path angles are computed from the vehicle's current translational state, in the body-fixed
	frame of body C. The angle of attack, sideslip angle and bank angle are to be defined by the user, through a single custom function that is passed to
	the ``angle_function`` argument of this functions
	
	
	:param central_body:
			Name of the central body C that is to be used.
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:param angle_function:
			Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\x07lpha,\x08eta,\\sigma]`. If this input is left empty, these angles are both fixed to 0.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.
	"""

def constant_rotation_model(base_frame: str, target_frame: str, initial_orientation: numpy.ndarray) -> RotationModelSettings:
    """Factory function for creating simple rotation model settings for target-frames with constant orientation.
	
	Factory function for settings object, defining simple rotation model setting objects with constant rotation matrix.
	These model settings are for target frames which do not have a rotational rate in the base frame and are fully defined by their initial orientation.
	
	
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:param initial_orientation:
			Rotation matrix from inertial to body-fixed (base to target) frame at initial time (constant throughout).
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class.
	"""

def custom_inertial_direction_based(inertial_body_axis_direction: typing.Callable[[float], numpy.ndarray], base_frame: str, target_frame: str, free_rotation_angle_function: typing.Callable[[float], float]=None) -> RotationModelSettings:
    """Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction
	
	Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction :math:`\\hat{\\mathbf{T}}_{I}`. Specifically, it ensures
	that the rotation matrix from body-fixed to inertial frame is set up such that :math:`\\hat{\\mathbf{T}}_{I}=R^{(I/B)}\\hat{\\mathbf{i}}` (where :math:`\\mathbf{i}` is the unit-vector in local x-direction).
	The complete rotation matrix requires an additional angle :math:`\\phi` (rotation of the body about its body-fixed x-axis), which is set to 0 by default.
	
	The full rotation matrix is computed from a 3-2-1 Euler angle rotation
	:math:`R^{(I/B)}=R_{z}(\\psi)R_{y}(	  heta)R_{x}(\\phi)`, with :math:`\\psi` and :math:`		heta` computed from the suitable decomposition of :math:`\\hat{\\mathbf{T}}_{I}`.
	This function is typically used for simulating the (guided) dynamics of a spacecraft under thrust, where the thrust is provided in the x-direction of the body-fixed frame. By providing a suitable
	``inertial_body_axis_direction``, this thrust can be defined to point in an arbitrary direction (typically defined by a guidance algorithm) in the inertial frame as a function of time.
	
	NOTE: this function may be extended in the future to allow an arbitrary body-fixed direction to align with an arbitrary inertial direction. At present, its functionality is limited to imposing the inertial direction of the body-fixed x-axis.
	
	
	:param inertial_body_axis_direction:
			Custom function defined by the user, which imposes the inertial orientation of the body-fixed x-axis, by providing :math:`\\hat{\\mathbf{T}}_{I}(t)`.
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:param free_rotation_angle_function:
			Custom function provided by the user, which returns a value for the free rotation angle :math:`\\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model.
	"""

def custom_rotation_model(base_frame: str, target_frame: str, custom_rotation_matrix_function: typing.Callable[[float], numpy.ndarray], finite_difference_time_step: float) -> RotationModelSettings:
    ...

def gcrs_to_itrs(precession_nutation_theory: IAUConventions=..., base_frame: str='GCRS', cio_interpolation_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorGenerationSettings=None, tdb_to_tt_interpolation_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorGenerationSettings=None, short_term_eop_interpolation_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorGenerationSettings=None) -> RotationModelSettings:
    """Factory function for creating high-accuracy Earth rotation model settings.
	
	Factory function for settings object, defining high-accuracy Earth rotation model according to the IERS 2010 Conventions.
	This settings class has various options to deviate from the default settings, typical applications will use default.
	Note that for this model the original frame must be J2000 or GCRS (in the case of the former, the frame bias between GCRS and J2000 is automatically corrected for). The target frame (e.g. body-fixed frame) name is ITRS.
	The precession-nutation theory may be any member of :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.IAUConventions` (``iau_2000a`` / ``iau_2000b`` or ``iau_2006``).
	Alternative options to modify the input (not shown here) include the EOP correction file, input time scale, short period UT1 and polar motion variations.
	The target frame (e.g. body-fixed frame) name is ITRS.
	
	
	:param precession_nutation_theory:
			Setting theory for modelling Earth nutation.
	:param base_frame:
			Base frame of rotation model
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.GcrsToItrsRotationModelSettings` class
	"""

def mars_high_accuracy(base_frame: str='ECLIPJ2000', target_frame: str='Mars_Fixed') -> RotationModelSettings:
    """Factory function for creating a high-accuracy Mars rotation model.
	
	Factory function for creating a high-accuracy Mars rotation model, using the default parameters of `Konopliv et al. (2016) <https://www.sciencedirect.com/science/article/abs/pii/S0019103516001305>`_
	and the mathematical model of ` Konopliv et al. (2006) <https://www.sciencedirect.com/science/article/pii/S0019103506000297>`_. The rotation matrix formulation is given in Eq. (13)-(19) of that paper.
	Note that, at the moment, all parameters in this rotation model are hard coded, and cannot be adapted by the user (except by estimating a number of its constituent parameters, see :ref:`\\`\\`parameter\\`\\`` module )
	As such, this model is at present applicable to Mars rotation only. If you require more fine-grained control of the parameters, please contact the Tudat support team
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.PlanetaryRotationModelSettings` class, which defines the required settings for the rotation model.
	"""

def orbital_state_direction_based(central_body: str, is_colinear_with_velocity: bool, direction_is_opposite_to_vector: bool, base_frame: str, target_frame: str='', free_rotation_angle_function: typing.Callable[[float], float]=None) -> RotationModelSettings:
    """Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector.
	
	Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector. This function is
	similar to the :func:`~custom_inertial_direction_based` function, with the exception that the :math:`\\hat{\\mathbf{T}}_{I}` vector is not defined by thee user, but is defined by the
	relative position vector :math:`\\mathbf{r}_{C}` or velocity vector :math:`\\mathbf{r}_{C}` of the vehicle w.r.t. some body C. The inputs to this function allow :math:`\\hat{\\mathbf{T}}_{I}` to
	be set to :math:`\\pm\\mathbf{r}_{C}` or :math:`\\pm\\mathbf{v}_{C}`, for any body C. It is typically used for simplified or preliminary thrust analyses.
	
	
	:param central_body:
			Name of central body w.r.t. which the position/velocity vector is to be computed
	:param is_colinear_with_velocity:
			Boolean defining whether :math:`\\hat{\\mathbf{T}}_{I}` is to be aligned with velocity (if true) or position (if false)
	:param direction_is_opposite_to_vector:
			Boolean defining whether :math:`\\hat{\\mathbf{T}}_{I}` is to be in the same direction as position/velocity (if false), or in the opposite direction (if true).
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:param free_rotation_angle_function:
			Custom function provided by the user, which returns a value for the free rotation angle :math:`\\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model.
	"""

def simple(base_frame: str, target_frame: str, initial_orientation: numpy.ndarray, initial_time: float, rotation_rate: float) -> RotationModelSettings:
    """Factory function for creating simple rotation model settings.
	
	Factory function for settings object, defining a basic rotation model with constant orientation of the rotation axis and constant rotation rate about this axis.
	Rotation from original (inertial) to target (body-fixed) frame at some reference time ``initial_time`` (:math:`t_{0}`) is defined by the ``initial_orientation`` (:math:`\\mathbf{R}^{(B/I)}(t_{0})`) rotation matrix.
	Rotation about the body-fixed z-axis is defined by the ``rotation_rate`` (:math:`\\omega`) float variable (in rad/s). The rotation matrix is computed from:
	
	.. math::
	   \\mathbf{R}^{(B/I)}(t)=\\mathbf{R}_{z}(\\omega(t-t_{0}))(t_{0})\\mathbf{R}^{(B/I)}(t_{0})
	
	where :math:`\\mathbf{R}^{(B/I)}` denotes the rotation matrix from inertial to body-fixed frame, and :math:`\\mathbf{R}_{z}` denotes a rotaion matrix about the z-axis.
	
	The matrix :math:`\\mathbf{R}^{(B/I)}(t_{0})` is sometimes parameterized by pole right ascension and declination (:math:`\x07lpha` and :math:`\\delta`), as well as the meridian of date :math:`W_{0}` with
	
	.. math::
	   \\mathbf{R}^{(B/I)}(t_{0})=\\mathbf{R}_{z}(W_{0})\\mathbf{R}_{x}(\\pi/2-\\delta)\\mathbf{R}_{z}(\\pi/2+\x07lpha)
	
	
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:param initial_orientation:
			Orientation of target frame in base frame at initial time.
	:param initial_time:
			Initial time (reference epoch for rotation matrices).
	:param rotation_rate:
			Constant rotation rate [rad/s] about rotational axis.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class
	"""

def simple_from_spice(base_frame: str, target_frame: str, target_frame_spice: str, initial_time: float) -> RotationModelSettings:
    """Factory function for creating simple rotation model settings using initial orientation and rotation rates from Spice.
	
	Factory function for settings object, defining a :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` rotation model with the added functionality that the initial orientation and rotation rate are extracted from Spice, as opposed to provided manually.
	Note that `only` the initial orientation and rotation rate ( at the time defined by ``initial_time`` ) are extracted from Spice - for
	the full Spice rotation model see :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.spice`.
	Also note the distinction between the ``target_frame`` and ``target_frame_spice`` parameters.
	
	
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Target frame of rotation model - name of frame that Tudat assigns to the body-fixed frame
	:param target_frame_spice:
			Spice reference of target frame - name of the frame in Spice for which the initial orientation and rotation rate are extracted.
	:param initial_time:
			Initial time (reference epoch for rotation matrices).
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class
	"""

def spice(base_frame: str, target_frame: str, spice_frame_name: str='') -> RotationModelSettings:
    """Factory function for creating rotation model settings from the Spice interface.
	
	Factory function for settings object, defining a rotation model directly (and entirely) from Spice interface.
	
	
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` class.
	"""

def synchronous(central_body_name: str, base_frame: str, target_frame: str) -> RotationModelSettings:
    """Factory function for creating synchronous rotational ephemeris settings.
	
	Factory function for settings object, defining a synchronous rotation model where rotation of a body is defined from its relative orbit w.r.t. some central body. Specifically
	- the body-fixed x-axis is *always* pointing towards the central body
	- the body-fixed z-axis is *always* perpendicular to the orbital plane (along the direction of :math:`\\mathbf{x}		imes\\mathbf{v}` )
	- the body-fixed y-axis completes the right-handed reference frame
	
	Such a model can be useful for, for instance, approximate rotation of tidally locked natural satellites or nadir-pointing spacecraft.
	
	
	:param central_body_name:
			Name of the base frame of rotation model.
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Spice reference of target frame.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SynchronousRotationModelSettings` class
	"""

def zero_pitch_moment_aerodynamic_angle_based(central_body: str, base_frame: str, target_frame: str, angle_funcion: typing.Callable[[float], numpy.ndarray]=None) -> RotationModelSettings:
    """Factory function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank.
	
	Factory function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank. This function is
	largely identical to the :func:`~aerodynamic_angle_based`, with the difference that the angle of attack :math:`\x07lpha` is not provided as a custom value by the user, but is
	calculated from the body's aerodynamic moment coefficients, such that we have :math:`C_{m}=0`. This requires aerodynamic moment coefficients to be defined for the vehicle that
	depend on (among others) the body's angle of attack
	
	
	:param central_body:
			Name of the central body C that is to be used.
	:param base_frame:
			Name of the base frame of rotation model.
	:param target_frame:
			Name of the target frame of rotation model.
	:param angle_funcion:
			Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\x08eta,\\sigma]`. If this input is left empty, these angles are both fixed to 0.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.
	"""
gcrs_to_itrs_rotation_model: RotationModelType
iau_2000_a: IAUConventions
iau_2000_b: IAUConventions
iau_2006: IAUConventions
planetary_rotation_model: RotationModelType
simple_rotational_model: RotationModelType
spice_rotation_model: RotationModelType
synchronous_rotation_model: RotationModelType