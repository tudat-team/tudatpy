import numpy
import tudatpy.math.interpolators.expose_interpolators
import typing
__all__ = ['ApproximateJplEphemerisSettings', 'ConstantEphemerisSettings', 'CustomEphemerisSettings', 'DirectSpiceEphemerisSettings', 'EphemerisSettings', 'InterpolatedSpiceEphemerisSettings', 'KeplerEphemerisSettings', 'ScaledEphemerisSettings', 'TabulatedEphemerisSettings', 'approximate_jpl_model', 'constant', 'create_ephemeris', 'custom', 'custom_ephemeris', 'direct_spice', 'interpolated_spice', 'keplerian', 'keplerian_from_spice', 'scaled_by_constant', 'scaled_by_vector', 'scaled_by_vector_function', 'tabulated', 'tabulated_from_existing']

class ApproximateJplEphemerisSettings(EphemerisSettings):
    """Class for creating settings of approximate ephemeris for major planets.
	
	`EphemerisSettings` derived class for approximate ephemeris for major planets as implemented in ApproximateJplEphemerisSettings class and derived class (described in `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).
	"""

    @property
    def body_name(self) -> str:
        ...

class ConstantEphemerisSettings(EphemerisSettings):
    """Class for defining settings of constant ephemerides.
	
	`EphemerisSettings` derived class for ephemerides producing a constant (time-independent) state.
	"""

class CustomEphemerisSettings(EphemerisSettings):
    """Class for defining settings of a custom ephemeris.
	
	`EphemerisSettings` derived class for ephemerides which represent an ideal Kepler orbit.
	"""

    @property
    def get_custom_state_function(self) -> typing.Callable[[float], numpy.ndarray]:
        ...

class DirectSpiceEphemerisSettings(EphemerisSettings):
    """Class for defining settings of an ephemeris linked directly to Spice.
	
	`EphemerisSettings` derived class for ephemeris which are directly linked to Spice.
	"""

    @property
    def converge_light_time_aberration(self) -> bool:
        """
        Boolean defining whether to use single iteration or max. 3 iterations for calculating light time correction.
        	
        """

    @property
    def correct_for_light_time_aberration(self) -> bool:
        """
        Boolean defining whether to correct for light time in retrieved values (of observed state).
        	
        """

    @property
    def correct_for_stellar_aberration(self) -> bool:
        """
        Boolean defining whether to correct for stellar aberrations in retrieved values (of observed state).
        	
        """

class EphemerisSettings:
    """Base class for providing settings for ephemeris model.
	
	Functional (base) class for settings of ephemeris models that require no information in addition to their type (and frame origin and orientation).
	Ephemeris model classes requiring additional information must be created using an object derived from this class.
	"""

    @property
    def ephemeris_type(self) -> ...:
        """
        Type of ephemeris that is to be created.
        	
        """

    @property
    def frame_orientation(self) -> str:
        """
        Orientation of frame in which ephemeris data is to be defined.
        	
        """

    @frame_orientation.setter
    def frame_orientation(self, arg1: str) -> None:
        ...

    @property
    def frame_origin(self) -> str:
        """
        Origin of frame in which ephemeris data is to be defined.
        	
        """

    @frame_origin.setter
    def frame_origin(self, arg1: str) -> None:
        ...

    @property
    def make_multi_arc_ephemeris(self) -> bool:
        """
        Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris.
        	
        """

    @make_multi_arc_ephemeris.setter
    def make_multi_arc_ephemeris(self, arg1: bool) -> None:
        ...

class InterpolatedSpiceEphemerisSettings(DirectSpiceEphemerisSettings):
    """Class for defining settings of an ephemeris interpolated from Spice data.
	
	`DirectSpiceEphemerisSettings` derived class for setting ephemerides to be created from interpolated Spice ephemeris data.
	"""

    @property
    def final_time(self) -> float:
        """
        Final time from which interpolated data from Spice should be created.
        	
        """

    @property
    def initial_time(self) -> float:
        """
        Initial time from which interpolated data from Spice should be created.
        	
        """

    @property
    def time_step(self) -> float:
        """
        Time step setting to be used for the state interpolation.
        	
        """

class KeplerEphemerisSettings(EphemerisSettings):
    """
		"""

    @property
    def central_body_gravitational_parameter(self) -> float:
        ...

    @property
    def epoch_of_initial_state(self) -> float:
        ...

    @property
    def initial_state_in_keplerian_elements(self) -> numpy.ndarray:
        ...

    @property
    def root_finder_absolute_tolerance(self) -> float:
        ...

    @property
    def root_finder_maximum_number_of_iterations(self) -> float:
        ...

class ScaledEphemerisSettings(EphemerisSettings):
    """Class for defining settings from scaling existing ephemeris settings.
	
	`EphemerisSettings` derived class for a new ephemeris created from scaling an existing ephemeris settings object. It allows the user to apply a scaling factor to the resulting Cartesian states (for instance for an uncertainty analysis).
	"""

class TabulatedEphemerisSettings(EphemerisSettings):
    """Class for defining settings of ephemeris to be created from tabulated data.
	
	`EphemerisSettings` derived class for ephemeris created from tabulated data. The provided data is interpolated into ephemerides.
	"""

    @property
    def body_state_history(self) -> dict[float, numpy.ndarray]:
        """
        Dictionary of the discrete state history data from which ephemeris is to be created.
        	
        """

def approximate_jpl_model(body_name: str) -> EphemerisSettings:
    """Factory function for creating approximate ephemeris model settings for major planets.
	
	Factory function for settings object, defining approximate ephemeris model for major planets.
	In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).
	Note that this option is only available for solar system planets. For the case of the Earth the approximate ephemeris of the Earth-Moon barycenter is returned.
	
	
	:param body_name:
			String that is attempted to be matched to an identifier for the body that the ephemeris is to be created for.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ApproximateJPLEphemerisSettings` class
	"""

def constant(constant_state: numpy.ndarray, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Factory function for creating constant ephemeris model settings.
	
	Factory function for settings object, defining ephemeris model with a constant, time-independent state.
	
	
	:param constant_state:
			Constant state that will be provided as output of the ephemeris at all times.
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.constantEphemerisSettings` class
	"""

def create_ephemeris(ephemeris_settings: EphemerisSettings, body_name: str) -> ...:
    ...

def custom(custom_state_function: typing.Callable[[float], numpy.ndarray], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    ...

def custom_ephemeris(custom_state_function: typing.Callable[[float], numpy.ndarray], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Factory function for creating custom ephemeris model settings.
	
	Factory function for settings object, defining ephemeris model with a custom state.
	This allows the user to provide an custom state function as ephemeris model.
	
	
	:param custom_state_function:
			Function returning the state as a function of time.
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.CustomEphemerisSettings` class
	"""

def direct_spice(frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', body_name_to_use: str='') -> EphemerisSettings:
    """Factory function for creating ephemeris model settings entirely from Spice.
	
	Factory function for settings object, defining ephemeris model directly and entirely from Spice.
	Requires an appropriate Spice kernel to be loaded.
	
	
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:param body_name_to_use:
			Body from which Spice ephemeris is to be created.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings` class
	"""

def interpolated_spice(initial_time: float, final_time: float, time_step: float, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=..., body_name_to_use: str='') -> EphemerisSettings:
    """Factory function for creating ephemeris model settings using interpolated Spice data.
	
	Factory function for settings object defining an ephemeris model from interpolated Spice data.
	Using this option the state of the body is retrieved from Spice at regular intervals `before` the environment propagation (as opposed to during the propagation).
	These data are then used to create an interpolator, which is put into the environment, and called during the propagation.
	This option has the downside of being applicable only during a limited time interval and requiring the tabulated data to be stored in RAM,
	but may for `some special cases <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/environment_setup/valid_time_range.html>`_
	offer an advantage over a direct Spice ephemeris (:func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.direct_spice`).
	
	
	:param initial_time:
			Initial time from which interpolated data from Spice should be created.
	:param final_time:
			Final time from which interpolated data from Spice should be created.
	:param time_step:
			Time step with which interpolated data from Spice should be created.
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:param interpolator_settings:
			Settings to be used for the state interpolation.
	:param body_name_to_use:
			Body from which Spice ephemeris is to be created.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.InterpolatedSpiceEphemerisSettings` class
	"""

def keplerian(initial_keplerian_state: numpy.ndarray, initial_state_epoch: float, central_body_gravitational_parameter: float, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', root_finder_absolute_tolerance: float=4.440892098500626e-14, root_finder_maximum_iterations: float=1000.0) -> EphemerisSettings:
    """Factory function for creating Keplerian ephemeris model settings.
	
	Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from the given Kepler elements.
	These are taken as the elements at the ``initial_state_epoch`` and propagated to any other time using the provided ``central_body_gravitational_parameter``.
	See `Frame/State Transformations <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/astrodynamics/transformations.html#frame-state-transformations>`_ and the :ref:`\\`\\`astro\\`\\`` module for more details on orbital elements in tudat.
	
	
	:param initial_state_in_keplerian_elements:
			Kepler elements at epoch given by ``initial_state_epoch``.
	
	:param initial_state_epoch:
			Epoch at which ``initial_state_epoch`` represents the Keplerian state.
	
	:param central_body_gravitational_parameter:
			Effective gravitational parameter of the central body that is used in the computations. Note that when
			the Keplerian orbit is to represent the relative state of two massive bodies, with one of these bodies as the origin
			this values should be the *sum* of the two bodies' gravitational parameters
	
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:param root_finder_absolute_tolerance:
			Convergence tolerance on iterative conversion from mean to eccentric anomaly;
			applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
	
	:param root_finder_maximum_number_of_iterations:
			Maximum iteration on iterative conversion from mean to eccentric anomaly;
			applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.KeplerEphemerisSettings` class
	"""

def keplerian_from_spice(body: str, initial_state_epoch: float, central_body_gravitational_parameter: float, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', root_finder_absolute_tolerance: float=4.440892098500626e-14, root_finder_maximum_iterations: float=1000.0) -> EphemerisSettings:
    """Factory function for creating Keplerian ephemeris model settings with initial state from Spice.
	
	Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from an initial state from Spice.
	The Kepler elements inferred from the initial state are propagated to any other time using the provided ``central_body_gravitational_parameter``.
	See `Frame/State Transformations <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/astrodynamics/transformations.html#frame-state-transformations>`_ and the :ref:`\\`\\`astro\\`\\`` module for more details on orbital elements in tudat.
	
	
	:param body:
			Name of body for which to create ephemeris settings and infer initial state from Spice.
	:param initial_state_epoch:
			Epoch at which ``initial_state_epoch`` represents the Keplerian state.
	
	:param central_body_gravitational_parameter:
			Gravitational parameter of the central body that is used in the computations.
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:param root_finder_absolute_tolerance:
			Convergence tolerance on iterative conversion from mean to eccentric anomaly;
			applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
	
	:param root_finder_maximum_number_of_iterations:
			Maximum iteration on iterative conversion from mean to eccentric anomaly;
			applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.KeplerEphemerisSettings` class
	"""

def scaled_by_constant(unscaled_ephemeris_settings: EphemerisSettings, scaling_constant: float, is_scaling_absolute: bool=False) -> EphemerisSettings:
    """Factory function for creating scaled ephemeris model settings.
	
	Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
	The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
	
	
	:param unscaled_ephemeris_settings:
			Sets base settings of ephemeris to be scaled.
	:param scaling_constant:
			Constant scaling factor to be applied to all elements of the Cartesian state.
	:param is_scaling_absolute:
			Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class
	"""

def scaled_by_vector(unscaled_ephemeris_settings: EphemerisSettings, scaling_vector: numpy.ndarray, is_scaling_absolute: bool=False) -> EphemerisSettings:
    """Factory function for creating scaled ephemeris model settings.
	
	Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
	The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
	
	
	:param unscaled_ephemeris_settings:
			Sets base settings of ephemeris to be scaled.
	:param scaling_vector:
			Vector containing scaling factors to be applied to each element of the Cartesian state.
	:param is_scaling_absolute:
			Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class
	"""

def scaled_by_vector_function(unscaled_ephemeris_settings: EphemerisSettings, scaling_vector_function: typing.Callable[[float], numpy.ndarray], is_scaling_absolute: bool=False) -> EphemerisSettings:
    """Factory function for creating scaled ephemeris model settings.
	
	Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
	The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
	
	
	:param unscaled_ephemeris_settings:
			Sets base settings of ephemeris to be scaled.
	:param scaling_vector_function:
			Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
	:param is_scaling_absolute:
			Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class
	"""

def tabulated(body_state_history: dict[float, numpy.ndarray], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Factory function for creating ephemeris model settings from tabulated data.
	
	Factory function for settings object, defining ephemeris model to be created from tabulated data.
	Currently the data that is provided gets interpolated by a 6th order Lagrange interpolator (hardcoded).
	At the edges of the interpolation interval a cubic spline interpolator is used to suppress the influence of Runge's phenomenon.
	
	
	:param body_state_history:
			Dictionary of the discrete state history data from which ephemeris is to be created. Keys representing the time (float) and values representing Cartesian states (numpy.ndarray).
	:param frame_origin:
			Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
			Orientation of frame in which ephemeris data is defined.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings` class
	"""

def tabulated_from_existing(ephemeris_settings: EphemerisSettings, start_time: float, end_time: float, time_step: float, interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=...) -> EphemerisSettings:
    """Factory function for creating tabulated ephemeris model settings from existing ephemeris.
	
	Factory function for creating tabulated ephemeris model settings from existing ephemeris.
	The ephemeris that is provided gets tabulated in a given time frame, for a given time step.
	When called, this tabulated ephemeris will use interpolation, when needed, from the specified interpolator.
	
	.. note:: Creating tabulated ephemeris from existing ephemeris can for instance be used when combined with estimation.
			  This is because estimation needs the ephemeris to be tabulated to work.
	
	
	:param ephemeris_settings:
			Existing ephemeris settings that have to be tabulated.
	:param start_time:
			Initial time for which to create the tabulated ephemeris.
	:param end_time:
			Final time for which to create the tabulated ephemeris.
	:param time_step:
			Time step to use to tabulate the existing ephemeris.
	:param interpolator_settings:
			Interpolator settings to use when interpolating between two tabulated ephemeris.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings` class
	"""