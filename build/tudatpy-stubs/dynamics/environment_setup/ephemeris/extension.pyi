import numpy
import pybind11_stubgen.typing_ext
from ....math import interpolators
import typing
__all__ = ['ApproximateJplEphemerisSettings', 'ConstantEphemerisSettings', 'CustomEphemerisSettings', 'DirectSpiceEphemerisSettings', 'EphemerisSettings', 'InterpolatedSpiceEphemerisSettings', 'KeplerEphemerisSettings', 'ScaledEphemerisSettings', 'TabulatedEphemerisSettings', 'approximate_jpl_model', 'constant', 'create_ephemeris', 'custom', 'custom_ephemeris', 'direct_spice', 'interpolated_spice', 'keplerian', 'keplerian_from_spice', 'scaled_by_constant', 'scaled_by_vector', 'scaled_by_vector_function', 'sgp4', 'tabulated', 'tabulated_from_existing']

class ApproximateJplEphemerisSettings(EphemerisSettings):
    """Class for creating settings of approximate ephemeris for major planets.
    
    `EphemerisSettings` derived class for approximate ephemeris for major planets as implemented in ApproximateJplEphemerisSettings class and derived class (described in `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_)."""

    @property
    def body_name(self) -> str:
        """
        No documentation found.
        """

class ConstantEphemerisSettings(EphemerisSettings):
    """Class for defining settings of constant ephemerides.
    
    `EphemerisSettings` derived class for ephemerides producing a constant (time-independent) state."""

class CustomEphemerisSettings(EphemerisSettings):
    """Class for defining settings of a custom ephemeris.
    
    `EphemerisSettings` derived class for ephemerides which represent an ideal Kepler orbit."""

    @property
    def get_custom_state_function(self) -> typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]]:
        """
        No documentation found.
        """

class DirectSpiceEphemerisSettings(EphemerisSettings):
    """Class for defining settings of an ephemeris linked directly to Spice.
    
    `EphemerisSettings` derived class for ephemeris which are directly linked to Spice."""

    @property
    def converge_light_time_aberration(self) -> bool:
        """
                 **read-only**
        
                 Boolean defining whether to use single iteration or max. 3 iterations for calculating light time correction.
        
                 :type: bool
        """

    @property
    def correct_for_light_time_aberration(self) -> bool:
        """
                 **read-only**
        
                 Boolean defining whether to correct for light time in retrieved values (of observed state).
        
                 :type: bool
        """

    @property
    def correct_for_stellar_aberration(self) -> bool:
        """
                 **read-only**
        
                 Boolean defining whether to correct for stellar aberrations in retrieved values (of observed state).
        
                 :type: bool
        """

class EphemerisSettings:
    """Base class for providing settings for ephemeris model.
    
    Functional (base) class for settings of ephemeris models that require no information in addition to their type (and frame origin and orientation).
    Ephemeris model classes requiring additional information must be created using an object derived from this class."""

    @property
    def ephemeris_type(self) -> ...:
        """
                 **read-only**
        
                 Type of ephemeris that is to be created.
        
                 :type: EphemerisType
        """

    @property
    def frame_orientation(self) -> str:
        """
                 Orientation of frame in which ephemeris data is to be defined.
        
                 :type: str
        """

    @frame_orientation.setter
    def frame_orientation(self, arg1: str) -> None:
        ...

    @property
    def frame_origin(self) -> str:
        """
                 Origin of frame in which ephemeris data is to be defined.
        
                 :type: str
        """

    @frame_origin.setter
    def frame_origin(self, arg1: str) -> None:
        ...

    @property
    def make_multi_arc_ephemeris(self) -> bool:
        """
                 Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris.
        
                 :type: bool
        """

    @make_multi_arc_ephemeris.setter
    def make_multi_arc_ephemeris(self, arg1: bool) -> None:
        ...

class InterpolatedSpiceEphemerisSettings(DirectSpiceEphemerisSettings):
    """Class for defining settings of an ephemeris interpolated from Spice data.
    
    `DirectSpiceEphemerisSettings` derived class for setting ephemerides to be created from interpolated Spice ephemeris data."""

    @property
    def final_time(self) -> float:
        """
                 **read-only**
        
                 Final time from which interpolated data from Spice should be created.
        
                 :type: float
        """

    @property
    def initial_time(self) -> float:
        """
                 **read-only**
        
                 Initial time from which interpolated data from Spice should be created.
        
                 :type: float
        """

    @property
    def time_step(self) -> float:
        """
                 **read-only**
        
                 Time step setting to be used for the state interpolation.
        
                 :type: float
        """

class KeplerEphemerisSettings(EphemerisSettings):
    """No documentation found."""

    @property
    def central_body_gravitational_parameter(self) -> float:
        """
        No documentation found.
        """

    @property
    def epoch_of_initial_state(self) -> float:
        """
        No documentation found.
        """

    @property
    def initial_state_in_keplerian_elements(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
        No documentation found.
        """

    @property
    def root_finder_absolute_tolerance(self) -> float:
        """
        No documentation found.
        """

    @property
    def root_finder_maximum_number_of_iterations(self) -> float:
        """
        No documentation found.
        """

class ScaledEphemerisSettings(EphemerisSettings):
    """Class for defining settings from scaling existing ephemeris settings.
    
    `EphemerisSettings` derived class for a new ephemeris created from scaling an existing ephemeris settings object. It allows the user to apply a scaling factor to the resulting Cartesian states (for instance for an uncertainty analysis)."""

class TabulatedEphemerisSettings(EphemerisSettings):
    """Class for defining settings of ephemeris to be created from tabulated data.
    
    `EphemerisSettings` derived class for ephemeris created from tabulated data. The provided data is interpolated into ephemerides."""

    @property
    def body_state_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]]:
        """
                 **read-only**
        
                 Dictionary of the discrete state history data from which ephemeris is to be created.
        
                 :type: Dict[[float], numpy.ndarray[numpy.float64[6, 1]]]
        """

def approximate_jpl_model(body_name: str) -> EphemerisSettings:
    """Function for creating approximate ephemeris model settings for major planets.
    
    Function for settings object, defining approximate ephemeris model for major planets.
    In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described in `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).
    Note that this option is only available for solar system planets. For the case of the Earth the approximate ephemeris of the Earth-Moon barycenter is returned.
    
    
    Parameters
    ----------
    body_name : str
        String that is attempted to be matched to an identifier for the body that the ephemeris is to be created for.
    Returns
    -------
    ApproximateJplEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.ApproximateJplEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for Jupiter using JPL's approximate planet position model:
    
    .. code-block:: python
    
       # create ephemeris settings and add to body settings of body "Jupiter"
       body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.approximate_jpl_model( "Jupiter" )
    
    
    Alternatively, we can assign the ApproximateJplEphemerisSettings of Jupiter (or any other body for which an approximate JPL ephemeris is available) to any custom body:
    
    .. code-block:: python
    
       # create ephemeris settings (normally used for Jupiter) and add to body settings of body "CustomBody"
       body_settings.get( "CustomBody" ).ephemeris_settings = environment_setup.ephemeris.approximate_jpl_model( "Jupiter" )
    
                                   ephemerides::ApproximatePlanetPositionsBase::jupiter, false );"""

def constant(constant_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Function for creating constant ephemeris model settings.
    
    Function for settings object, defining ephemeris model with a constant, time-independent state.
    
    
    Parameters
    ----------
    constant_state : numpy.ndarray[numpy.float64[6, 1]]
        Constant state that will be provided as output of the ephemeris at all times.
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    Returns
    -------
    ConstantEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.ConstantEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for a time-independent, constant state of Jupiter:
    
    .. code-block:: python
    
       # Define the constant cartesian state
       constant_cartesian_state = [100.0e9, 100.0e9, 100.0e9, 10.0e3, 10.0e3, 10.0e3]
       # Define the ephemeris frame
       frame_origin = "SSB"
       frame_orientation = "J2000"
       # Make the ephemeris settings
       body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.constant(
         constant_cartesian_state,
         frame_origin, frame_orientation)"""

def create_ephemeris(ephemeris_settings: EphemerisSettings, body_name: str) -> ...:
    """No documentation found."""

def custom(custom_state_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    ...

def custom_ephemeris(custom_state_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Function for creating custom ephemeris model settings.
    
    Function for settings object, defining ephemeris model with a custom state.
    This allows the user to provide a custom state function as ephemeris model.
    
    
    Parameters
    ----------
    custom_state_function : Callable[[float], numpy.ndarray[numpy.float64[6, 1]]]
        Function returning the state as a function of time.
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    Returns
    -------
    CustomEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.CustomEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for Earth from a custom state history function:
    
    .. code-block:: python
    
      # Define the custom state function for Earth
      def custom_state_function(time):
          # Compute what fraction of the year it is
          frac_year = (time - 2451545) % (365.25*24*3600)
          # Distance and velocity of the Earth w.r.t. the Sun
          AU, v_E = 1.496e11, 30e3
          # Compute the position and velocity of the Earth in a 2D circle
          x_pos = np.sin(frac_year*np.pi) * AU
          y_pos = np.cos(frac_year*np.pi) * AU
          x_vel = np.cos(frac_year*np.pi) * v_E
          y_vel = np.sin(frac_year*np.pi) * v_E
          return [x_pos, y_pos, 0, x_vel, y_vel, 0]
      # Define the ephemeris frame
      frame_origin = "SSB"
      frame_orientation = "J2000"
      # Make the ephemeris settings
      body_settings.get("Earth").ephemeris_settings = environment_setup.ephemeris.custom(
          custom_state_function,
          frame_origin,
          frame_orientation)"""

def direct_spice(frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', body_name_to_use: str='') -> EphemerisSettings:
    """Function for creating ephemeris model settings entirely from Spice.
    
    Function for settings object, defining ephemeris model directly and entirely from Spice.
    Requires an appropriate Spice kernel to be loaded.
    
    
    Parameters
    ----------
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    body_name_to_use : str, default = ""
        Body from which Spice ephemeris is to be created (if empty, it uses the name of the body to which the settings are assigned, see below for example).
    Returns
    -------
    DirectSpiceEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.DirectSpiceEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create barycentric (origin: SSB) :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` with axes along J2000, using data directly from spice:
    
    .. code-block:: python
    
       frame_origin = "SSB"
       frame_orientation = "J2000"
       body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.direct_spice(
           frame_origin, frame_orientation)
    
    
    Alternatively, we can assign the DirectSpiceEphemerisSettings of Jupiter (or any other body for which a direct Spice ephemeris is available) to any custom body:
    
    .. code-block:: python
    
       frame_origin = "SSB"
       frame_orientation = "J2000"
       body_name_to_use =  "Jupiter"
       # create ephemeris settings from "Jupiter" spice data and add to body settings of body "CustomBody"
       body_settings.get( "CustomBody" ).ephemeris_settings = environment_setup.ephemeris.direct_spice(
           frame_origin, frame_orientation, body_name_to_use )"""

def interpolated_spice(initial_time: float, final_time: float, time_step: float, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', interpolator_settings: interpolators.InterpolatorSettings=..., body_name_to_use: str='') -> EphemerisSettings:
    """Function for creating ephemeris model settings using interpolated Spice data.
    
    Function for settings object defining an ephemeris model from interpolated Spice data.
    Using this option the state of the body is retrieved from Spice at regular intervals `before` the environment propagation (as opposed to during the propagation).
    These data are then used to create an interpolator, which is put into the environment, and called during the propagation.
    This option has the downside of being applicable only during a limited time interval and requiring the tabulated data to be stored in RAM,
    but may for `some special cases <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models/default_bodies_limited_time_range.html>`_
    offer an advantage over a direct Spice ephemeris (:func:`~tudatpy.dynamics.environment_setup.ephemeris.direct_spice`).
    
    
    Parameters
    ----------
    initial_time : float
        Initial time from which interpolated data from Spice should be created.
    final_time : float
        Final time from which interpolated data from Spice should be created.
    time_step : float
        Time step with which interpolated data from Spice should be created.
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    interpolator_settings : std::make_shared< interpolators::InterpolatorSettings >, default=std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 )
        Settings to be used for the state interpolation.
    body_name_to_use : str, default = ""
        Body from which Spice ephemeris is to be created.
    Returns
    -------
    InterpolatedSpiceEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.DirectSpiceEphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.InterpolatedSpiceEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we define :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for Jupiter by retrieving ephemeris data from Spice at 3600 s intervals between t=0 and t=1.0E8:
    
    .. code-block:: python
    
      # Define the interpolation settings
      initial_time = 0.0
      final_time = 1.0E8
      time_step = 3600.0
      # Define the ephemeris frame
      frame_origin = "SSB"
      frame_orientation = "J2000"
      # create ephemeris settings and add to body settings of body "Jupiter"
      body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.interpolated_spice(
        initial_time, final_time, time_step, frame_origin, frame_orientation )
    
    
    By default, a 6th order Lagrange interpolator is used (NOTE: the Lagrange interpolator is not reliable at the edges of the interpolation interval, as discussed here: :func:`~tudatpy.math.interpolators.lagrange_interpolation`).
    Settings for an alternative interpolator can be use by specifying the optional input argument.
    Additionally, as is the case for the :func:`~tudatpy.dynamics.environment_setup.ephemeris.direct_spice` and :func:`~tudatpy.dynamics.environment_setup.ephemeris.approximate_jpl_model` functions, an optional input argument ``body_name_to_use`` allows to use an ephemeris model from Spice for some body and assign it to a custom body."""

def keplerian(initial_keplerian_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], initial_state_epoch: float, central_body_gravitational_parameter: float, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', root_finder_absolute_tolerance: float=4.440892098500626e-14, root_finder_maximum_iterations: float=1000.0) -> EphemerisSettings:
    """Function for creating Keplerian ephemeris model settings.
    
    Function for settings object, defining ephemeris model which represents an ideal Kepler orbit from the given Kepler elements.
    These are taken as the elements at the ``initial_state_epoch`` and propagated to any other time using the provided ``central_body_gravitational_parameter``.
    See `Element Types <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#element-types>`_ and the :ref:`astro` module for more details on orbital elements in tudat.
    
    
    Parameters
    ----------
    initial_state_in_keplerian_elements : numpy.ndarray[numpy.float64[6, 1]]
        Kepler elements at epoch given by ``initial_state_epoch``.
    
    initial_state_epoch : float
        Epoch at which ``initial_state_epoch`` represents the Keplerian state.
    
    central_body_gravitational_parameter : float
        Effective gravitational parameter of the central body that is used in the computations. Note that when
        the Keplerian orbit is to represent the relative state of two massive bodies, with one of these bodies as the origin
        this values should be the *sum* of the two bodies' gravitational parameters
    
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    root_finder_absolute_tolerance : float
        Convergence tolerance on iterative conversion from mean to eccentric anomaly;
        applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
    
    root_finder_maximum_number_of_iterations : float
        Maximum iteration on iterative conversion from mean to eccentric anomaly;
        applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
    
    Returns
    -------
    KeplerEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.KeplerEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for a simple, barycentric (SSB) Kepler orbit of Jupiter:
    
    .. code-block:: python
    
      # Define the computation of the Kepler orbit ephemeris
      initial_state_in_keplerian_elements = [100e9, 0.7, 1.0, 2.0, 2.0, 2.0]
      initial_state_epoch = 12345
      central_body_gravitational_parameter = 1.3284e20 # (sum of Sun and Jupiter)
      # Define the ephemeris frame
      frame_origin = "SSB"
      frame_orientation = "J2000"
      # Create ephemeris settings and add to body settings of "Jupiter"
      body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.keplerian(
          initial_state_in_keplerian_elements,
          initial_state_epoch,
          central_body_gravitational_parameter,
          frame_origin, frame_orientation )"""

def keplerian_from_spice(body: str, initial_state_epoch: float, central_body_gravitational_parameter: float, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000', root_finder_absolute_tolerance: float=4.440892098500626e-14, root_finder_maximum_iterations: float=1000.0) -> EphemerisSettings:
    """Function for creating Keplerian ephemeris model settings with initial state from Spice.
    
    Function for settings object, defining ephemeris model which represents an ideal Kepler orbit from an initial state from Spice.
    The Kepler elements inferred from the initial state are propagated to any other time using the provided ``central_body_gravitational_parameter``.
    See `Element Types <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#element-types>`_ and the :ref:`astro` module for more details on orbital elements in tudat.
    
    
    Parameters
    ----------
    body : str
        Name of body for which to create ephemeris settings and infer initial state from Spice.
    initial_state_epoch : float
        Epoch at which ``initial_state_epoch`` represents the Keplerian state.
    
    central_body_gravitational_parameter : float
        Gravitational parameter of the central body that is used in the computations.
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    root_finder_absolute_tolerance : float
        Convergence tolerance on iterative conversion from mean to eccentric anomaly;
        applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
    
    root_finder_maximum_number_of_iterations : float
        Maximum iteration on iterative conversion from mean to eccentric anomaly;
        applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.
    
    Returns
    -------
    KeplerEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.KeplerEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for a simple, barycentric (SSB) Kepler orbit of Jupiter.
    The initial keplerian state is extracted from Spice as the state of ``body_name`` w.r.t. ``frame_origin``
    
    .. code-block:: python
    
      # Define the parameters for retrieval of the initial Kepler orbit elements from spice
      body_name = "Jupiter"
      initial_state_epoch = 12345
      central_body_gravitational_parameter = 1.3284e20 # (sum of Sun and Jupiter)
      # Define the ephemeris frame
      frame_origin = "SSB"
      frame_orientation = 'J2000'
      # Make ephemeris the settings and add to body settings of "Jupiter"
      body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.keplerian_from_spice(
          body_name,
          initial_state_epoch,
          central_body_gravitational_parameter,
          frame_origin,
          frame_orientation )
    
    
    Additionally, as is the case for the :func:`~tudatpy.dynamics.environment_setup.ephemeris.direct_spice`, :func:`~tudatpy.dynamics.environment_setup.ephemeris.approximate_jpl_model` and :func:`~tudatpy.dynamics.environment_setup.ephemeris.interpolated_spice` functions, the ephemeris model from Spice can be retrieved for some body and assigned to a custom body."""

def scaled_by_constant(unscaled_ephemeris_settings: EphemerisSettings, scaling_constant: float, is_scaling_absolute: bool=False) -> EphemerisSettings:
    """Function for creating scaled ephemeris model settings.
    
    Function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
    The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
    
    
    Parameters
    ----------
    unscaled_ephemeris_settings : EphemerisSettings
        Sets base settings of ephemeris to be scaled.
    scaling_constant : float
        Constant scaling factor to be applied to all elements of the Cartesian state.
    is_scaling_absolute : bool, default=false
        Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
    Returns
    -------
    ScaledEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.ScaledEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettingsObject` with a constant factor:
    
    .. code-block:: python
    
       # define variable for scaling factor
       scaling_constant = 1.001
       # define variables containing the existing ephemeris settings
       unscaled_ephemeris_settings = body_settings.get( "Jupiter" ).ephemeris_settings
       # make new ephemeris settings
       body_settings.get( "Jupiter" ).ephemeris_settings =  environment_setup.ephemeris.scaled_by_constant(
              unscaled_ephemeris_settings, scaling_constant )
    
    In the above case, the original Jupiter ephemeris setting is taken and each state element (x,y,z position and velocity) from the original ephemeris is multiplied by a factor 1.001."""

def scaled_by_vector(unscaled_ephemeris_settings: EphemerisSettings, scaling_vector: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], is_scaling_absolute: bool=False) -> EphemerisSettings:
    """Function for creating scaled ephemeris model settings.
    
    Function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
    The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
    
    
    Parameters
    ----------
    unscaled_ephemeris_settings : EphemerisSettings
        Sets base settings of ephemeris to be scaled.
    scaling_vector : numpy.ndarray[numpy.float64[6, 1]]
        Vector containing scaling factors to be applied to each element of the Cartesian state.
    is_scaling_absolute : bool, default=false
        Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
    Returns
    -------
    ScaledEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.ScaledEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettingsObject` with the constant elements of a vector:
    
    .. code-block:: python
    
      # Define the scaling vector
      scaling_vector = [1.01, 0.99, 1, 1, 1, 0]
      # Extract the unscaled ephemeris settings from Jupiter
      unscaled_ephemeris_settings = body_settings.get( "Jupiter" ).ephemeris_settings
      # Create the scaled ephemeris settings and apply to the body "Jupiter"
      body_settings.get( "Jupiter" ).ephemeris_settings =  environment_setup.ephemeris.scaled_by_vector(
          unscaled_ephemeris_settings,
          scaling_vector)
    
    In the above case, the original Jupiter ephemeris setting is taken and each state element (x,y,z position and velocity) from the original ephemeris is multiplied by the corresponding scaling factor in ``scaling_vector``."""

def scaled_by_vector_function(unscaled_ephemeris_settings: EphemerisSettings, scaling_vector_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]], is_scaling_absolute: bool=False) -> EphemerisSettings:
    """Function for creating scaled ephemeris model settings.
    
    Function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
    The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
    
    
    Parameters
    ----------
    unscaled_ephemeris_settings : EphemerisSettings
        Sets base settings of ephemeris to be scaled.
    scaling_vector_function : callable[[float], numpy.ndarray[numpy.float64[6, 1]]]
        Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
    is_scaling_absolute : bool, default=false
        Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
    Returns
    -------
    ScaledEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.ScaledEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` object with factors from a custom function:
    
    .. code-block:: python
    
      # Define the scaling vector function
      def scaling_vector_function(time):
          # Add a wobble in the x and y coordinates
          wobble = 1 + 0.1 * np.cos(time/50)
          return [wobble, wobble, 1, 1, 1, 1]
      # Extract the existing unscaled ephemeris settings
      unscaled_ephemeris_settings = body_settings.get( "Jupiter" ).ephemeris_settings
      # Create the scaled ephemeris settings and apply to the body "Jupiter"
      body_settings.get( "Jupiter" ).ephemeris_settings =  environment_setup.ephemeris.scaled_by_vector_function(
          unscaled_ephemeris_settings,
          scaling_vector_function )
    
    In the above case, the original Jupiter ephemeris setting is taken and each state element (x,y,z position and velocity) from the original ephemeris is multiplied by the corresponding scaling factor in the vector returned by ``vector_scaling_function``."""

def sgp4(tle_line_1: str, tle_line_2: str, frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Function for creating ephemeris model settings for an SGP4-propagated TLE.
    
    Function for creating ephemeris model settings for an SGP4-propagated two-line element (TLE). Our implementation uses the `evsgp4_c`_ function of the SPICE library
    to perform the SGP4 propagation, and the :func:`~tudatpy.astro.element_conversion.teme_to_j2000` function to rotate the resulting state from the TEME frame to the J2000 frame
    (and, if required for this ephemeris model, a subsequent different inertial frame).
    
    .. _`evsgp4_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/evsgp4_c.html
    
    Parameters
    ----------
    tle_line_1 : str
        First line of the two-line element set
    tle_line_2 : str
        second line of the two-line element set
    frame_origin : str, default="Earth"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="J2000"
        Orientation of frame in which ephemeris data is defined.
    Returns
    -------
    DirectTleEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.DirectTleEphemerisSettings` class
    
    
    
    Examples
    --------
    In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` object with factors from a custom function:
    
    .. code-block:: python
    
      # Create ephemeris settings for Delfi C-3 spacecraft using launch TLE
       body_settings.get( "DelfiC-3" ).ephemeris_settings =  environment_setup.ephemeris.sgp4(
            '1 32789U 07021G   08119.60740078 -.00000054  00000-0  00000+0 0  9999',
            '2 32789 098.0082 179.6267 0015321 307.2977 051.0656 14.81417433    68' )
     )
    
    In the above case, the given TLE will be used as input to an SGP4 propagation, with the output given from this ephemeris model in the (default) Earth-centered J2000 frame."""

def tabulated(body_state_history: dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]], frame_origin: str='SSB', frame_orientation: str='ECLIPJ2000') -> EphemerisSettings:
    """Function for creating ephemeris model settings from tabulated data.
    
    Function for settings object, defining ephemeris model to be created from tabulated data.
    Currently the data that is provided gets interpolated by a 6th order Lagrange interpolator (hardcoded).
    At the edges of the interpolation interval a cubic spline interpolator is used to suppress the influence of Runge's phenomenon.
    
    
    Parameters
    ----------
    body_state_history : dict
        Dictionary of the discrete state history data from which ephemeris is to be created. Keys representing the time (float) and values representing Cartesian states (numpy.ndarray).
    frame_origin : str, default="SSB"
        Origin of frame in which ephemeris data is defined.
    frame_orientation : str, default="ECLIPJ2000"
        Orientation of frame in which ephemeris data is defined.
    Returns
    -------
    TabulatedEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.TabulatedEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for Jupiter from tabulated state history data:
    
    .. code-block:: python
    
      # Define the Dict containing Jupiter's tabulated state history
      body_state_history = {
          0: [7.4713e11, 0, 0, 13.5e3, 0, 0],
          1000: [7.4711e11, 5e9, 0, 13.4998e3, 75, 0],
          2150: [7.4671e11, 2.5e10, 0, 13.498e3, 200, 0],
          # ... truncated
          15650: [7.3899e11, 1.1e11, 0, 13.416e3, 1.5e3, 0]
      }
      # Define the ephemeris frame
      frame_origin = "SSB"
      frame_orientation = "J2000"
      # Create the tabulated ephemeris settings and add them to the body "Jupiter"
      body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.tabulated( body_state_history,
          frame_origin,
          frame_orientation )"""

def tabulated_from_existing(ephemeris_settings: EphemerisSettings, start_time: float, end_time: float, time_step: float, interpolator_settings: interpolators.InterpolatorSettings=...) -> EphemerisSettings:
    """Function for creating tabulated ephemeris model settings from existing ephemeris.
    
    Function for creating tabulated ephemeris model settings from existing ephemeris.
    The ephemeris that is provided gets tabulated in a given time frame, for a given time step.
    When called, this tabulated ephemeris will use interpolation, when needed, from the specified interpolator.
    
    .. note:: Creating tabulated ephemeris from existing ephemeris can for instance be used when combined with estimation.
              This is because estimation needs the ephemeris to be tabulated to work.
    
    
    Parameters
    ----------
    ephemeris_settings : tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings
        Existing ephemeris settings that have to be tabulated.
    start_time : float
        Initial time for which to create the tabulated ephemeris.
    end_time : float
        Final time for which to create the tabulated ephemeris.
    time_step : float
        Time step to use to tabulate the existing ephemeris.
    interpolator_settings : tudatpy.math.interpolators.InterpolatorSettings, default=tudatpy.math.interpolators.lagrange_interpolation(8)
        Interpolator settings to use when interpolating between two tabulated ephemeris.
    Returns
    -------
    TabulatedEphemerisSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.dynamics.environment_setup.ephemeris.TabulatedEphemerisSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.ephemeris.EphemerisSettings` for Io.
    First, we extract the existing ephemeris. Then, we define new tabulated ephemeris settings, from the original settings.
    
    .. code-block:: python
    
      # Get the original ephemeris settings
      original_io_ephemeris_settings = body_settings.get( "Io" ).ephemeris_settings
      # Apply new tabulated ephemeris settings
      body_settings.get( "Io" ).ephemeris_settings =  environment_setup.ephemeris.tabulated_from_existing(
        original_io_ephemeris_settings,
        initial_time,
        final_time,
        time_step )"""