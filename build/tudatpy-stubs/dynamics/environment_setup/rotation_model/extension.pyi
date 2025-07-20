import numpy
import pybind11_stubgen.typing_ext
from ....math import interpolators
import typing
__all__ = ['IAUConventions', 'IAURotationModelSettings', 'PlanetaryRotationModelSettings', 'RotationModelSettings', 'RotationModelType', 'SimpleRotationModelSettings', 'aerodynamic_angle_based', 'constant_rotation_model', 'custom_inertial_direction_based', 'custom_rotation_model', 'gcrs_to_itrs', 'gcrs_to_itrs_rotation_model', 'iau_2000_a', 'iau_2000_b', 'iau_2006', 'iau_rotation_model', 'mars_high_accuracy', 'orbital_state_direction_based', 'planetary_rotation_model', 'simple', 'simple_from_spice', 'simple_rotational_model', 'spice', 'spice_rotation_model', 'synchronous', 'synchronous_rotation_model', 'zero_pitch_moment_aerodynamic_angle_based']

class IAUConventions:
    """Enumeration of IAU conventions for Earth rotation.
    
             Enumeration of IAU conventions for Earth rotation supported by tudat.
    
    
    
    
    
          
    
    Members:
    
      iau_2000_a : 
          
    
      iau_2000_b : 
          
    
      iau_2006 : 
          """
    __members__: typing.ClassVar[dict[str, IAUConventions]]
    iau_2000_a: typing.ClassVar[IAUConventions]
    iau_2000_b: typing.ClassVar[IAUConventions]
    iau_2006: typing.ClassVar[IAUConventions]

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

class IAURotationModelSettings(RotationModelSettings):
    """No documentation found."""

class PlanetaryRotationModelSettings(RotationModelSettings):
    """No documentation found."""

class RotationModelSettings:
    """Base class for providing settings for automatic rotation model creation.
    
    This class is a functional base class for settings of rotation models that require no information in addition to their type.
    Basic rotation model has constant orientation of the rotation axis (body-fixed z-axis) and constant rotation rate about this axis.
    Rotation models requiring additional information must be created using the functions which create the specific object derived from this base class."""

    @property
    def base_frame(self) -> str:
        """
                 Name of the base frame of rotation model.
        
                 :type: str
        """

    @base_frame.setter
    def base_frame(self, arg1: str) -> None:
        ...

    @property
    def rotation_type(self) -> RotationModelType:
        """
                 **read-only**
        
                 Type of rotation model that is to be created.
        
                 :type: RotationModelType
        """

    @property
    def target_frame(self) -> str:
        """
                 **read-only**
        
                 Name of the target frame of rotation model.
        
                 :type: str
        """

class RotationModelType:
    """Enumeration of rotation model types.
    
             Enumeration of rotation model types supported by tudat.
    
    
    
    
    
          
    
    Members:
    
      simple_rotational_model : No documentation found.
    
      spice_rotation_model : 
          
    
      gcrs_to_itrs_rotation_model : 
          
    
      synchronous_rotation_model : 
          
    
      planetary_rotation_model : 
          """
    __members__: typing.ClassVar[dict[str, RotationModelType]]
    gcrs_to_itrs_rotation_model: typing.ClassVar[RotationModelType]
    planetary_rotation_model: typing.ClassVar[RotationModelType]
    simple_rotational_model: typing.ClassVar[RotationModelType]
    spice_rotation_model: typing.ClassVar[RotationModelType]
    synchronous_rotation_model: typing.ClassVar[RotationModelType]

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
    """No documentation found."""

def aerodynamic_angle_based(central_body: str, base_frame: str, target_frame: str, angle_funcion: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=None) -> RotationModelSettings:
    """Function for creating rotation model settings based on custom aerodynamic angles (attack, sideslip, bank).
    
    Function for creating rotation model settings based on custom aerodynamic angles:
    angle of attack :math:`\\alpha`, sideslip angle :math:`\\beta` and bank angle :math:`\\sigma`. The use of this function is typical for
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
    
    
    Parameters
    ----------
    central_body : str
        Name of the central body C that is to be used.
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    angle_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]], default = None
        Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\\alpha,\\beta,\\sigma]`. If this input is left empty, these angles are both fixed to 0.
    Returns
    -------
    CustomRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model."""

def constant_rotation_model(base_frame: str, target_frame: str, initial_orientation: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]) -> RotationModelSettings:
    """Function for creating simple rotation model settings for target-frames with constant orientation.
    
    Function for settings object, defining simple rotation model setting objects with constant rotation matrix.
    These model settings are for target frames which do not have a rotational rate in the base frame and are fully defined by their initial orientation.
    
    
    Parameters
    ----------
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    initial_orientation : numpy.ndarray[numpy.float64[3, 3]]
        Rotation matrix from inertial to body-fixed (base to target) frame at initial time (constant throughout).
    Returns
    -------
    SimpleRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.SimpleRotationModelSettings` class.
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` for Earth,
    using a constant rotation matrix between Earth-fixed and inertial frame:
    
    .. code-block:: python
    
      # define parameters describing the constant orientation between frames
      original_frame = "ECLIPJ2000"
      target_frame = "Earth_fixed"
      constant_orientation = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
      # create rotation model settings and assign to body settings of "Earth"
      body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.constant(
          original_frame,
          target_frame,
          constant_orientation )"""

def custom_inertial_direction_based(inertial_body_axis_direction: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]], base_frame: str, target_frame: str, free_rotation_angle_function: typing.Callable[[float], float]=None) -> RotationModelSettings:
    """Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction
    
    Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction :math:`\\hat{\\mathbf{T}}_{I}`. Specifically, it ensures
    that the rotation matrix from body-fixed to inertial frame is set up such that :math:`\\hat{\\mathbf{T}}_{I}=R^{(I/B)}\\hat{\\mathbf{i}}` (where :math:`\\mathbf{i}` is the unit-vector in local x-direction).
    The complete rotation matrix requires an additional angle :math:`\\phi` (rotation of the body about its body-fixed x-axis), which is set to 0 by default.
    
    The full rotation matrix is computed from a 3-2-1 Euler angle rotation
    :math:`R^{(I/B)}=R_{z}(\\psi)R_{y}(\\theta)R_{x}(\\phi)`, with :math:`\\psi` and :math:`\\theta` computed from the suitable decomposition of :math:`\\hat{\\mathbf{T}}_{I}`.
    This function is typically used for simulating the (guided) dynamics of a spacecraft under thrust, where the thrust is provided in the x-direction of the body-fixed frame. By providing a suitable
    ``inertial_body_axis_direction``, this thrust can be defined to point in an arbitrary direction (typically defined by a guidance algorithm) in the inertial frame as a function of time.
    
    NOTE: this function may be extended in the future to allow an arbitrary body-fixed direction to align with an arbitrary inertial direction. At present, its functionality is limited to imposing the inertial direction of the body-fixed x-axis.
    
    
    Parameters
    ----------
    inertial_body_axis_direction : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
        Custom function defined by the user, which imposes the inertial orientation of the body-fixed x-axis, by providing :math:`\\hat{\\mathbf{T}}_{I}(t)`.
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    free_rotation_angle_function : Callable[[float], float], default = None
        Custom function provided by the user, which returns a value for the free rotation angle :math:`\\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
    Returns
    -------
    BodyFixedDirectionBasedRotationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model."""

def custom_rotation_model(base_frame: str, target_frame: str, custom_rotation_matrix_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]], finite_difference_time_step: float) -> RotationModelSettings:
    """Function for creating rotation model settings based on custom definition of rotation matrix.
    
    Function for creating rotation model settings based on custom definition of rotation matrix. The user provides a custom function that computes the rotation matrix
    from body-fixed to inertial frame as a function of time. This function can
    depend on any quantities of the user's choosing, for details on how to link the properties of the environment to this function, see `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/custom_models.html>`_.
    Since this function only computes the rotation matrix directly, the rotation matrix time derivative (and consequently, the angular velocity) are computed numerically, using a second
    order finite-difference method. Note that this computation of time derivative will only take into account the explicit time-dependence of thh custom rotation matrix.
    
    Parameters
    ----------
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    custom_rotation_matrix_function: Callable[[float], numpy.ndarray[numpy.float64[3, 3]]]
        Function computing the body-fixed to inertial rotation matrix as a function of time
    finite_difference_time_step: float
        Step size to use when computing the rotation matrix derivative numerically
    Returns
    -------
    :class:`~tudatpy.dynamics.environment_setup.rotation_model.CustomRotationModelSettings`
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model."""

def gcrs_to_itrs(precession_nutation_theory: IAUConventions=..., base_frame: str='GCRS', cio_interpolation_settings: interpolators.InterpolatorGenerationSettings=None, tdb_to_tt_interpolation_settings: interpolators.InterpolatorGenerationSettings=None, short_term_eop_interpolation_settings: interpolators.InterpolatorGenerationSettings=None) -> RotationModelSettings:
    """Function for creating high-accuracy Earth rotation model settings.
    
    Function for settings object, defining high-accuracy Earth rotation model according to the IERS Conventions 2010.  The model
    computes the rotation from ITRS to GCRS (with rotation matrix :math:`\\mathbf{R}^{(\\text{GCRS}/\\text{ITRS})}`) and its inverse from:
    
    .. math::
       \\mathbf{R}^{(\\text{GCRS}/\\text{ITRS})} = \\mathbf{R}^{(\\text{GCRS}/\\text{CIRS})}(X,Y,s)\\mathbf{R}^{(\\text{CIRS}/\\text{TIRS})}(\\theta_{E})\\mathbf{R}^{(\\text{TIRS}/\\text{ITRS})}(x_{p}, y_{p}, s')
    
    using the intermediate frames TIRS (Terrestial Intermediate Reference System) and CIRS (Celestial Intermediate Reference System) where (with equations referring to the IERS 2010 Conventions) :math:`\\mathbf{R}^{(\\text{GCRS}/\\text{CIRS})}` implements Eq. (5.10), :math:`\\mathbf{R}^{(\\text{CIRS}/\\text{TIRS})}` implements Eq. (5.5), and
    :math:`\\mathbf{R}^{(\\text{TIRS}/\\text{ITRS})}` implements Eq. (5.3). The inputs to these rotation matrices are :
    
    * :math:`X`, :math:`Y`: Celestial pole position elements
    * :math:`s`: The CIO (Celestial Intermediate Origin)
    * :math:`\\theta_{E}`: Earth rotation angle (denoted as :math:`ERA` in IERS Conventions)
    * :math:`x_{p}`, :math:`y_{p}`: Polar motion components
    * :math:`s'`: The TIO (Terrestial Intermediate Origin)
    
    Depending on the selected ``precession_nutation_theory`` input, the SOFA function ``iauXys00a``, ``iauXys00b`` or ``iauXys06a`` is used to compute :math:`X,Y,s`, when selecting
    :class:`~tudatpy.dynamics.environment_setup.rotation_model.IAUConventions` ``iau_2000a``, ``iau_2000b`` or ``iau_2006``, respectively. For epoch 01-01-1962 and later, corrections to the nominal values of :math:`X,Y`
    are applied using linear interpolation of daily corrections for :math:`X,Y` from the eopc04_14_IAU2000.62-now.txt file. The quantity :math:`s'` is computed from Eq. (5.13) (implemented in SOFA's ``iauSp00`` function).
    
    The value of :math:`\\theta_{E}` is computed directly from UTC-UT1, which is computed using settings given in :func:`~tudatpy.astro.time_representation.default_time_scale_converter`, the computation of
    :math:`\\theta_{E}` from this quantity follows from Eq. (5.15), implemented by SOFA's ``iauEra00`` function.
    
    The polar motion components :math:`x_{p}`, :math:`y_{p}` are computed from:
    
    * Corrections for semi-diurnal variations due to libration for a non-rigid Earth as per Table 5.1a (with :math:`n=2`) of IERS Conventions 2010
    * Corrections diurnal and semidiurnal variations due to ocean tides as per Tables 8.1a and 8.1b of the IERS Conventions 2010
    * Linear interpolation (correcting for discontinuities during days with leap seconds) of daily corrections for :math:`x_{p}, y_{p}`: from the eopc04_14_IAU2000.62-now.txt file in the tudat-resources directory (for epoch 01-01-1962 and later, zero otherwise)
    
    Note that for this model the original frame must be J2000 or GCRS (in the case of the former, the frame bias between GCRS and J2000 is automatically corrected for). The target frame (e.g. body-fixed frame) name is ITRS.
    The target frame (e.g. body-fixed frame) name is ITRS.
    
    Alternative options to modify the input (not exposed here) include the EOP correction file, input time scale, short period UT1 and polar motion variations.
    
    Parameters
    ----------
    precession_nutation_theory : IAUConventions, default=tba::iau_2006
        Setting theory for modelling Earth nutation.
    base_frame : str, default='GCRS'
        Base frame of rotation model
    Returns
    -------
    GcrsToItrsRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.GcrsToItrsRotationModelSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` for Earth,
    using a high-accuracy Earth rotation model as defined by IERS Conventions 2010:
    
    
    .. code-block:: python
    
       # define parameters describing the rotation between frames
       precession_nutation_theory = environment_setup.rotation_model.IAUConventions.iau_2006
       original_frame = "J2000"
       # create rotation model settings and assign to body settings of "Earth"
       body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs(
       precession_nutation_theory, original_frame)"""

def iau_rotation_model(base_frame: str, target_frame: str, nominal_meridian: float, nominal_pole: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(2, 1)], rotation_rate: float, pole_precession: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(2, 1)], merdian_periodic_terms: dict[float, tuple[float, float]], pole_periodic_terms: dict[float, tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(2, 1)], float]]) -> IAURotationModelSettings:
    ...

def mars_high_accuracy(base_frame: str='ECLIPJ2000', target_frame: str='Mars_Fixed') -> RotationModelSettings:
    """Function for creating a high-accuracy Mars rotation model.
    
    Function for creating a high-accuracy Mars rotation model, using the default parameters of `Konopliv et al. (2016) <https://www.sciencedirect.com/science/article/abs/pii/S0019103516001305>`_
    and the mathematical model of `Konopliv et al. (2006) <https://www.sciencedirect.com/science/article/pii/S0019103506000297>`_ . The rotation matrix formulation is given in Eq. (13)-(19) of that paper.
    Note that, at the moment, all parameters in this rotation model are hard coded, and cannot be adapted by the user (except by estimating a number of its constituent parameters, see :ref:`parameter` module )
    As such, this model is at present applicable to Mars rotation only. If you require more fine-grained control of the parameters, please contact the Tudat support team
    
    
    Parameters
    ----------
    base_frame : str, default = "ECLIPJ2000"
        Name of the base frame of rotation model.
    target_frame : str, default = "Mars_Fixed"
        Name of the target frame of rotation model.
    Returns
    -------
    RotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.PlanetaryRotationModelSettings` class, which defines the required settings for the rotation model."""

def orbital_state_direction_based(central_body: str, is_colinear_with_velocity: bool, direction_is_opposite_to_vector: bool, base_frame: str, target_frame: str='', free_rotation_angle_function: typing.Callable[[float], float]=None) -> RotationModelSettings:
    """Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector.
    
    Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector. This function is
    similar to the :func:`~custom_inertial_direction_based` function, with the exception that the :math:`\\hat{\\mathbf{T}}_{I}` vector is not defined by thee user, but is defined by the
    relative position vector :math:`\\mathbf{r}_{C}` or velocity vector :math:`\\mathbf{r}_{C}` of the vehicle w.r.t. some body C. The inputs to this function allow :math:`\\hat{\\mathbf{T}}_{I}` to
    be set to :math:`\\pm\\mathbf{r}_{C}` or :math:`\\pm\\mathbf{v}_{C}`, for any body C. It is typically used for simplified or preliminary thrust analyses.
    
    
    Parameters
    ----------
    central_body : str
        Name of central body w.r.t. which the position/velocity vector is to be computed
    is_colinear_with_velocity : bool
        Boolean defining whether :math:`\\hat{\\mathbf{T}}_{I}` is to be aligned with velocity (if true) or position (if false)
    direction_is_opposite_to_vector : bool
        Boolean defining whether :math:`\\hat{\\mathbf{T}}_{I}` is to be in the same direction as position/velocity (if false), or in the opposite direction (if true).
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    free_rotation_angle_function : Callable[[float], float], default = None
        Custom function provided by the user, which returns a value for the free rotation angle :math:`\\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
    Returns
    -------
    BodyFixedDirectionBasedRotationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model."""

def simple(base_frame: str, target_frame: str, initial_orientation: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)], initial_time: float, rotation_rate: float) -> RotationModelSettings:
    """Function for creating simple rotation model settings.
    
    Function for settings object, defining a basic rotation model with constant orientation of the rotation axis and constant rotation rate about this axis.
    Rotation from original (inertial) to target (body-fixed) frame at some reference time ``initial_time`` (:math:`t_{0}`) is defined by the ``initial_orientation`` (:math:`\\mathbf{R}^{(B/I)}(t_{0})`) rotation matrix.
    Rotation about the body-fixed z-axis is defined by the ``rotation_rate`` (:math:`\\omega`) float variable (in rad/s). The rotation matrix is computed from:
    
    .. math::
       \\mathbf{R}^{(B/I)}(t)=\\mathbf{R}_{z}(\\omega(t-t_{0}))(t_{0})\\mathbf{R}^{(B/I)}(t_{0})
    
    where :math:`\\mathbf{R}^{(B/I)}` denotes the rotation matrix from inertial to body-fixed frame, and :math:`\\mathbf{R}_{z}` denotes a rotation matrix about the z-axis.
    
    The matrix :math:`\\mathbf{R}^{(B/I)}(t_{0})` is sometimes parameterized by pole right ascension and declination (:math:`\\alpha` and :math:`\\delta`), as well as the meridian of date :math:`W_{0}` with
    
    .. math::
       \\mathbf{R}^{(B/I)}(t_{0})=\\mathbf{R}_{z}(W_{0})\\mathbf{R}_{x}(\\pi/2-\\delta)\\mathbf{R}_{z}(\\pi/2+\\alpha)
    
    
    Parameters
    ----------
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    initial_orientation : numpy.ndarray[numpy.float64[3, 3]]
        Orientation of target frame in base frame at initial time.
    initial_time : float
        Initial time (reference epoch for rotation matrices).
    rotation_rate : float
        Constant rotation rate [rad/s] about rotational axis.
    Returns
    -------
    SimpleRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.SimpleRotationModelSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` for Earth,
    using a simple rotation model with constant orientation of the rotation axis (body-fixed z-axis), and constant rotation rate about this axis:
    
    .. code-block:: python
    
      # Set parameters describing the rotation between the two frames
      initial_orientation = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
      initial_time = 12345 # [sec since J2000]
      rotation_rate = 2e-5 # [rad/s]
      original_frame = "J2000"
      target_frame = "Earth_Fixed_Simplified"
      # Create the rotation model settings and assign to body settings of "Earth"
      body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.simple(
          original_frame,
          target_frame,
          initial_orientation,
          initial_time,
          rotation_rate)"""

def simple_from_spice(base_frame: str, target_frame: str, target_frame_spice: str, initial_time: float) -> RotationModelSettings:
    """Function for creating simple rotation model settings using initial orientation and rotation rates from Spice.
    
    Function for settings object, defining a :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple` rotation model with the added functionality that the initial orientation and rotation rate are extracted from Spice, as opposed to provided manually.
    Note that `only` the initial orientation and rotation rate ( at the time defined by ``initial_time`` ) are extracted from Spice - for
    the full Spice rotation model see :func:`~tudatpy.dynamics.environment_setup.rotation_model.spice`.
    Also note the distinction between the ``target_frame`` and ``target_frame_spice`` parameters.
    
    
    Parameters
    ----------
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Target frame of rotation model - name of frame that Tudat assigns to the body-fixed frame
    target_frame_spice : str
        Spice reference of target frame - name of the frame in Spice for which the initial orientation and rotation rate are extracted.
    initial_time : float
        Initial time (reference epoch for rotation matrices).
    Returns
    -------
    SimpleRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.SimpleRotationModelSettings` class
    
    
    
    Notes
    -----
    In order to create a :class:`~tudatpy.dynamics.environment_setup.rotation_model.SimpleRotationModelSettings` object which describes a synchronous rotation w.r.t. some ``central_body``,
    we require an ``ephemeris_settings`` attribute to the :class:`~tudatpy.dynamics.environment_setup.BodySettings` object of the ``central_body``.
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` for Earth,
    using a simple rotation model with constant orientation of the rotation axis (body-fixed z-axis), and constant rotation rate about this axis.
    The initial orientation and rotation rate are extracted from Spice at the time defined by ``initial_time``:
    
    .. code-block:: python
    
       # set parameters for time at which initial data is extracted from spice
       initial_time = 12345
       # set parameters for defining the rotation between frames
       original_frame = "J2000"
       target_frame = "IAU_Earth_Simplified"
       target_frame_spice = "IAU_Earth"
       # create rotation model settings and assign to body settings of "Earth"
       body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.simple_from_spice(
       original_frame, target_frame, target_frame_spice, initial_time)"""

def spice(base_frame: str, target_frame: str, spice_frame_name: str='') -> RotationModelSettings:
    """Function for creating rotation model settings from the Spice interface.
    
    Function for settings object, defining a rotation model directly (and entirely) from Spice interface.
    
    
    Parameters
    ----------
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    spice_frame_name : str, default = ""
        Name of the spice reference frame name that will be used to compute the rotation to the target frame. For instance, if target_frame is set to "IAU_Earth", and ``spice_frame_name`` is set to "IAU_Mars", Tudat will extract the rotation to the IAU_Mars frame from Spice, and assign this rotation to the "IAU_Earth" frame in Tudat. By default, this input is left empty, which corresponds to it being equal to the ``target_frame``.
    Returns
    -------
    RotationModelSettings
        Instance of :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` class.
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` for Earth,
    using full rotation model data from Spice:
    
    .. code-block:: python
    
       # define parameters describing the rotation between frames
       original_frame = "J2000"
       target_frame = "IAU_Earth"
       # create rotation model settings and assign to body settings of "Earth"
       body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.spice(
       original_frame, target_frame)"""

def synchronous(central_body_name: str, base_frame: str, target_frame: str) -> RotationModelSettings:
    """Function for creating synchronous rotational ephemeris settings.
    
    Function for settings object, defining a synchronous rotation model where rotation of a body is defined from its relative orbit w.r.t. some central body. Specifically
    - the body-fixed x-axis is *always* pointing towards the central body
    - the body-fixed z-axis is *always* perpendicular to the orbital plane (along the direction of :math:`\\mathbf{x}\\times\\mathbf{v}` )
    - the body-fixed y-axis completes the right-handed reference frame
    
    Such a model can be useful for, for instance, approximate rotation of tidally locked natural satellites or nadir-pointing spacecraft.
    
    
    Parameters
    ----------
    central_body_name : str
        Name of the central body of synchronous rotation.
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Spice reference of target frame.
    Returns
    -------
    SynchronousRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.SynchronousRotationModelSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` for the martian moon Phobos,
    We do so by assigning a synchronous rotation model to the rotation model settings of Phobos, using in this case ``"ECLIPJ2000"`` as the base frame,
    and ``"Phobos_Fixed"`` as the target frame.
    
    .. code-block:: python
    
       # define parameters describing the synchronous rotation model
       central_body_name = "Mars"
       original_frame = "ECLIPJ2000"
       target_frame = "Phobos_Fixed"
       # create rotation model settings for target frame and assign to body settings of "Phobos"
       body_settings.get( "Phobos" ).rotation_model_settings = environment_setup.rotation_model.synchronous(
       central_body_name, original_frame, target_frame)"""

def zero_pitch_moment_aerodynamic_angle_based(central_body: str, base_frame: str, target_frame: str, angle_funcion: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(2, 1)]]=None) -> RotationModelSettings:
    """Function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank.
    
    Function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank. This function is
    largely identical to the :func:`~aerodynamic_angle_based`, with the difference that the angle of attack :math:`\\alpha` is not provided as a custom value by the user, but is
    calculated from the body's aerodynamic moment coefficients, such that we have :math:`C_{m}=0`. This requires aerodynamic moment coefficients to be defined for the vehicle that
    depend on (among others) the body's angle of attack
    
    
    Parameters
    ----------
    central_body : str
        Name of the central body C that is to be used.
    base_frame : str
        Name of the base frame of rotation model.
    target_frame : str
        Name of the target frame of rotation model.
    angle_funcion : Callable[[float], numpy.ndarray[numpy.float64[2, 1]]], default = None
        Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\\beta,\\sigma]`. If this input is left empty, these angles are both fixed to 0.
    Returns
    -------
    CustomRotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model."""
gcrs_to_itrs_rotation_model: RotationModelType
iau_2000_a: IAUConventions
iau_2000_b: IAUConventions
iau_2006: IAUConventions
planetary_rotation_model: RotationModelType
simple_rotational_model: RotationModelType
spice_rotation_model: RotationModelType
synchronous_rotation_model: RotationModelType