import numpy
import pybind11_stubgen.typing_ext
from ...math import root_finders
import typing
__all__ = ['KeplerianElementIndices', 'PositionElementTypes', 'SphericalOrbitalStateElementIndices', 'argument_of_periapsis_index', 'cartesian_position_type', 'cartesian_to_keplerian', 'cartesian_to_mee', 'cartesian_to_mee_manual_singularity', 'cartesian_to_spherical', 'cartesian_to_usm_6', 'cartesian_to_usm_7', 'cartesian_to_usm_em', 'convert_cartesian_to_geodetic_coordinates', 'convert_geographic_to_geodetic_latitude', 'convert_position_elements', 'delta_mean_anomaly_to_elapsed_time', 'eccentric_to_mean_anomaly', 'eccentric_to_true_anomaly', 'eccentricity_index', 'eclipj2000_to_j2000', 'elapsed_time_to_delta_mean_anomaly', 'exponential_map_to_quaternion', 'flight_path_index', 'flip_mee_singularity', 'geodetic_position_type', 'heading_angle_index', 'inclination_index', 'j2000_to_eclipj2000', 'j2000_to_teme', 'keplerian_to_cartesian', 'keplerian_to_cartesian_elementwise', 'keplerian_to_mee', 'keplerian_to_mee_manual_singularity', 'latitude_index', 'longitude_index', 'longitude_of_ascending_node_index', 'mean_motion_to_semi_major_axis', 'mean_to_eccentric_anomaly', 'mean_to_true_anomaly', 'mee_to_cartesian', 'mee_to_keplerian', 'modified_rodrigues_parameters_to_quaternion', 'quaternion_entries_to_rotation_matrix', 'quaternion_to_exponential_map', 'quaternion_to_modified_rodrigues_parameters', 'radius_index', 'rotation_matrix_to_quaternion_entries', 'semi_latus_rectum_index', 'semi_major_axis_index', 'semi_major_axis_to_mean_motion', 'speed_index', 'spherical_position_type', 'spherical_to_cartesian', 'spherical_to_cartesian_elementwise', 'teme_to_j2000', 'true_anomaly_index', 'true_to_eccentric_anomaly', 'true_to_mean_anomaly', 'usm_6_to_cartesian', 'usm_7_to_cartesian', 'usm_em_to_cartesian']

class KeplerianElementIndices:
    """Enumeration for indices of Keplerian elements
    
    
    
    Members:
    
      semi_major_axis_index : 
    
    Element 0 in vector of Keplerian elements (for eccentricity not equal to 1.0)
    
    
    
      semi_latus_rectum_index : 
    
    Element 0 in vector of Keplerian elements (for eccentricity equal to 1.0)
    
    
    
      eccentricity_index : 
    
    Element 1 in vector of Keplerian elements
    
    
    
      inclination_index : 
    
    Element 2 in vector of Keplerian elements
    
    
    
      argument_of_periapsis_index : 
    
    Element 3 in vector of Keplerian elements
    
    
    
      longitude_of_ascending_node_index : 
    
    Element 4 in vector of Keplerian elements
    
    
    
      true_anomaly_index : 
    
    Element 5 in vector of Keplerian elements"""
    __members__: typing.ClassVar[dict[str, KeplerianElementIndices]]
    argument_of_periapsis_index: typing.ClassVar[KeplerianElementIndices]
    eccentricity_index: typing.ClassVar[KeplerianElementIndices]
    inclination_index: typing.ClassVar[KeplerianElementIndices]
    longitude_of_ascending_node_index: typing.ClassVar[KeplerianElementIndices]
    semi_latus_rectum_index: typing.ClassVar[KeplerianElementIndices]
    semi_major_axis_index: typing.ClassVar[KeplerianElementIndices]
    true_anomaly_index: typing.ClassVar[KeplerianElementIndices]

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

class PositionElementTypes:
    """Enumeration describing different types of position element types (typically used for body-centered, body-fixed position)
          
    
    Members:
    
      cartesian_position_type
    
      spherical_position_type
    
      geodetic_position_type"""
    __members__: typing.ClassVar[dict[str, PositionElementTypes]]
    cartesian_position_type: typing.ClassVar[PositionElementTypes]
    geodetic_position_type: typing.ClassVar[PositionElementTypes]
    spherical_position_type: typing.ClassVar[PositionElementTypes]

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

class SphericalOrbitalStateElementIndices:
    """Enumeration for indices of spherical orbital state elements"
    
    
    
    Members:
    
      radius_index : 
    
    Element 0 in vector of spherical orbital state elements
    
    
    
      latitude_index : 
    
    Element 1 in vector of spherical orbital state elements
    
    
    
      longitude_index : 
    
    Element 2 in vector of spherical orbital state elements
    
    
    
      speed_index : 
    
    Element 3 in vector of spherical orbital state elements
    
    
    
      flight_path_index : 
    
    Element 4 in vector of spherical orbital state elements
    
    
    
      heading_angle_index : 
    
    Element 5 in vector of spherical orbital state elements"""
    __members__: typing.ClassVar[dict[str, SphericalOrbitalStateElementIndices]]
    flight_path_index: typing.ClassVar[SphericalOrbitalStateElementIndices]
    heading_angle_index: typing.ClassVar[SphericalOrbitalStateElementIndices]
    latitude_index: typing.ClassVar[SphericalOrbitalStateElementIndices]
    longitude_index: typing.ClassVar[SphericalOrbitalStateElementIndices]
    radius_index: typing.ClassVar[SphericalOrbitalStateElementIndices]
    speed_index: typing.ClassVar[SphericalOrbitalStateElementIndices]

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

def cartesian_to_keplerian(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Cartesian to Keplerian elements.
    
    .. note:: See module level documentation for the standard ordering
              convention of Keplerian elements used.
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian state that is to be converted to Keplerian elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        Keplerian elements, as computed from Cartesian element input."""

def cartesian_to_mee(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Cartesian to Modified equinoctial elements.
    
    Convert cartesian to Modified equinoctial elements. The singularity-flipping
    element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)
    
    .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian elements that are to be converted to Modified equinoctial elements
    gravitational_parameter : float
        Gravitational parameter of central body
    Returns
    -------
    numpy.ndarray
        Modified equinoctial elements, as computed from Cartesian element input."""

def cartesian_to_mee_manual_singularity(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float, singularity_at_zero_inclination: bool) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Cartesian to Modified equinoctial elements.
    
    Convert cartesian to Modified equinoctial elements. The singularity-flipping
    element :math:`I` is to be provided manually for this function
    
    .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian elements that are to be converted to Modified equinoctial elements
    gravitational_parameter : float
        Gravitational parameter of central body
    singularity_at_zero_inclination : bool
        Singularity at 0 degrees inclination if false, 180 degrees if true
    Returns
    -------
    numpy.ndarray
        Modified equinoctial elements, as computed from Cartesian element input."""

def cartesian_to_spherical(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Cartesian to spherical elements.
    
    .. note:: See module level documentation for the standard ordering  convention of spherical state elements used.
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian state that is to be converted to spherical elements
    Returns
    -------
    numpy.ndarray
        Spherical elements, as computed from Cartesian element input."""

def cartesian_to_usm_6(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(7, 1)]:
    """Convert Cartesian elements to Unified State Model (USM) elements with Modified Rodrigues parameters map for rotational coordinates.
    
    .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian state that is to be converted to USM elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        USM elements using Modified Rodrigues parameters, as computed from Cartesian element input."""

def cartesian_to_usm_7(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(7, 1)]:
    """Convert Cartesian elements to Unified State Model (USM) elements with quaternion for rotational coordinates.
    
    .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian state that is to be converted to USM elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        USM elements using quaternion, as computed from Cartesian element input."""

def cartesian_to_usm_em(cartesian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(7, 1)]:
    """Convert Cartesian elements to Unified State Model (USM) elements with Exponential map for rotational coordinates.
    
    .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements
    
    
    Parameters
    ----------
    cartesian_elements : numpy.ndarray
        Cartesian state that is to be converted to USM elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        USM elements using exponential map, as computed from Cartesian element input."""

def convert_cartesian_to_geodetic_coordinates(cartesian_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], equatorial_radius: float, flattening: float, tolerance: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
    """Convert cartesian position to geodetic position, given equatorial radius and flattening (based on the Body Shape Model).
    
    Parameters
    ----------
    cartesian_coordinates : numpy.ndarray
        Position from which the conversion is to be performed
    equatorial_radius : float
    flattening : float
    tolerance : float
        Tolerance (in meters) used as convergence criterion for converting to/from geodetic altitude
    Returns
    -------
    numpy.ndarray
        Geodetic coordinates, as computed from Cartesian element input."""

def convert_geographic_to_geodetic_latitude(geographic_latitude: float, equatorial_radius: float, flattening: float, altitude: float, tolerance: float, maximum_number_of_iterations: int) -> float:
    """Convert geographic latitude to geodetic latitude.
    
    This function converts a geographic latitude to a geodetic latitude, given the equatorial radius, flattening, altitude, 
    and a convergence tolerance. The conversion is performed iteratively, with a maximum number of iterations allowed.
    
    Parameters
    ----------
    geographic_latitude : float
        Geographic latitude (in radians) to be converted.
    equatorial_radius : float
        Equatorial radius of the reference ellipsoid (in meters).
    flattening : float
        Flattening of the reference ellipsoid.
    altitude : float
        Altitude above the reference ellipsoid (in meters).
    tolerance : float
        Convergence tolerance for the iterative conversion (in radians).
    maximum_number_of_iterations : int
        Maximum number of iterations allowed for the conversion process.
    
    Returns
    -------
    float
        Geodetic latitude (in radians), as computed from the geographic latitude input.
    
    Raises
    ------
    RuntimeError
        If the maximum number of iterations is exceeded without achieving the desired tolerance."""

def convert_position_elements(original_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], original_element_type: PositionElementTypes, new_element_type: PositionElementTypes, shape_model: ..., tolerance: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
    """Convert position from one element set to another
    
    Parameters
    ----------
    original_elements : numpy.ndarray
        Position from which the conversion is to be performed
    original_element_type : PositionElementTypes
        Element type in which ``original_elements`` is provided
    new_element_type : PositionElementTypes
        Element type to which the ``original_elements`` are to be converted
    shape_model : BodyShapeModel
        Shape model of body sed to transform altitudes (w.r.t. this shape)
        to/from distances from the center of the body (can be set to ``None`` if requested conversion does not require this model)
    tolerance : float
        Tolerance (in meters) used as convergence criterion for converting to/from geodetic altitude
    Returns
    -------
    numpy.ndarray
        Keplerian elements, as computed from Cartesian element input."""

def delta_mean_anomaly_to_elapsed_time(mean_anomaly_change: float, gravitational_parameter: float, semi_major_axis: float) -> float:
    """Convert change in mean anomaly along a Keplerian orbit to the corresponding elapsed time.
    
    
    Parameters
    ----------
    mean_anomaly_change : float
        Total change in mean anomaly along the Kepler orbit
    gravitational_parameter : float
        Gravitational parameter of central body
    semi_major_axis : float
        Semi-major axis of orbit
    Returns
    -------
    float
        Time required for the provided mean anomaly change to be accumulated"""

def eccentric_to_mean_anomaly(eccentric_anomaly: float, eccentricity: float) -> float:
    """Convert eccentric to mean anomaly.
    
    
    Parameters
    ----------
    eccentric_anomaly : float
        Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1
    eccentricity : float
        Value of the orbital eccentricity
    Returns
    -------
    float
        Value of the mean anomaly"""

def eccentric_to_true_anomaly(eccentric_anomaly: float, eccentricity: float) -> float:
    """Convert eccentric to true anomaly.
    
    
    Parameters
    ----------
    eccentric_anomaly : float
        Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1
    eccentricity : float
        Value of the orbital eccentricity
    Returns
    -------
    float
        Value of the true anomaly"""

def eclipj2000_to_j2000() -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Provides the (constant) rotation matrix from the ECLIPJ2000 to the J2000 frame, as defined in the SPICE library (see :ref:`spice` for more details on our interface with this library).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from ECLIPJ2000 to J2000 frame"""

def elapsed_time_to_delta_mean_anomaly(elapsed_time: float, gravitational_parameter: float, semi_major_axis: float) -> float:
    """Convert elapsed time to the corresponding change in mean anomaly along a Keplerian orbit.
    
    
    Parameters
    ----------
    elapsed_time : float
        Elapsed time (in seconds)
    gravitational_parameter : float
        Gravitational parameter of central body
    semi_major_axis : float
        Semi-major axis of orbit
    Returns
    -------
    float
        Total change in mean anomaly along the Kepler orbit, accumulated in the provided time."""

def exponential_map_to_quaternion(exponential_map: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]:
    """Converts modified Rodrigues parameters to the equivalent exponential map (for rotation representation).
    
    Parameters
    ----------
    exponential_map : numpy.ndarray
        Exponential map rotation elements
    Returns
    -------
    numpy.ndarray
        Equivalent quaternion elements"""

def flip_mee_singularity(keplerian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> bool:
    """Function to determine 'optimal' location of the singularity-flipping modified equinoctial element.
    
    Function to determine 'optimal' location of the singularity-flipping modified equinoctial element :math:`I`, if orbit inclination is less than
    90 degrees, it puts the singularity at 180 degrees, if it is larger than 90 degrees, it puts it at 0 degrees.
    
    
    Parameters
    ----------
    keplerian_elements : numpy.ndarray
        Keplerian elements that are to be converted to Modified equinoctial elements
    Returns
    -------
    bool
        Singularity at 0 degrees inclination if false, 180 degrees if true"""

def j2000_to_eclipj2000() -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Provides the (constant) rotation matrix from the J2000 to the ECLIPJ2000 frame, as defined in the SPICE library (see :ref:`spice` for more details on our interface with this library).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from J2000 to ECLIPJ2000 frame"""

def j2000_to_teme(epoch: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from the J2000 to the TEME (True Equator Mean Equinox) frame, which is the inverse of the :func:`~teme_to_j2000` function.
    
    Parameters
    ----------
    epoch : float
        Time as TDB seconds since J2000
    Returns
    -------
    numpy.ndarray
        Rotation matrix from J2000 to TEME frame"""

def keplerian_to_cartesian(keplerian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Keplerian elements to Cartesian.
    
    .. note:: See module level documentation for the standard ordering
              convention of Keplerian elements used.
    
    
    Parameters
    ----------
    keplerian_elements : numpy.ndarray
        Keplerian state that is to be converted to Cartesian elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from Keplerian element input."""

def keplerian_to_cartesian_elementwise(semi_major_axis: float, eccentricity: float, inclination: float, argument_of_periapsis: float, longitude_of_ascending_node: float, true_anomaly: float, gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Keplerian elements to Cartesian, with elementwise input.
    
    .. note:: The final Keplerian element is always the true anomaly.
    
    
    Parameters
    ----------
    semi_major_axis : float
        Semi-major axis (except if eccentricity = 1.0, then represents semi-latus rectum)
    eccentricity : float
        Eccentricity
    inclination : float
        Inclination
    argument_of_periapsis : float
        Argument of periapsis
    longitude_of_ascending_node : float
        Longitude of ascending node
    true_anomaly : float
        True anomaly
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from Keplerian element input."""

def keplerian_to_mee(keplerian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Keplerian to Modified equinoctial elements.
    
    Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
    element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)
    
    .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.
    
    
    Parameters
    ----------
    keplerian_elements : numpy.ndarray
        Keplerian elements that are to be converted to Modified equinoctial elements
    Returns
    -------
    numpy.ndarray
        Modified equinoctial elements, as computed from Keplerian element input (with element :math:`I` defined by :func:`flip_mee_singularity`)."""

def keplerian_to_mee_manual_singularity(keplerian_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], singularity_at_zero_inclination: bool) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Keplerian to Modified equinoctial elements.
    
    Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
    element :math:`I` is to be provided manually for this function
    
    .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.
    
    
    Parameters
    ----------
    keplerian_elements : numpy.ndarray
        Keplerian elements that are to be converted to Modified equinoctial elements
    singularity_at_zero_inclination : bool
        Singularity at 0 degrees inclination if ``True``, 180 degrees if ``False``
    Returns
    -------
    numpy.ndarray
        Modified equinoctial elements, as computed from Keplerian element input."""

def mean_motion_to_semi_major_axis(mean_motion: float, gravitational_parameter: float) -> float:
    """Convert mean motion to corresponding semi-major axis (in a Keplerian orbit).
    
    
    Parameters
    ----------
    mean_motion : float
        Orbital mean motion
    gravitational_parameter : float
        Gravitational parameter of central body
    Returns
    -------
    float
        Semi-major axis corresponding to mean motion"""

def mean_to_eccentric_anomaly(eccentricity: float, mean_anomaly: float, use_default_initial_guess: bool=True, non_default_initial_guess: float=..., root_finder: root_finders.RootFinderCore=None) -> float:
    """Convert mean to eccentric anomaly.
    
    
    Parameters
    ----------
    eccentricity : float
        Value of the orbital eccentricity
    mean_anomaly : float
        Value of the mean anomaly
    use_default_initial_guess : bool, default = True
        Boolean to determine whether the user-defined initial guess is used for conversion, or an automatically generated one.
    non_default_initial_guess : float, default = NaN
        User-defined initial guess for conversion, to be used only if ``use_default_initial_guess`` is set to ``False``.
    root_finder : RootFinder, default = None
        User-defined root finder, overriding default root-finding algorithm for conversion (default is used if this input is left empty)
    Returns
    -------
    float
        Value of the eccentric anomaly"""

def mean_to_true_anomaly(eccentricity: float, mean_anomaly: float, use_default_initial_guess: bool=True, non_default_initial_guess: float=..., root_finder: root_finders.RootFinderCore=None) -> float:
    """Convert mean to true anomaly.
    
    Convert the mean anomaly of the orbit to its true anomaly. This conversion first converts mean to eccentric anomaly
    (hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1), and subsequently to true anomaly.
    
    
    Parameters
    ----------
    eccentricity : float
        Value of the orbital eccentricity
    mean_anomaly : float
        Value of the mean anomaly
    use_default_initial_guess : bool, default = True
        Boolean to determine whether the user-defined initial guess (for mean-to-eccentric anomaly conversion) is used, or an automatically generated one.
    non_default_initial_guess : float, default = NaN
        User-defined initial guess for mean-to-eccentric anomaly conversion, to be used only if ``use_default_initial_guess`` is set to ``False``.
    root_finder : RootFinder, default = None
        User-defined root finder, overriding default root-finding algorithm for mean-to-eccentric anomaly conversion (default is used if this input is left empty)
    Returns
    -------
    float
        Value of the true anomaly"""

def mee_to_cartesian(modified_equinoctial_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], gravitational_parameter: float, singularity_at_zero_inclination: bool) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Modified equinoctial to Cartesian elements.
    
    
    .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.
    
    
    Parameters
    ----------
    modified_equinoctial_elements : numpy.ndarray
        Modified equinoctial elements that are to be converted to Cartesian elements
    gravitational_parameter : float
        Gravitational parameter of central body
    singularity_at_zero_inclination : bool
        Singularity at 0 degrees inclination if false, 180 degrees if true
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from Modified equinoctial element input."""

def mee_to_keplerian(modified_equinoctial_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], singularity_at_zero_inclination: bool) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Modified equinoctial to Keplerian elements.
    
    Modified equinoctial elements to Keplerian (without intermediate step to Cartesian elements).
    
    .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.
    
    
    Parameters
    ----------
    modified_equinoctial_elements : numpy.ndarray
        Modified equinoctial elements that are to be converted to Keplerian elements
    singularity_at_zero_inclination : bool
        Singularity at 0 degrees inclination if false, 180 degrees if true
    Returns
    -------
    numpy.ndarray
        Keplerian elements, as computed from Modified equinoctial element input."""

def modified_rodrigues_parameters_to_quaternion(modified_rodrigues_parameters: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]:
    """Converts modified Rodrigues parameters to the equivalent array of four quaternion elements (for rotation representation).
    
    Parameters
    ----------
    modified_rodrigues_parameters : numpy.ndarray
        Modified Rodrigues parameters
    Returns
    -------
    numpy.ndarray
        Equivalent quaternion elements"""

def quaternion_entries_to_rotation_matrix(quaternion_entries: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Converts an array of four quaternion elements to the equivalent rotation matrix.
    
    Function to convert an array of four quaternion elements to the equivalent rotation matrix. These quaternion elements
    are for instance used when propagating rotational dynamics in Tudat, and this function can be used to convert the
    numerical results to a usable rotation matrix. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html?highlight=rotational%20states#rotational-states>`_ for more details.
    
    
    Parameters
    ----------
    quaternion_entries : numpy.ndarray
        Quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
    Returns
    -------
    numpy.ndarray
        Rotation matrix defining the equivalent rotation."""

def quaternion_to_exponential_map(quaternion_entries: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]:
    """Converts quaternion elements to the equivalent exponential map (for rotation representation).
    
    Parameters
    ----------
    quaternion_entries : numpy.ndarray
        Quaternion elements
    Returns
    -------
    numpy.ndarray
        Equivalent exponential map rotation elements"""

def quaternion_to_modified_rodrigues_parameters(quaternion_entries: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]:
    """Converts quaternion elements to the equivalent modified Rodrigues parameters (for rotation representation).
    
    Parameters
    ----------
    quaternion_entries : numpy.ndarray
        Quaternion elements
    Returns
    -------
    numpy.ndarray
        Equivalent modified Rodrigues parameters"""

def rotation_matrix_to_quaternion_entries(rotation_matrix: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(4, 1)]:
    """Converts a rotation matrix to the equivalent array of four quaternion elements.
    
    Inverse function of :func:`quaternion_entries_to_rotation_matrix`.
    
    
    Parameters
    ----------
    rotation_matrix : numpy.ndarray
        Rotation matrix
    Returns
    -------
    numpy.ndarray
        Equivalent quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_"""

def semi_major_axis_to_mean_motion(semi_major_axis: float, gravitational_parameter: float) -> float:
    """Convert semi-major axis to corresponding mean motion (along a Keplerian orbit).
    
    
    Parameters
    ----------
    semi_major_axis : float
        Semi-major axis of orbit
    gravitational_parameter : float
        Gravitational parameter of central body
    Returns
    -------
    float
        Semi-major axis corresponding to mean motion"""

def spherical_to_cartesian(spherical_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert spherical elements to Cartesian.
    
    .. note:: See module level documentation for the standard ordering convention of spherical state elements used.
    
    
    Parameters
    ----------
    spherical_elements : numpy.ndarray
        Spherical state that is to be converted to Cartesian elements
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from spherical element input."""

def spherical_to_cartesian_elementwise(radial_distance: float, latitude: float, longitude: float, speed: float, flight_path_angle: float, heading_angle: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Spherical elements to Cartesian, with elementwise input.
    
    
    Parameters
    ----------
    radial_distance : float
        Distance from origin of central body
    latitude : float
        Central body-fixed latitude
    longitude : float
        Central body-fixed longitude
    speed : float
        Central body-fixed speed (norm of velocity vector). Note that this is *not* the norm of the inertial velocity
    flight_path_angle : float
        Flight-path angle (of central body-fixed velocity vector)
    heading_angle : float
        Heading angle (of central body-fixed velocity vector)
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from spherical element input."""

def teme_to_j2000(epoch: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from the TEME (True Equator Mean Equinox) frame to the J2000 frame, using the following:
    
    .. math::
        \\mathbf{R}^{(\\text{J2000}/\\text{TEME})}=\\mathbf{PN}(t)\\mathbf{R}_{z}(-\\theta(t)))
    
    where :math:`\\theta` is the difference between the actual and mean position of the first point of Aries (or 'equation of the equinoxes`), computes using the ``iauEe00b`` function of the SOFA
    library (which computes this angle compatible with IAU 2000 resolutions but using the truncated nutation model IAU 2000B), and the precession-nutation matrix :math:`\\mathbf{PN}` is computes using
    the function ``iauPnm80`` of the Sofa library, which uses th IAU 1976 precession model, and the IAU 1980 nutation model. The choice of slightly inconsistent IAU conventions is made for computational
    efficiency, as the 1976/1980 model includes fewer terms that the newer resolutions, and is combined with the truncated nutation model. Since the definition of the TEME frame is slightly ambiguous,
    the possible computational error that is incurred is insignificant for the application of the TEME frame: the transformation of SGP4-propagated TLEs to other epochs.
    
    Parameters
    ----------
    epoch : float
        Time as TDB seconds since J2000
    Returns
    -------
    numpy.ndarray
        Rotation matrix from TEME to J2000 frame"""

def true_to_eccentric_anomaly(true_anomaly: float, eccentricity: float) -> float:
    """Convert true to eccentric anomaly.
    
    
    Parameters
    ----------
    eccentricity : float
        Value of the orbital eccentricity
    true_anomaly : float
        Value of the true anomaly
    Returns
    -------
    float
        Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1"""

def true_to_mean_anomaly(eccentricity: float, true_anomaly: float) -> float:
    """Convert true to mean anomaly.
    
    Convert the true anomaly of the orbit to its mean anomaly. This conversion first converts true to eccentric anomaly
    (hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1),
    and subsequently to mean anomaly.
    
    
    Parameters
    ----------
    eccentricity : float
        Value of the orbital eccentricity
    true_anomaly : float
        Value of the true anomaly
    Returns
    -------
    float
        Value of the mean anomaly"""

def usm_6_to_cartesian(usm_6_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(7, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Unified State Model (USM) elements with Modified Rodrigues parameters for rotational coordinates to Cartesian elements.
    
    .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements
    
    
    Parameters
    ----------
    usm_6_elements : numpy.ndarray
        USM elements using Modified Rodrigues parameters that is to be converted to Cartesian elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from USM element input."""

def usm_7_to_cartesian(usm_7_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(7, 1)], gravitational_parameter: float, normalize_quaternion: bool=True) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Unified State Model (USM) elements with quaternion for rotational coordinates to Cartesian elements.
    
    .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements
    
    
    Parameters
    ----------
    usm_7_elements : numpy.ndarray
        USM elements using quaternion that is to be converted to Cartesian elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from USM element input."""

def usm_em_to_cartesian(usm_em_elements: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(7, 1)], gravitational_parameter: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Convert Unified State Model (USM) elements with Exponential map for rotational coordinates to Cartesian elements.
    
    .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements
    
    
    Parameters
    ----------
    usm_em_elements : numpy.ndarray
        USM elements using exponential map that is to be converted to Cartesian elements
    gravitational_parameter : float
        Gravitational parameter of central body used for conversion
    Returns
    -------
    numpy.ndarray
        Cartesian elements, as computed from USM element input."""
argument_of_periapsis_index: KeplerianElementIndices
cartesian_position_type: PositionElementTypes
eccentricity_index: KeplerianElementIndices
flight_path_index: SphericalOrbitalStateElementIndices
geodetic_position_type: PositionElementTypes
heading_angle_index: SphericalOrbitalStateElementIndices
inclination_index: KeplerianElementIndices
latitude_index: SphericalOrbitalStateElementIndices
longitude_index: SphericalOrbitalStateElementIndices
longitude_of_ascending_node_index: KeplerianElementIndices
radius_index: SphericalOrbitalStateElementIndices
semi_latus_rectum_index: KeplerianElementIndices
semi_major_axis_index: KeplerianElementIndices
speed_index: SphericalOrbitalStateElementIndices
spherical_position_type: PositionElementTypes
true_anomaly_index: KeplerianElementIndices