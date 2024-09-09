import numpy
import typing
__all__ = ['KeplerianElementIndices', 'PositionElementTypes', 'SphericalOrbitalStateElementIndices', 'argument_of_periapsis_index', 'cartesian_position_type', 'cartesian_to_keplerian', 'cartesian_to_mee', 'cartesian_to_mee_manual_singularity', 'cartesian_to_spherical', 'convert_position_elements', 'delta_mean_anomaly_to_elapsed_time', 'eccentric_to_mean_anomaly', 'eccentric_to_true_anomaly', 'eccentricity_index', 'eclipj2000_state_to_teme', 'elapsed_time_to_delta_mean_anomaly', 'flight_path_index', 'flip_mee_singularity', 'geodetic_position_type', 'heading_angle_index', 'inclination_index', 'j2000_state_to_teme', 'keplerian_to_cartesian', 'keplerian_to_cartesian_elementwise', 'keplerian_to_mee', 'keplerian_to_mee_manual_singularity', 'latitude_index', 'longitude_index', 'longitude_of_ascending_node_index', 'mean_motion_to_semi_major_axis', 'mean_to_eccentric_anomaly', 'mean_to_true_anomaly', 'mee_to_cartesian', 'mee_to_keplerian', 'quaternion_entries_to_rotation_matrix', 'radius_index', 'rotation_matrix_to_quaternion_entries', 'semi_latus_rectum_index', 'semi_major_axis_index', 'semi_major_axis_to_mean_motion', 'speed_index', 'spherical_position_type', 'spherical_to_cartesian', 'spherical_to_cartesian_elementwise', 'teme_state_to_eclipj2000', 'teme_state_to_j2000', 'true_anomaly_index', 'true_to_eccentric_anomaly', 'true_to_mean_anomaly']

class KeplerianElementIndices:
    """Members:
	
	semi_major_axis_index
	
	eccentricity_index
	
	inclination_index
	
	argument_of_periapsis_index
	
	longitude_of_ascending_node_index
	
	true_anomaly_index
	
	semi_latus_rectum_index
	"""
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
    """Members:
	
	cartesian_position_type
	
	spherical_position_type
	
	geodetic_position_type
	"""
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
    """Members:
	
	radius_index
	
	latitude_index
	
	longitude_index
	
	speed_index
	
	flight_path_index
	
	heading_angle_index
	"""
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

def cartesian_to_keplerian(cartesian_elements: numpy.ndarray, gravitational_parameter: float) -> numpy.ndarray:
    """Convert Cartesian to Keplerian elements.
	
	.. note:: See module level documentation for the standard ordering
			  convention of Keplerian elements used.
	
	
	:param cartesian_elements:
			Cartesian state that is to be converted to Keplerian elements
	:param gravitational_parameter:
			Gravitational parameter of central body used for conversion
	:return:
			Keplerian elements, as computed from Cartesian element input.
	"""

def cartesian_to_mee(cartesian_elements: numpy.ndarray, gravitational_parameter: float) -> numpy.ndarray:
    """Convert Cartesian to Modified equinoctial elements.
	
	Convery cartesian to Modified equinoctial elements. The singularity-flipping
	element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)
	.. note:: See module level documentation for the standard ordering
			  convention of Modified Equinoctial elements used.
	
	
	:param cartesian_elements:
			Cartesian elements that are to be converted to Modified equinoctial elements
	:param gravitational_parameter:
			Gravitational parameter of central body
	:return:
			Modified equinoctial elements, as computed from Cartesian element input.
	"""

def cartesian_to_mee_manual_singularity(cartesian_elements: numpy.ndarray, gravitational_parameter: float, singularity_at_zero_inclination: bool) -> numpy.ndarray:
    """Convert Cartesian to Modified equinoctial elements.
	
	Convery cartesian to Modified equinoctial elements. The singularity-flipping
	element :math:`I` is to be provided manually for this function
	.. note:: See module level documentation for the standard ordering
			  convention of Modified Equinoctial elements used.
	
	
	:param cartesian_elements:
			Cartesian elements that are to be converted to Modified equinoctial elements
	:param gravitational_parameter:
			Gravitational parameter of central body
	:param singularity_at_zero_inclination:
			Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
			Modified equinoctial elements, as computed from Cartesian element input.
	"""

def cartesian_to_spherical(cartesian_elements: numpy.ndarray) -> numpy.ndarray:
    """Convert Cartesian to spherical elements.
	
	.. note:: See module level documentation for the standard ordering
			  convention of spherical state elements used.
	
	
	:param cartesian_elements:
			Cartesian state that is to be converted to spherical elements
	:return:
			Spherial elements, as computed from Cartesian element input.
	"""

def convert_position_elements(originalElements: numpy.ndarray, original_elemet_types: PositionElementTypes, new_element_types: PositionElementTypes, shape_model: ..., tolerance: float) -> numpy.ndarray:
    ...

def delta_mean_anomaly_to_elapsed_time(mean_anomaly_change: float, gravitational_parameter: float, semi_major_axis: float) -> float:
    """Convert change in mean anomaly along a Keplerian orbit to the corresponding elapsed time.
	
	:param mean_anomaly_change:
			Total change in mean anomaly along the Kepler orbit
	:param gravitational_parameter:
			Gravitational parameter of central body
	:param semi_major_axis:
			Semi-major axis of orbit
	:return:
			Time required for the provided mean anomaly change to be accumulated
	"""

def eccentric_to_mean_anomaly(eccentric_anomaly: float, eccentricity: float) -> float:
    """Convert eccentric to mean anomaly.
	
	:param eccentricity:
			Value of the orbital eccentricity
	:param eccentric_anomaly:
			Hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1
	:return:
			Mean of the true anomaly
	"""

def eccentric_to_true_anomaly(eccentric_anomaly: float, eccentricity: float) -> float:
    """Convert eccentric to true anomaly.
	
	:param eccentricity:
			Value of the orbital eccentricity
	:param eccentric_anomaly:
			Hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1
	:return:
			Value of the true anomaly
	"""

def eclipj2000_state_to_teme(epoch: float, eclipj2000_state: numpy.ndarray) -> numpy.ndarray:
    ...

def elapsed_time_to_delta_mean_anomaly(elapsed_time: float, gravitational_parameter: float, semi_major_axis: float) -> float:
    """Convert elapsed time to the corresponding change in mean anomaly along a Keplerian orbit.
	
	:param elapsed_time:
			Elapsed time (in seconds)
	:param gravitational_parameter:
			Gravitational parameter of central body
	:param semi_major_axis:
			Semi-major axis of orbit
	:return:
			Total change in mean anomaly along the Kepler orbit, accumulated in the provided time.
	"""

def flip_mee_singularity(keplerian_elements: numpy.ndarray) -> bool:
    """Function to determine 'optimal' location of the singularity-flipping modified equinoctial element.
	
	Function to determine 'optimal' location of the singularity-flipping modified equinoctial element :math:`I`, if orbit inclination is less than
	90 degrees, it puts the singularity at 180 degrees, if it is larger than 90 degrees, it puts it at 0 degrees.
	
	
	:param keplerian_elements:
			Keplerian elements that are to be converted to Modified equinoctial elements
	:return:
			Singularity at 0 degrees inclination if false, 180 degrees if true
	"""

def j2000_state_to_teme(epoch: float, j2000_state: numpy.ndarray) -> numpy.ndarray:
    ...

def keplerian_to_cartesian(keplerian_elements: numpy.ndarray, gravitational_parameter: float) -> numpy.ndarray:
    """Convert Keplerian elements to Cartesian.
	
	.. note:: See module level documentation for the standard ordering
			  convention of Keplerian elements used.
	
	
	:param keplerian_elements:
			Keplerian state that is to be converted to Cartesian elements
	:param gravitational_parameter:
			Gravitational parameter of central body used for conversion
	:return:
			Cartesian elements, as computed from Keplerian element input.
	"""

def keplerian_to_cartesian_elementwise(semi_major_axis: float, eccentricity: float, inclination: float, argument_of_periapsis: float, longitude_of_ascending_node: float, true_anomaly: float, gravitational_parameter: float) -> numpy.ndarray:
    """Convert Keplerian elements to Cartesian, with elementwise input.
	
	.. note:: The final Keplerian element is always the true anomaly.
	
	
	:param semi_major_axis:
			Semi-major axis (except if eccentricity = 1.0, then represents semi-latus rectum)
	:param eccentricity:
			Eccentricity
	:param inclination:
			Inclination
	:param argument_of_periapsis:
			Argument of periapsis
	:param longitude_of_ascending_node:
			Longitude of ascending node
	:param true_anomaly:
			True anomaly
	:param gravitational_parameter:
			Gravitational parameter of central body used for conversion
	:return:
			Cartesian elements, as computed from Keplerian element input.
	"""

def keplerian_to_mee(keplerian_elements: numpy.ndarray) -> numpy.ndarray:
    """Convert Keplerian to Modified equinoctial elements.
	
	Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
	element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)
	.. note:: See module level documentation for the standard ordering
			  convention of Modified Equinoctial elements used.
	
	
	:param keplerian_elements:
			Keplerian elements that are to be converted to Modified equinoctial elements
	:return:
			Modified equinoctial elements, as computed from Keplerian element input (with element :math:`I` defined by :func:`flip_mee_singularity`).
	"""

def keplerian_to_mee_manual_singularity(keplerian_elements: numpy.ndarray, singularity_at_zero_inclination: bool) -> numpy.ndarray:
    """Convert Keplerian to Modified equinoctial elements.
	
	Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
	element :math:`I` is to be provided manually for this function
	.. note:: See module level documentation for the standard ordering
			  convention of Modified Equinoctial elements used.
	
	
	:param keplerian_elements:
			Keplerian elements that are to be converted to Modified equinoctial elements
	:param singularity_at_zero_inclination:
			Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
			Modified equinoctial elements, as computed from Keplerian element input.
	"""

def mean_motion_to_semi_major_axis(mean_motion: float, gravitational_parameter: float) -> float:
    """Convert mean motion to corresponding semi-major axis (in a Keplerian orbit).
	
	:param mean_motion:
			Orbital mean motion
	:param gravitational_parameter:
			Gravitational parameter of central body
	:return:
			Semi-major axis corresponding to mean motion
	"""

def mean_to_eccentric_anomaly(eccentricity: float, mean_anomaly: float, use_default_initial_guess: bool=True, non_default_initial_guess: float=..., root_finder: ...=None) -> float:
    """Convert mean to eccentric anomaly.
	
	:param eccentricity:
			Value of the orbital eccentricity
	:param mean_anomaly:
			Value of the mean anomaly
	:param use_default_initial_guess:
			Boolean to determine whether the user-defined initial guess is used for conversion, or an automatically generated one.
	:param non_default_initial_guess:
			User-defined initial guess for conversion, to be used only if ``use_default_initial_guess`` is set to ``True``.
	:param root_finder:
			User-defined root finder, overriding default root-finding algorithm for conversion (default is used if this input is left empty)
	:return:
			Value of the eccentric anomaly
	"""

def mean_to_true_anomaly(eccentricity: float, mean_anomaly: float, use_default_initial_guess: bool=True, non_default_initial_guess: float=..., root_finder: ...=None) -> float:
    """Convert mean to true anomaly.
	
	Convert the mean anomaly of the orbit to its true anomaly. This conversion first converts mean to eccentric anomaly
	(hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1), and subsequently to true anomaly.
	
	
	:param eccentricity:
			Value of the orbital eccentricity
	:param mean_anomaly:
			Value of the mean anomaly
	:param use_default_initial_guess:
			Boolean to determine whether the user-defined initial guess (for mean-to-eccentric anomaly conversion) is used, or an automatically generated one.
	:param non_default_initial_guess:
			User-defined initial guess for mean-to-eccentric anomaly conversion, to be used only if ``use_default_initial_guess`` is set to ``True``.
	:param root_finder:
			User-defined root finder, overriding default root-finding algorithm for mean-to-eccentric anomaly conversion (default is used if this input is left empty)
	:return:
			Value of the true anomaly
	"""

def mee_to_cartesian(modified_equinoctial_elements: numpy.ndarray, gravitational_parameter: float, singularity_at_zero_inclination: bool) -> numpy.ndarray:
    """Convert Modified equinoctial to Cartesian elements.
	
	Convert Modified equinoctial to Cartesian elements
	.. note:: See module level documentation for the standard ordering
			  convention of Modified Equinoctial elements used.
	
	
	:param modified_equinoctial_elements:
			Modified equinoctial elements that are to be converted to Cartesian elements
	:param gravitational_parameter:
			Gravitational parameter of central body
	:param singularity_at_zero_inclination:
			Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
			Cartesian elements, as computed from Modified equinoctial element input.
	"""

def mee_to_keplerian(modified_equinoctial_elements: numpy.ndarray, singularity_at_zero_inclination: bool) -> numpy.ndarray:
    """Convert Modified equinoctial to Keplerian elements.
	
	Modified equinoctial elements to Keplerian (without intermediate step to Cartesian elements).
	.. note:: See module level documentation for the standard ordering
			  convention of Modified Equinoctial elements used.
	
	
	:param modified_equinoctial_elements:
			Modified equinoctial elements that are to be converted to Keplerian elements
	:param singularity_at_zero_inclination:
			Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
			Keplerian elements, as computed from Modified equinoctial element input.
	"""

def quaternion_entries_to_rotation_matrix(quaternion_entries: numpy.ndarray) -> numpy.ndarray:
    """Converts an array of four quaternion elements to the equivalent rotation matrix.
	
	Function to convert an array of four quaternion elements to the equivalent rotation matrix. These quaternion elements
	are for instance used when propagating rotational dynamics in Tudat, and this function can be used to convert the
	numerical results to a usable rotation matrix. See `our user guide <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/use_of_reference_frames.html#rotational-states>`_ for more details.
	
	
	:param quaternion_entries:
			Quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
	:return:
			Rotation matrix defining the equivalent rotation.
	"""

def rotation_matrix_to_quaternion_entries(rotation_matrix: numpy.ndarray) -> numpy.ndarray:
    """Converts a rotation matrix to the equivalent array of four quaternion elements.
	
	Inverse function of :func:`quaternion_entries_to_rotation_matrix`.
	
	
	:param rotation_matrix:
			Rotation matrix
	:return:
			Equivalent quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
	"""

def semi_major_axis_to_mean_motion(semi_major_axis: float, gravitational_parameter: float) -> float:
    """Convert semi-major axis to corresponding mean motion (along a Keplerian orbit).
	
	:param semi_major_axis:
			Semi-major axis of orbit
	:param gravitational_parameter:
			Gravitational parameter of central body
	:return:
			Semi-major axis corresponding to mean motion
	"""

def spherical_to_cartesian(spherical_elements: numpy.ndarray) -> numpy.ndarray:
    """Convert spherical elements to Cartesian.
	
	.. note:: See module level documentation for the standard ordering
			  convention of spherical state elements used.
	
	
	:param spherical_elements:
			Spherical state that is to be converted to Cartesian elements
	:return:
			Cartesian elements, as computed from spherical element input.
	"""

def spherical_to_cartesian_elementwise(radial_distance: float, latitude: float, longitude: float, speed: float, flight_path_angle: float, heading_angle: float) -> numpy.ndarray:
    """Convert Spherical elements to Cartesian, with elementwise input.
	
	:param radial_distance:
			Distance from origin of central body
	:param latitude:
			Central body-fixed latitude
	:param longitude:
			Central body-fixed longitude
	:param speed:
			Central body-fixed speed (norm of velocity vector). Note that this is *not* the norm of the inertial velocity
	:param flight_path_angle:
			Flight-path angle (of central body-fixed velocity vector)
	:param heading_angle:
			Heading angle (of central body-fixed velocity vector)
	:return:
			Cartesian elements, as computed from spherical element input.
	"""

def teme_state_to_eclipj2000(epoch: float, teme_state: numpy.ndarray) -> numpy.ndarray:
    ...

def teme_state_to_j2000(epoch: float, teme_state: numpy.ndarray) -> numpy.ndarray:
    ...

def true_to_eccentric_anomaly(true_anomaly: float, eccentricity: float) -> float:
    """Convert true to eccentric anomaly.
	
	:param eccentricity:
			Value of the orbital eccentricity
	:param true_anomaly:
			Value of the true anomaly
	:return:
			Hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1
	"""

def true_to_mean_anomaly(eccentricity: float, true_anomaly: float) -> float:
    """Convert true to mean anomaly.
	
	Convert the true anomaly of the orbit to its mean anomaly. This conversion first converts true to eccentric anomaly
	(hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1),
	and subsequently to mean anomaly.
	
	
	:param eccentricity:
			Value of the orbital eccentricity
	:param true_anomaly:
			Value of the true anomaly
	:return:
			Value of the mean anomaly
	"""
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