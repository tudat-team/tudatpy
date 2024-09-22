import typing
import numpy
__all__ = ['body_fixed_to_inertial_rotation_matrix', 'inertial_to_body_fixed_rotation_matrix', 'inertial_to_rsw_rotation_matrix', 'inertial_to_tnw_rotation_matrix', 'rsw_to_inertial_rotation_matrix', 'tnw_to_inertial_rotation_matrix', 'transform_cartesian_state_to_frame']

def body_fixed_to_inertial_rotation_matrix(pole_declination: float, pole_right_ascension: float, pole_meridian: float) -> numpy.ndarray:
    """Computes the rotation matrix from body-fixed to inertial frame.
	
	
	Function to compute the rotation matrix from body-fixed to inertial
	frame, using typical pole right ascension (:math:`\x07lpha`), pole
	declination (:math:`\\delta`), and prime meridian longitude
	(:math:`W`) angles.
	
	
	:param pole_declination:
			Declination of body pole in inertial frame (:math:`\\delta`).
	
	:param pole_right_ascension:
			Right ascension of body pole in inertial frame (:math:`\x07lpha`).
	
	:param prime_meridian_longitude:
			Longitude of prime meridian w.r.t. intermediate frame
			(:math:`W`).
	
	:return:
			Rotation matrix from body-fixed to inertial frame.
	"""

def inertial_to_body_fixed_rotation_matrix(pole_declination: float, pole_right_ascension: float, prime_meridian_longitude: float) -> numpy.ndarray:
    """Computes the rotation matrix from inertial to body-fixed frame.
	
	
	Function to compute the rotation matrix from inertial to body-fixed
	frame, using typical pole right ascension (:math:`\x07lpha`), pole
	declination (:math:`\\delta`), and prime meridian longitude
	(:math:`W`) angles.
	
	
	:param pole_declination:
			Declination of body pole in inertial frame (:math:`\\delta`).
	
	:param pole_right_ascension:
			Right ascension of body pole in inertial frame (:math:`\x07lpha`).
	
	:param prime_meridian_longitude:
			Longitude of prime meridian w.r.t. intermediate frame
			(:math:`W`).
	
	:return:
			Rotation matrix from inertial to body-fixed frame
	"""

def inertial_to_rsw_rotation_matrix(inertial_cartesian_state: numpy.ndarray) -> numpy.ndarray:
    """Computes the rotation matrix from inertial to RSW frame.
	
	
	Function to compute the rotation matrix from inertial to RSW frame.
	The RSW frame is defined  by the state of a body w.r.t. to some
	central body. The x-axis of the RSW frame points away from the
	origin, and the y-axis lies in the orbital plane, and is positive
	for in the direction of the velocity vector (but is not colinear
	with the velocity vector, except for circular orbits). The z-axis
	is perpendicular to the orbital plane, and completes the
	right-handed coordinate system.
	
	
	:param inertial_cartesian_state:
			Cartesian state, in an inertial frame, for which the rotation
			matrix is to be calculated. Note that the RSW frame is defined
			w.r.t. some central body, and this Cartesian state must be
			defined w.r.t. that central body (e.g. central body at the
			origin).
	
	:return:
			Rotation matrix from inertial to RSW frame.
	"""

def inertial_to_tnw_rotation_matrix(inertial_cartesian_state: numpy.ndarray, n_axis_points_away_from_central_body: bool=True) -> numpy.ndarray:
    """Computes the rotation matrix from inertial to TNW frame.
	
	
	Function to compute the rotation matrix from inertial to TNW frame.
	The TNW frame is defined by the state of a body w.r.t. to some
	central body. The x-axis of the TNW frame points along the velocity
	vector, and the y-axis lies in the orbital plane, and is positive
	in the direction away from the central body (or positive **towards**
	the central body if the ``n_axis_points_away_from_central_body``
	variable is set to false, see below). The z-axis is perpendicular
	to the orbital plane, and completes the right-handed coordinate
	system.
	
	
	:param inertial_cartesian_state:
			Cartesian state, in an inertial frame, for which the rotation
			matrix is to be calculated. Note that the RSW frame is defined
			w.r.t. some central body, and this Cartesian state must be
			defined w.r.t. that central body (e.g. central body at the
			origin).
	
	:param n_axis_points_away_from_central_body:
			Boolean (default is ``True``) defining whether the N axis of the
			TNW frame points away from the central body (if ``True``) or
			towards the central body (if ``False``).
	
	:return:
			Rotation matrix from inertial to TNW frame.
	"""

def rsw_to_inertial_rotation_matrix(inertial_cartesian_state: numpy.ndarray) -> numpy.ndarray:
    """Computes the rotation matrix from RSW to inertial frame.
	
	
	Function to compute the rotation matrix from RSW to inertial. The
	RSW frame is defined  by the state of a body w.r.t. to some central
	body. The x-axis of the RSW frame points away from the origin, and
	the y-axis lies in the orbital plane, and is positive for in the
	direction of the velocity vector (but is not colinear with the
	velocity vector, except for circular orbits). The z-axis is
	perpendicular to the orbital plane, and completes the right-handed
	coordinate system.
	
	
	:param inertial_cartesian_state:
			Cartesian state, in an inertial frame, for which the rotation
			matrix is to be calculated. Note that the RSW frame is defined
			w.r.t. some central body, and this Cartesian state must be
			defined w.r.t. that central body (e.g. central body at the
			origin).
	
	:return:
			Rotation matrix from RSW to inertial frame.
	"""

def tnw_to_inertial_rotation_matrix(inertial_cartesian_state: numpy.ndarray, n_axis_points_away_from_central_body: bool=True) -> numpy.ndarray:
    """Computes the rotation matrix from TNW to inertial frame.
	
	
	Function to compute the rotation matrix from TNW to inertial frame.
	The TNW frame is defined by the state of a body w.r.t. to some
	central body. The x-axis of the TNW frame points along the velocity
	vector, and the y-axis lies in the orbital plane, and is positive
	in the direction away from the central body (or positive **towards**
	the central body if the ``n_axis_points_away_from_central_body``
	variable is set to false, see below). The z-axis is perpendicular
	to the orbital plane, and completes the right-handed coordinate
	system.
	
	
	:param inertial_cartesian_state:
			Cartesian state, in an inertial frame, for which the rotation
			matrix is to be calculated. Note that the TNW frame is defined
			w.r.t. some central body, and this Cartesian state must be
			defined w.r.t. that central body (e.g. central body at the
			origin).
	
	:param n_axis_points_away_from_central_body:
			Boolean (default=``True``) defining whether the N axis of the
			TNW frame points away from the central body (if ``True``) or
			towards the central body (if ``False``).
	
	:return:
			Rotation matrix from TNW to inertial frame
	"""

def transform_cartesian_state_to_frame(original_state: numpy.ndarray, rotation_matrix: numpy.ndarray, transform_cartesian_state_to_frame: numpy.ndarray) -> numpy.ndarray:
    ...