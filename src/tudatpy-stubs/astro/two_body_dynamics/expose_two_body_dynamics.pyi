import typing
import numpy
import tudatpy.math.root_finders.expose_root_finders
__all__ = ['EccentricityFindingFunctions', 'LambertTargeter', 'LambertTargeterGooding', 'LambertTargeterIzzo', 'MultiRevolutionLambertTargeterIzzo', 'PericenterFindingFunctions', 'ZeroRevolutionLambertTargeterIzzo', 'compute_escape_or_capture_delta_v', 'propagate_kepler_orbit']

class EccentricityFindingFunctions:

    def __init__(self, absolute_incoming_semi_major_axis: float, absolute_outgoing_semi_major_axis: float, bending_angle: float) -> None:
        ...

    def compute_derivative_incoming_eccentricity_fn(self, arg0: float) -> float:
        ...

    def compute_incoming_eccentricity_fn(self, arg0: float) -> float:
        ...

class LambertTargeter:

    def __init__(self, departure_position: numpy.ndarray, arrival_position: numpy.ndarray, time_of_flight: float, gravitational_parameter: float) -> None:
        ...

    def get_arrival_velocity(self) -> numpy.ndarray:
        ...

    def get_departure_velocity(self) -> numpy.ndarray:
        ...

    def get_velocity_vectors(self) -> tuple[numpy.ndarray, numpy.ndarray]:
        ...

class LambertTargeterGooding(LambertTargeter):

    def __init__(self, departure_position: numpy.ndarray, arrival_position: numpy.ndarray, time_of_flight: float, gravitational_parameter: float, root_finder: tudatpy.math.root_finders.expose_root_finders.RootFinderCore=None) -> None:
        ...

    def get_radial_arrival_velocity(self) -> float:
        ...

    def get_radial_departure_velocity(self) -> float:
        ...

    def get_semi_major_axis(self) -> float:
        ...

    def get_transverse_arrival_velocity(self) -> float:
        ...

    def get_transverse_departure_velocity(self) -> float:
        ...

class LambertTargeterIzzo(LambertTargeter):

    def __init__(self, departure_position: numpy.ndarray, arrival_position: numpy.ndarray, time_of_flight: float, gravitational_parameter: float, is_retrograde: bool=False, tolerance: float=1e-09, max_iter: int=50) -> None:
        ...

    def get_radial_arrival_velocity(self) -> float:
        ...

    def get_radial_departure_velocity(self) -> float:
        ...

    def get_semi_major_axis(self) -> float:
        ...

    def get_transverse_arrival_velocity(self) -> float:
        ...

    def get_transverse_departure_velocity(self) -> float:
        ...

class MultiRevolutionLambertTargeterIzzo(ZeroRevolutionLambertTargeterIzzo):

    def __init__(self, departure_position: numpy.ndarray, arrival_position: numpy.ndarray, time_of_flight: float, gravitational_parameter: float, n_revolutions: int=0, is_right_branch: bool=False, is_retrograde: bool=False, tolerance: float=1e-09, max_iter: int=50) -> None:
        ...

    def compute_for_revolutions_and_branch(self) -> float:
        ...

    def get_max_n_revolutions(self) -> float:
        ...

class PericenterFindingFunctions:

    def __init__(self, absolute_incoming_semi_major_axis: float, absolute_outgoing_semi_major_axis: float, bending_angle: float) -> None:
        ...

    def compute_derivative_pericenter_radius_fn(self, arg0: float) -> float:
        ...

    def compute_pericenter_radius_fn(self, arg0: float) -> float:
        ...

class ZeroRevolutionLambertTargeterIzzo(LambertTargeter):

    def __init__(self, departure_position: numpy.ndarray, arrival_position: numpy.ndarray, time_of_flight: float, gravitational_parameter: float, is_retrograde: bool=False, tolerance: float=1e-09, max_iter: int=50) -> None:
        ...

    def get_radial_arrival_velocity(self) -> float:
        ...

    def get_radial_departure_velocity(self) -> float:
        ...

    def get_semi_major_axis(self) -> float:
        ...

    def get_transverse_arrival_velocity(self) -> float:
        ...

    def get_transverse_departure_velocity(self) -> float:
        ...

def compute_escape_or_capture_delta_v(gravitational_param: float, semi_major_axis: float, eccentricity: float, excess_velocity: float) -> float:
    ...

def propagate_kepler_orbit(initial_kepler_elements: numpy.ndarray, propagation_time: float, gravitational_parameter: float, root_finder: tudatpy.math.root_finders.expose_root_finders.RootFinderCore=None) -> numpy.ndarray:
    """Function to propagate Keplerian elements to a later epoch, assuming an unperturbed system.
	
	Function to propagate Keplerian elements to a later epoch, assuming an unperturbed system. This function will
	take the initial Keplerian elements, and propagate the true anomaly in time as per the requested input. This
	is done by converting true anomaly to mean anomaly, apply the constant rate in mean motion for the requested
	time, and converting the result back to true anomaly. Currently both elliptic and hyperbolic orbits are supported.
	Parabolic orbits are not supported and will result in an error message.
	
	
	:param initial_kepler_elements:
			Keplerian elements that are to be propagated (see :ref:`\\`\\`element_conversion\\`\\`` for order)
	:param propagation_time:
			Time for which the elements are to be propagated w.r.t. the initial elements
	:param gravitational_parameter:
			Gravitational parameter of central body used for propagation
	:param root_finder:
			Root finder used to solve Kepler's equation when converting mean to eccentric anomaly. When no root finder is specified, the default option of the mean to eccentric anomaly function is used (see :func:`~mean_to_eccentric_anomaly').
	:return:
			Keplerian elements, propagated in time from initial elements assuming unperturbed dynamics. Note that the true anomaly is returned within the -PI to PI spectrum. If the user desires a different spectrum (possibly including the number of revolutions), these should be added by the user a posteriori.
	"""