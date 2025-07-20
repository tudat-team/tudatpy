import numpy
import pybind11_stubgen.typing_ext
from ...math import root_finders
from ...trajectory_design import shape_based_thrust
import typing
__all__ = ['CaptureAndInsertionNodeSettings', 'DEFAULT_MINIMUM_PERICENTERS', 'EscapeAndDepartureNodeSettings', 'HodographicShapingLeg', 'SphericalShapingLeg', 'SwingbyNodeSettings', 'TransferLeg', 'TransferLegSettings', 'TransferLegTypes', 'TransferNodeSettings', 'TransferTrajectory', 'capture_node', 'create_transfer_trajectory', 'departure_node', 'dsm_position_based_leg', 'dsm_position_based_leg_type', 'dsm_velocity_based_leg', 'dsm_velocity_based_leg_type', 'hodographic_low_thrust_leg', 'hodographic_shaping_leg', 'mga_settings_dsm_position_based_legs', 'mga_settings_dsm_velocity_based_legs', 'mga_settings_hodographic_shaping_legs', 'mga_settings_hodographic_shaping_legs_with_recommended_functions', 'mga_settings_spherical_shaping_legs', 'mga_settings_unpowered_unperturbed_legs', 'print_parameter_definitions', 'set_low_thrust_acceleration', 'spherical_shaping_leg', 'spherical_shaping_low_thrust_leg', 'swingby_node', 'unpowered_leg', 'unpowered_unperturbed_leg_type']

class CaptureAndInsertionNodeSettings(TransferNodeSettings):
    """Class for defining settings of capture and insertion node.
    
    `TransferNodeSettings` derived class for providing settings for capture and insertion nodes, which consist of the
    capture semi-major axis and eccentricity."""

class EscapeAndDepartureNodeSettings(TransferNodeSettings):
    """Class for defining settings of escape and departure node.
    
    `TransferNodeSettings` derived class for providing settings for escape and departure nodes, which consist of the
    departure semi-major axis and eccentricity."""

class HodographicShapingLeg(TransferLeg):
    """Class for defining low-thrust hodographic-shaping leg.
    
    `TransferLeg` derived class for defining a low-thrust leg described using hodographic shaping :cite:p:`gondelach2012`."""

class SphericalShapingLeg(TransferLeg):
    """Class for defining low-thrust spherical-shaping leg.
    
    `TransferLeg` derived class for defining a low-thrust leg described using spherical shaping :cite:p:`roegiers2014`."""

class SwingbyNodeSettings(TransferNodeSettings):
    """Class for defining settings of swingby node.
    
    `TransferNodeSettings` derived class for providing settings for swingby nodes, which consist of the minimum periapsis
    radius."""

class TransferLeg:
    """Base class for defining a transfer leg.
    
    Functional (base) class for transfer legs, requiring the leg type, departure body ephemeris and arrival body ephemeris.
    Transfer node classes requiring additional information must be created using an object derived from this class."""

    def state_along_trajectory(self, time_since_leg_beginning: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
        """
        No documentation found.
        """

class TransferLegSettings:
    """Base class for providing settings for transfer legs.
    
    Functional (base) class for settings of transfer legs that require no information in addition to their type."""

class TransferLegTypes:
    """Enumeration of available leg types.
    
    
    
    
    
          
    
    Members:
    
      unpowered_unperturbed_leg_type : 
          
    
      dsm_position_based_leg_type : 
          
    
      dsm_velocity_based_leg_type : 
          
    
      spherical_shaping_low_thrust_leg : 
          
    
      hodographic_low_thrust_leg : 
          """
    __members__: typing.ClassVar[dict[str, TransferLegTypes]]
    dsm_position_based_leg_type: typing.ClassVar[TransferLegTypes]
    dsm_velocity_based_leg_type: typing.ClassVar[TransferLegTypes]
    hodographic_low_thrust_leg: typing.ClassVar[TransferLegTypes]
    spherical_shaping_low_thrust_leg: typing.ClassVar[TransferLegTypes]
    unpowered_unperturbed_leg_type: typing.ClassVar[TransferLegTypes]

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

class TransferNodeSettings:
    """Base class for providing settings for transfer nodes.
    
    Functional (base) class for settings of transfer nodes that require no information in addition to their type.
    Transfer node classes requiring additional information must be created using an object derived from this class."""

class TransferTrajectory:
    """Class defining a transfer trajectory constituted by transfer legs and nodes.
    
    Class defining a transfer trajectory constituted by transfer legs and nodes. The object is tipically created using the `create_transfer_trajectory` function."""

    def evaluate(self, node_times: list[float], leg_parameters: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], node_parameters: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
                 Evaluate transfer trajectory.
        
                 Function to evaluate the transfer trajectory, which consists of computing the transfer Delta V and time of
                 flight, for the specified set of parameters.
        
        
                 Parameters
                 ----------
                 node_times : list[float]
                     List of the time at each node.
                 leg_parameters : list[list[float]]
                     List of lists with the parameters characterizing each leg. Each inner list corresponds to the
                     parameters of one leg; if a leg does not require any parameter, its list can contain any value(s),
                     therefore it is recommended to leave it empty.
        
                 node_parameters : list[list[float]]
                     List of lists with the parameters characterizing each node. Each inner list corresponds to the
                     parameters of one node; if a node does not require any parameter, its list can contain any value(s),
                     therefore it is recommended to leave it empty.
        
                 Returns
                 -------
                 None
                     None
        """

    def inertial_thrust_accelerations_along_trajectory(self, number_of_data_points_per_leg: int) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]:
        """
                 Returns the inertial thrust acceleration history throughout the trajectory.
        
                 Function that returns the inertial thrust acceleration history throughout the trajectory, using the same number of data points in
                 each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
                 For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.
        
        
                 Parameters
                 ----------
                 number_of_data_points_per_leg : int
                     Number of data points used to describe each leg.
                 Returns
                 -------
                 tuple[numpy.ndarray,numpy.ndarray]
                     Tuple of (state history, time history).
        """

    def rsw_thrust_accelerations_along_trajectory(self, number_of_data_points_per_leg: int) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]:
        """
                 Returns the thrust acceleration history in the RSW frame throughout the trajectory.
        
                 Function that returns the thrust acceleration history in the RSW frame throughout the trajectory, using the same number of data points in
                 each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
                 For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.
        
        
                 Parameters
                 ----------
                 number_of_data_points_per_leg : int
                     Number of data points used to describe each leg.
                 Returns
                 -------
                 tuple[numpy.ndarray,numpy.ndarray]
                     Tuple of (state history, time history).
        """

    def single_leg_delta_v(self, leg_index: int) -> float:
        """
                 Retrieves the Delta V applied in the specified leg.
        
        
                 Parameters
                 ----------
                 leg_index : int
                     Index of the leg for which the Delta V is to be retrieved.
                 Returns
                 -------
                 float
                     Delta V for the specified leg.
        """

    def single_node_delta_v(self, node_index: int) -> float:
        """
                 Retrieves the Delta V applied in the specified node.
        
        
                 Parameters
                 ----------
                 node_index : int
                     Index of the node for which the Delta V to be is retrieved.
                 Returns
                 -------
                 float
                     Delta V for the specified node.
        """

    def states_along_trajectory(self, number_of_data_points_per_leg: int) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]]:
        """
                 Returns the state history throughout the trajectory.
        
                 Function that returns the state history throughout the trajectory, using the same number of data points in
                 each leg. For each leg, the retrieved states are equally spaced in time.
        
        
                 Parameters
                 ----------
                 number_of_data_points_per_leg : int
                     Number of data points used to describe each leg.
                 Returns
                 -------
                 tuple[numpy.ndarray,numpy.ndarray]
                     Tuple of (state history, time history).
        """

    def tnw_thrust_accelerations_along_trajectory(self, number_of_data_points_per_leg: int) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]:
        """
                 Returns the thrust acceleration history in the TNW frame throughout the trajectory.
        
                 Function that returns the thrust acceleration history in the TNW frame throughout the trajectory, using the same number of data points in
                 each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
                 For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.
        
        
                 Parameters
                 ----------
                 number_of_data_points_per_leg : int
                     Number of data points used to describe each leg.
                 Returns
                 -------
                 tuple[numpy.ndarray,numpy.ndarray]
                     Tuple of (state history, time history).
        """

    @property
    def delta_v(self) -> float:
        """
                 **read-only**
        
                 Total Delta V used in the transfer trajectory.
        
                 :type: float
        """

    @property
    def delta_v_per_leg(self) -> list[float]:
        """
                 **read-only**
        
                 List of the Delta V applied in each leg.
        
                 :type: list[float]
        """

    @property
    def delta_v_per_node(self) -> list[float]:
        """
                 **read-only**
        
                 List of the Delta V applied in each node.
        
                 :type: list[float]
        """

    @property
    def legs(self) -> list[TransferLeg]:
        """
        No documentation found.
        """

    @property
    def number_of_legs(self) -> int:
        """
                 **read-only**
        
                 Number of legs in the transfer trajectory.
        
                 :type: float
        """

    @property
    def number_of_nodes(self) -> int:
        """
                 **read-only**
        
                 Number of nodes in the transfer trajectory.
        
                 :type: float
        """

    @property
    def time_of_flight(self) -> float:
        """
                 **read-only**
        
                 Total time of flight of the transfer trajectory.
        
                 :type: float
        """

def capture_node(capture_semi_major_axis: float, capture_eccentricity: float) -> TransferNodeSettings:
    """Function for creating the settings of a capture or insertion node.
    
    Function for creating the settings of a capture or insertion node. The settings consist of the
    capture orbit eccentricity and semi-major axis.
    Given the the arrival velocity and the final orbit, the node computes the Delta V that needs to be applied
    at the periapsis of the final orbit to exit the capture trajectory.
    The calculations performed in this node do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`musegaas2012`.
    
    
    Parameters
    ----------
    capture_semi_major_axis : float
        Capture orbit semi-major axis.
    capture_eccentricity : float
        Capture orbit eccentricity.
    Returns
    -------
    CaptureAndInsertionNodeSettings
        Capture or insertion node settings."""

def create_transfer_trajectory(bodies: ..., leg_settings: list[TransferLegSettings], node_settings: list[TransferNodeSettings], node_names: list[str], central_body: str) -> TransferTrajectory:
    """Function for creating a transfer trajectory consisting of the specified sequence of transfer nodes and
    transfer legs.
    
    
    The function creates a transfer trajectory based on the provided transfer nodes settings and transfer legs
    settings. The number of nodes should be equal to the number of legs plus 1.
    This function creates an instance of the `TransferTrajectory` class.
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        System of bodies to be used in the transfer trajectory.
    leg_settings : list[TransferLegSettings]
        List of transfer leg settings.
    node_settings : list [TransferNodeSettings]
        List of transfer node settings.
    node_names : list [str]
        Sequence of bodies used as transfer nodes.
    central_body : str
        Central body with respect to which the two-body trajectory of the spacecraft
        is calculated.
    
    Returns
    -------
    TransferTrajectory
        Transfer trajectory object."""

def departure_node(*args, **kwargs) -> TransferNodeSettings:
    """Function for creating the settings of an escape or departure node.
    
    Function for creating the settings of an escape or departure node. The settings consist of the
    departure orbit eccentricity and semi-major axis.
    Given the initial orbit and the departure velocity, the node computes the Delta V that needs to be applied
    at the periapsis of the initial orbit to enter the escape trajectory.
    The calculations performed in this node do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`musegaas2012`.
    
    
    Parameters
    ----------
    departure_semi_major_axis : float
        Departure orbit semi-major axis.
    departure_eccentricity : float
        Departure orbit eccentricity.
    Returns
    -------
    EscapeAndDepartureNodeSettings
        Escape or departure node settings."""

def dsm_position_based_leg() -> TransferLegSettings:
    """Function for creating the settings of a transfer leg with 1 impulsive deep space maneuver (DSM) described using
    the position formulation.
    
    
    Function for creating the settings of a transfer leg with 1 position-based DSM; the settings consist of just the leg type.
    Given the departure position and the DSM position this leg uses a Lambert targeter to compute the departure
    velocity and the velocity before the DSM. Given the DSM position and the arrival position, the leg also uses
    a Lambert targeter to compute the velocity after the DSM and the arrival velocity. The Delta V applied in the
    DSM is computed using the velocity before and after the DSM.
    The calculations performed in this leg do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`musegaas2012`.
    
    Returns
    -------
    TransferLegSettings
        Transfer leg settings."""

def dsm_velocity_based_leg() -> TransferLegSettings:
    """Function for creating the settings of a transfer leg with 1 impulsive deep space maneuver (DSM) described using
    the velocity formulation.
    
    
    Function for creating the settings of a transfer leg with 1 velocity-based DSM; the settings consist of just the leg type.
    Given the departure position and velocity this leg "propagates" the Kepler elements until the instant of application
    of the DSM (giving the position at the DSM and the velocity before the DSM). Given the position of the DSM
    and the arrival position, it computes the velocity after the DSM
    (which is used to compute the Delta V applied in the DSM) and the arrival velocity using a Lambert targeter.
    The calculations performed in this leg do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`musegaas2012`.
    
    Returns
    -------
    TransferLegSettings
        Transfer leg settings."""

def hodographic_shaping_leg(radial_velocity_function_components: list[shape_based_thrust.BaseFunctionHodographicShaping], normal_velocity_function_components: list[shape_based_thrust.BaseFunctionHodographicShaping], axial_velocity_function_components: list[shape_based_thrust.BaseFunctionHodographicShaping]) -> TransferLegSettings:
    """Function for creating the settings of a low-thrust hodographic shaping leg.
    
    
    Function for creating the settings of a low-thrust hodographic shaping leg; the settings consist of
    the functions used to shape the velocity throughout the transfer.
    Note that shape functions with at least 3 terms (i.e. 3 degrees of freedom) must be provided for each velocity
    component (radial, normal and axial); this is required to obtain a trajectory that satisfies the boundary
    conditions.
    The calculations performed in this leg do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`gondelach2012`.
    
    
    Parameters
    ----------
    radial_velocity_function_components : list[BaseFunctionHodographicShaping]
        List with components of the radial velocity shaping function, which determine the radial velocity profile
        throughout the transfer.
    
    normal_velocity_function_components : list[BaseFunctionHodographicShaping]
        List with components of the normal velocity shaping function, which determine the normal velocity profile
        throughout the transfer.
    
    axial_velocity_function_components : list[BaseFunctionHodographicShaping]
        List with components of the axial velocity shaping function, which determine the axial velocity profile
        throughout the transfer.
    
    Returns
    -------
    TransferLegSettings
        Transfer leg settings."""

def mga_settings_dsm_position_based_legs(body_order: list[str], departure_orbit: tuple[float, float]=..., arrival_orbit: tuple[float, float]=..., minimum_pericenters: dict[str, float]={'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}) -> tuple[list[TransferLegSettings], list[TransferNodeSettings]]:
    """Function to get the legs and nodes settings of a transfer constituted by legs with 1 impulsive deep space maneuver (DSM)
    described using the position formulation.
    
    
    Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
    one initial node (departure or swingby), position-based DSM transfer leg(s) connected by swingby nodes, one final
    (capture or swingby) node.
    If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
    final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
    node.
    
    
    Parameters
    ----------
    body_order : list[str]
        List of bodies to visit, including departure body, swingby bodies and arrival body.
    departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
        node as a swingby node (instead of a departure node).
    
    arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
        node as a swingby node (instead of a capture node).
    
    minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
        Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
        value. Default values from :cite:t:`izzo2007`.
    
    Returns
    -------
    tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
        Tuple specifying the settings of each transfer leg and node."""

def mga_settings_dsm_velocity_based_legs(body_order: list[str], departure_orbit: tuple[float, float]=..., arrival_orbit: tuple[float, float]=..., minimum_pericenters: dict[str, float]={'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}) -> tuple[list[TransferLegSettings], list[TransferNodeSettings]]:
    """Function to get the legs and nodes settings of a transfer constituted by legs with 1 impulsive deep space maneuver (DSM)
    described using the velocity formulation.
    
    
    Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
    one initial node (departure or swingby), velocity-based DSM transfer leg(s) connected by swingby nodes, one final
    (capture or swingby) node.
    If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
    final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
    node.
    
    
    Parameters
    ----------
    body_order : list[str]
        List of bodies to visit, including departure body, swingby bodies and arrival body.
    departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
        node as a swingby node (instead of a departure node).
    
    arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
        node as a swingby node (instead of a capture node).
    
    minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
        Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
        value. Default values from :cite:t:`izzo2007`.
    
    Returns
    -------
    tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
        Tuple specifying the settings of each transfer leg and node."""

def mga_settings_hodographic_shaping_legs(body_order: list[str], radial_velocity_function_components_per_leg: list[list[shape_based_thrust.BaseFunctionHodographicShaping]], normal_velocity_function_components_per_leg: list[list[shape_based_thrust.BaseFunctionHodographicShaping]], axial_velocity_function_components_per_leg: list[list[shape_based_thrust.BaseFunctionHodographicShaping]], departure_orbit: tuple[float, float]=..., arrival_orbit: tuple[float, float]=..., minimum_pericenters: dict[str, float]={'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}) -> tuple[list[TransferLegSettings], list[TransferNodeSettings]]:
    """Function to get the legs and nodes settings of a transfer constituted by low-thrust hodographic shaping legs,
    with user-provided velocity shaping functions.
    
    
    Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
    one initial node (departure or swingby), hodographic shaping leg(s) connected by swingby nodes, one final
    (capture or swingby) node.
    If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
    final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
    node.
    
    
    Parameters
    ----------
    body_order : list[str]
        List of bodies to visit, including departure body, swingby bodies and arrival body.
    radial_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
        List with the lists of radial velocity function components used in each leg.
    
    normal_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
        List with the lists of normal velocity function components used in each leg.
    
    axial_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
        List with the lists of axial velocity function components used in each leg.
    
    departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
        node as a swingby node (instead of a departure node).
    
    arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
        node as a swingby node (instead of a capture node).
    
    minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
        Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
        value. Default values from :cite:t:`izzo2007`.
    
    Returns
    -------
    tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
        Tuple specifying the settings of each transfer leg and node."""

def mga_settings_hodographic_shaping_legs_with_recommended_functions(body_order: list[str], time_of_flight_per_leg: list[float], number_of_revolutions_per_leg: list[float], departure_orbit: tuple[float, float]=..., arrival_orbit: tuple[float, float]=..., minimum_pericenters: dict[str, float]={'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}) -> tuple[list[TransferLegSettings], list[TransferNodeSettings]]:
    """Function to get the legs and nodes settings of a transfer constituted by low-thrust hodographic shaping legs,
    with user-provided velocity shaping functions.
    
    
    Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
    one initial node (departure or swingby), hodographic shaping leg(s) connected by swingby nodes, one final
    (capture or swingby) node.
    If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
    final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
    node.
    
    
    Parameters
    ----------
    body_order : list[str]
        List of bodies to visit, including departure body, swingby bodies and arrival body.
    radial_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
        List with the lists of radial velocity function components used in each leg.
    
    normal_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
        List with the lists of normal velocity function components used in each leg.
    
    axial_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
        List with the lists of axial velocity function components used in each leg.
    
    departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
        node as a swingby node (instead of a departure node).
    
    arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
        node as a swingby node (instead of a capture node).
    
    minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
        Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
        value. Default values from :cite:t:`izzo2007`.
    
    Returns
    -------
    tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
        Tuple specifying the settings of each transfer leg and node."""

def mga_settings_spherical_shaping_legs(body_order: list[str], root_finder_settings: root_finders.RootFinderSettings, departure_orbit: tuple[float, float]=..., arrival_orbit: tuple[float, float]=..., lower_bound_free_coefficient: float=..., upper_bound_free_coefficient: float=..., initial_value_free_coefficient: float=..., minimum_pericenters: dict[str, float]={'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}) -> tuple[list[TransferLegSettings], list[TransferNodeSettings]]:
    """Function to get the legs and nodes settings of a transfer constituted by low-thrust spherical shaping legs.
    
    
    Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
    one initial node (departure or swingby), spherical shaping leg(s) connected by swingby nodes, one final
    (capture or swingby) node.
    If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
    final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
    node.
    
    
    Parameters
    ----------
    body_order : list[str]
        List of bodies to visit, including departure body, swingby bodies and arrival body.
    root_finder_settings : RootFinderSettings
        Settings of the root finder used by the spherical shaping leg when computing the value of the free coefficient
        that allows meeting the desired time of flight.
    
    departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
        node as a swingby node (instead of a departure node).
    
    arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
        node as a swingby node (instead of a capture node).
    
    lower_bound_free_coefficient : float, default=TUDAT_NAN
        Lower bound of the possible values for the free coefficient. Parameter is potentially used by the root finder:
        it must be specified if the selected root finder requires the definition of a lower bound.
    
    upper_bound_free_coefficient : float, default=TUDAT_NAN
        Upper bound of the possible values for the free coefficient. Parameter is potentially used by the root finder:
        it must be specified if the selected root finder requires the definition of an upper bound.
    
    initial_value_free_coefficient : float, default=TUDAT_NAN
        Initial guess for the free coefficient. Parameter is potentially used by the root finder:
        it must be specified if the selected root finder requires the definition of an initial guess.
    
    minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
        Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
        value. Default values from :cite:t:`izzo2007`.
    
    Returns
    -------
    tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
        Tuple specifying the settings of each transfer leg and node."""

def mga_settings_unpowered_unperturbed_legs(body_order: list[str], departure_orbit: tuple[float, float]=..., arrival_orbit: tuple[float, float]=..., minimum_pericenters: dict[str, float]={'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}) -> tuple[list[TransferLegSettings], list[TransferNodeSettings]]:
    """Function to get the legs and nodes settings of a transfer with just upowered legs.
    
    
    Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
    one initial node (departure or swingby), unpowered transfer leg(s) connected by swingby nodes, one final
    (capture or swingby) node.
    If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
    final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
    node.
    
    
    Parameters
    ----------
    body_order : list[str]
        List of bodies to visit, including departure body, swingby bodies and arrival body.
    departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
        node as a swingby node (instead of a departure node).
    
    arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
        Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
        node as a swingby node (instead of a capture node).
    
    minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
        Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
        value. Default values from :cite:t:`izzo2007`.
    
    Returns
    -------
    tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
        Tuple specifying the settings of each transfer leg and node."""

def print_parameter_definitions(leg_settings: list[TransferLegSettings], node_settings: list[TransferNodeSettings]) -> None:
    """Prints the list of parameters required to define the transfer trajectory, according to the
    specified node and leg settings.
    
    
    
    Parameters
    ----------
    leg_settings : list[TransferLegSettings]
        List of transfer leg settings.
    node_settings : list [TransferNodeSettings]
        List of transfer node settings.
    Returns
    -------
    None
        None"""

def set_low_thrust_acceleration(transfer_leg: TransferLeg, bodies: ..., body_name: str, engine_name: str) -> None:
    """No documentation found."""

def spherical_shaping_leg(root_finder_settings: root_finders.RootFinderSettings, lower_bound_free_coefficient: float=..., upper_bound_free_coefficient: float=..., initial_value_free_coefficient: float=..., time_to_azimuth_interpolator_step_size: float=86400.0) -> TransferLegSettings:
    """Function for creating the settings of a low-thrust spherical shaping leg.
    
    
    Function for creating the settings of a low-thrust spherical shaping leg; the settings consist of
    variables necessary for setting up the root finder and variable to set up an interpolator.
    The trajectory is determined via spherical shaping, which shapes the position and time history throughout the
    transfer. The trajectory depends on a single parameter, which is selected using the root finder in order to meet
    a user-specified time of flight.
    The calculations performed in this leg do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`roegiers2014`.
    
    
    Parameters
    ----------
    root_finder_settings : RootFinderSettings
        Settings of the root finder used by the spherical shaping leg when computing the value of the free coefficient
        that allows meeting the desired time of flight.
    
    lower_bound_free_coefficient : float, default=TUDAT_NAN
        Lower bound of the possible values for the free coefficient. Parameter is potentially used by the root finder:
        it must be specified if the selected root finder requires the definition of a lower bound.
    
    upper_bound_free_coefficient : float, default=TUDAT_NAN
        Upper bound of the possible values for the free coefficient. Parameter is potentially used by the root finder:
        it must be specified if the selected root finder requires the definition of an upper bound.
    
    initial_value_free_coefficient : float, default=TUDAT_NAN
        Initial guess for the free coefficient. Parameter is potentially used by the root finder:
        it must be specified if the selected root finder requires the definition of an initial guess.
    
    time_to_azimuth_interpolator_step_size : float, default=physical_constants::JULIAN_DAY
        Time step size used as reference to define the azimuth values at which the epoch is computed, when defining an
        interpolator to convert between epoch and azimuth.
    
    Returns
    -------
    TransferLegSettings
        Transfer leg settings."""

def swingby_node(minimum_periapsis: float=...) -> TransferNodeSettings:
    """Function for creating the settings of a swingby node.
    
    Function for creating the settings of a swingby node. The settings consist consist of the minimum
    allowed periapsis radius.
    The minimum periapsis radius can be set to infinity. In that case, the swingby does not affect the velocity of the
    spacecraft (e.g. might be relevant for swingbys of small bodies).
    The exact behavior of this node depends on the types of legs that precede and follow it. Given a known incoming
    and unknown outgoing velocity, the node forward propagates the gravity assist, possibly with a Delta V at
    the a periapsis. Given an unknown incoming and known outgoing velocity, the node backward propagates the
    gravity assist, possibly with a Delta V at the a periapsis. Given known incoming and outgoing velocities,
    the node computes the Delta V required to meet those.
    The calculations performed in this node do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`musegaas2012`.
    
    
    Parameters
    ----------
    minimum_periapsis : float, default=TUDAT_NAN
        Minimum periapsis radius. The minimum periapsis only needs to be specified if the types of swignby nodes that
        requires it is used. If that is the case and no minimum periapsis was selected an error is thrown.
    
    Returns
    -------
    SwingbyNodeSettings
        Swingby node settings."""

def unpowered_leg() -> TransferLegSettings:
    """Function for creating the settings of an unpowered leg.
    
    
    Function for creating the settings of an unpowered leg; the settings consist of just the leg type.
    Given the departure and arrival position, this leg computes the departure and arrival velocity using a Lambert
    targeter.
    The calculations performed in this leg do not involve any numerical integration, and are solved by
    (semi-)analytical models. For details see :cite:t:`musegaas2012`.
    
    Returns
    -------
    TransferLegSettings
        Transfer leg settings."""
DEFAULT_MINIMUM_PERICENTERS: dict = {'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0}
dsm_position_based_leg_type: TransferLegTypes
dsm_velocity_based_leg_type: TransferLegTypes
hodographic_low_thrust_leg: TransferLegTypes
spherical_shaping_low_thrust_leg: TransferLegTypes
unpowered_unperturbed_leg_type: TransferLegTypes