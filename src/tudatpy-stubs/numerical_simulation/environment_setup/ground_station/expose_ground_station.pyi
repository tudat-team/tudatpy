import numpy
import tudatpy.astro.element_conversion.expose_element_conversion
import typing
__all__ = ['CustomGroundStationMotionSettings', 'GroundStationMotionSettings', 'GroundStationSettings', 'LinearGroundStationMotionSettings', 'PiecewiseConstantGroundStationMotionSettings', 'basic_station', 'custom_station_motion', 'dsn_stations', 'linear_station_motion', 'piecewise_constant_station_motion']

class CustomGroundStationMotionSettings(GroundStationMotionSettings):
    """Class for defining custom time-dependent motion of a ground station.
	
	`CustomGroundStationMotionSettings` derived class for custom time-dependent motion of a ground station
	"""

class GroundStationMotionSettings:
    """Base class for providing settings for the motion of a single ground station.
	
	Non-functional base class for settings for the motion of a single ground station
	Station motion settings requiring additional information must be defined using an object derived from this class.
	"""

class GroundStationSettings:
    """Base class for providing settings for the creation of a ground station.
	"""

class LinearGroundStationMotionSettings(GroundStationMotionSettings):
    """Class for defining linear motion (in an Earth-fixed frame) in time of a ground station.
	
	`GroundStationMotionSettings` derived class for time-linear station motion
	"""

class PiecewiseConstantGroundStationMotionSettings(GroundStationMotionSettings):
    """Class for defining piecewise-constant position (e.g. instantaneous change in position at given epochs) of a ground station.
	
	`GroundStationMotionSettings` derived class for piecewise-constant position of a ground station
	"""

def basic_station(station_name: str, station_nominal_position: numpy.ndarray, station_position_element_type: tudatpy.astro.element_conversion.expose_element_conversion.PositionElementTypes=..., station_motion_settings: list[GroundStationMotionSettings]=[]) -> GroundStationSettings:
    """Factory function for creating settings for a ground station
	
	Factory function for creating settings for a ground station, defining only its name, body-fixed position, and (optionally) time-variations of its position
	
	
	:param station_name:
			Name (unique identifier) by which the station is to be known.
	:param station_nominal_position:
			Nominal position of the station in a body-fixed frame. Depending on the choice of ``station_position_element_type`` input, this vector must contain
			* Cartesian - :math:`[x,y,z]`, denoting :math:`x-`, :math:`y-` and :math:`z-` components of body-fixed position (w.r.t body-fixed frame origin, typically center of mass) * Spherical - :math:`[r,\\phi',		heta]`, denoting distance from body-fixed frame origin (typically center of mass), latitude and longitude * Geodetic - :math:`[h,\\phi,  heta]`, denoting the altitude w.r.t. the body shape model, geodetic latitude and longitude
	:param station_position_element_type:
			Type of elements for ``station_nominal_position``
	:param station_motion_settings:
			List of settings defining time-variations of the individual ground station
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.GroundStationSettings` defining settings of the to be created ground station
	"""

def custom_station_motion(custom_displacement_function: typing.Callable[[float], numpy.ndarray]) -> GroundStationMotionSettings:
    """Factory function for creating settings for a custom ground station position variation
	
	Factory function for creating settings for a custom ground station position. An arbitrary user-defined function of the signature :math:`\\Delta\\mathbf{r}=\\Delta\\mathbf{r}(t)` is provided and
	applied to the station position
	
	
	:param custom_displacement_function:
			Function returning :math:`\\Delta\\mathbf{r}`, with the time :math:`t` as input.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.CustomGroundStationMotionSettings` class
	"""

def dsn_stations() -> list[GroundStationSettings]:
    """Factory function for creating settings for all DSN stations
	
	Factory function for creating settings for all DSN stations, defined by nominal positions and linear velocities, as defined
	by Cartesian elements in *DSN No. 810-005, 301, Rev. K*,  see `this link <https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf>`_.
	Note that calling these settings will use the Cartesian elements provided in this document (in ITRF93) and apply them to the Earth-fixed
	station positions, regardless of the selected Earth rotation model.
	
	:return:
			List of settings to create DSN stationss
	"""

def linear_station_motion(linear_velocity: numpy.ndarray, reference_epoch: float=0.0) -> GroundStationMotionSettings:
    """Factory function for creating settings for a linear station motion
	
	Factory function for creating settings for a linear station motion, implementing :math:`\\Delta \\mathbf{r}=\\dot{\\mathbf{r}}(t-t_{0})`.
	
	
	:param linear_velocity:
			Linear velocity :math:`\\dot{\\mathbf{r}}` of the station (in m/s)
	:param reference_epoch:
			Reference epoch :math:`t_{0}`, in seconds since J2000 epoch
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.LinearGroundStationMotionSettings` class
	"""

def piecewise_constant_station_motion(displacement_list: dict[float, numpy.ndarray]) -> GroundStationMotionSettings:
    """Factory function for creating settings for a piecewise constant ground station position variation
	
	Factory function for creating settings for a piecewise constant ground station position. Using this model, the added station velocity in a body-fixed frame :math:`\\dot{\\mathbf{r}}` is
	always zero, but its displacement :math:`\\Delta\\mathbf{r}` is set according to the input list, which contains a list of times and displacments :math:`[t_{i},\\Delta\\mathbf{r}_{i}]`.
	When the resulting model is queried at a given time :math:`t`, the nearest lower neighbour :math:`t_{i}` from this list is found, and the associated :math:`\\Delta\\mathbf{r}_{i}` is applied.
	
	
	:param displacement_list:
			Dictionary with the epochs :math:`t_{i}` as values, and the associated displacement :math:`\\Delta\\mathbf{r}_{i}` as value
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings` class
	"""