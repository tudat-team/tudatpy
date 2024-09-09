import numpy
import tudatpy.numerical_simulation.environment.expose_environment
import typing
__all__ = ['AtmosphereDependentVariables', 'AtmosphereSettings', 'ExponentialAtmosphereSettings', 'NRLMSISE00Atmosphere', 'WindModelSettings', 'constant_wind_model', 'custom_constant_temperature', 'custom_four_dimensional_constant_temperature', 'custom_wind_model', 'exponential', 'exponential_predefined', 'nrlmsise00', 'scaled_by_constant', 'scaled_by_function', 'tabulated', 'tabulated_density', 'tabulated_gas_constant', 'tabulated_molar_mass', 'tabulated_pressure', 'tabulated_specific_heat_ratio', 'tabulated_temperature', 'us76']

class AtmosphereDependentVariables:
    """Members:
	
	tabulated_density
	
	tabulated_pressure
	
	tabulated_temperature
	
	tabulated_gas_constant
	
	tabulated_specific_heat_ratio
	
	tabulated_molar_mass
	"""
    __members__: typing.ClassVar[dict[str, AtmosphereDependentVariables]]
    tabulated_density: typing.ClassVar[AtmosphereDependentVariables]
    tabulated_gas_constant: typing.ClassVar[AtmosphereDependentVariables]
    tabulated_molar_mass: typing.ClassVar[AtmosphereDependentVariables]
    tabulated_pressure: typing.ClassVar[AtmosphereDependentVariables]
    tabulated_specific_heat_ratio: typing.ClassVar[AtmosphereDependentVariables]
    tabulated_temperature: typing.ClassVar[AtmosphereDependentVariables]

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

class AtmosphereSettings:
    """Base class for providing settings for atmosphere model.
	
	Functional (base) class for settings of atmosphere models that require no information in addition to their type.
	Atmosphere model classes requiring additional information must be created using an object derived from this class.
	"""

    @property
    def wind_settings(self) -> WindModelSettings:
        """
        Wind model settings for the atmosphere model settings object.
        	
        """

    @wind_settings.setter
    def wind_settings(self, arg1: WindModelSettings) -> None:
        ...

class ExponentialAtmosphereSettings(AtmosphereSettings):
    """Class for providing settings for exponential atmosphere model.
	
	`AtmosphereSettings` derived class for a defining the settings of an exponential atmosphere model.
	"""

class NRLMSISE00Atmosphere:
    """NRLMSISE00 atmosphere model.
	
	This class uses the NRLMSISE00 model to compute the atmospheric density and temperature. The GTD7 function is used: Neutral Atmosphere Empirical Model from the surface to the lower exosphere.
	
	Currently, the ideal gas law is used to compute the speed of sound and the specific heat ratio is assumed to be constant and equal to 1.4.
	
	:param solar_activity_data: Solar activity data for a range of epochs as produced by tudatpy.io.read_solar_activity_data.
	
	"""

    def __init__(self, solar_activity_data: dict[float, ...]) -> None:
        ...

    def get_density(self, altitude: float, longitude: float, latitude: float, time: float) -> float:
        """
        Get local density
        
                                    Returns the local density at the given altitude,
                                    longitude, latitude and time.
        
                                    :param altitude: Altitude at which to get the density. [m]
                                    :param longitude: Longitude at which to get the density [rad].
                                    :param latitude: Latitude at which to get the density [rad].
                                    :param time: Time at which density is to be computed [seconds since J2000].
                                    :return: Local density. [kg/m^3]
        """

class WindModelSettings:
    """Class for providing settings for wind model.
	
	Functional (base) class for settings of wind models that require no information in addition to their type.
	Wind model classes requiring additional information must be created using an object derived from this class.
	"""

def constant_wind_model(wind_velocity: numpy.ndarray, associated_reference_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicsReferenceFrames=...) -> WindModelSettings:
    """Factory function for creating wind model settings with constant wind velocity.
	
	Factory function for settings object, defining wind model entirely from constant wind velocity in a given reference frame.
	
	
	:param wind_velocity:
			Constant wind velocity in the specified reference frame.
	
	:param associated_reference_frame:
			Reference frame in which constant wind velocity is defined.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ConstantWindModelSettings` class
	"""

def custom_constant_temperature(density_function: typing.Callable[[float], float], constant_temperature: float, specific_gas_constant: float=287.0, ratio_of_specific_heats: float=1.4) -> AtmosphereSettings:
    """Factory function for creating atmospheric model settings from custom density profile.
	
	Factory function for settings object, defining constant temperature atmosphere model from custom density profile.
	The user is specifying the density profile as a function of altitude.
	The value of pressure is computed by assuming hydrostatic equilibrium, temperature, gas constant and the ratio of specific heats are modelled as constants.
	
	
	:param density_function:
			Function to retrieve the density at the current altitude.
	
	:param constant_temperature:
			Constant atmospheric temperature.
	:param specific_gas_constant:
			Specific gas constant for (constant) atmospheric chemical composition.
	:param ratio_specific_heats:
			Ratio of specific heats for (constant) atmospheric chemical composition.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class
	"""

def custom_four_dimensional_constant_temperature(density_function: typing.Callable[[float, float, float, float], float], constant_temperature: float, specific_gas_constant: float=287.0, ratio_of_specific_heats: float=1.4) -> AtmosphereSettings:
    """Factory function for creating atmospheric model settings from custom density profile.
	
	Factory function for settings object, defining constant temperature atmosphere model from custom density profile.
	The user is specifying the density profile as a function of altitude, longitude, latitude and time.
	
	.. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
			  The altitude is in meters, and the time is a Julian date in seconds since J2000.
	
	
	:param density_function:
			Function to retrieve the density at the current altitude, longitude, latitude and time.
	
	:param constant_temperature:
			Constant atmospheric temperature.
	:param specific_gas_constant:
			Specific gas constant for (constant) atmospheric chemical composition.
	:param ratio_specific_heats:
			Ratio of specific heats for (constant) atmospheric chemical composition.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class
	"""

def custom_wind_model(wind_function: typing.Callable[[float, float, float, float], numpy.ndarray], associated_reference_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicsReferenceFrames=...) -> WindModelSettings:
    """Factory function for creating wind model settings with custom wind velocity.
	
	Factory function for settings object, defining wind model entirely from custom wind velocity function in a given reference frame.
	The custom wind velocity has to be given as a function of altitude, longitude, latitude and time.
	
	.. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
			  The altitude is in meters, and the time is a Julian date in seconds since J2000.
	
	
	:param wind_velocity:
			Custom wind velocity function (w.r.t. altitude, longitude, latitude and time) in the specified reference frame.
	
	:param associated_reference_frame:
			Reference frame in which wind velocity is defined.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomWindModelSettings` class
	"""

def exponential(scale_height: float, surface_density: float, constant_temperature: float=288.15, specific_gas_constant: float=287.0, ratio_specific_heats: float=1.4) -> AtmosphereSettings:
    """Factory function for creating atmospheric model settings from fully parametrized exponential model.
	
	Factory function for settings object, defining exponential atmosphere model.
	The model is solely based on an exponentially decaying density profile with a constant temperature and composition
	(i.e. independent of time, latitude and longitude).
	
	The user has access to a fully parametrized model, meaning that in addition to the required input parameters ``scale_height`` and ``surface_density`` (ground-level air density),
	the user can specify non-standard values for constant temperature, gas constant and specific heats ratio.
	
	
	:param scale_height:
			Scale height for density profile of atmosphere.
	:param surface_density:
			Atmospheric density at ground level.
	:param constant_temperature:
			Constant atmospheric temperature.
	:param specific_gas_constant:
			Specific gas constant for (constant) atmospheric chemical composition.
	:param ratio_specific_heats:
			Ratio of specific heats for (constant) atmospheric chemical composition.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class
	"""

def exponential_predefined(body_name: str) -> AtmosphereSettings:
    """Factory function for creating atmospheric model settings from pre-defined exponential model.
	
	Factory function for settings object, defining atmosphere model from pre-defined exponential model.
	The pre-encoded properties are available for Earth and Mars, as can be seen on the table below.
	This function creates an instance of an `AtmosphereSettings` derived `ExponentialAtmosphereSettings` object.
	
	.. list-table:: Pre-defined exponential atmosphere model properties
	  :widths: 25 25 25 25
	  :header-rows: 1
	
	  * - Property
		- Earth
		- Mars
		- Units
	  * - Scale Height
		- 7.2
		- 11.1
		- km
	  * - Density at Zero Altitude
		- 1.225
		- 0.02
		- kg/m :math:`{}^3`
	  * - Constant Temperature
		- 246.0
		- 215.0
		- K
	  * - Specific Gas Constant
		- 287.0
		- 197.0
		- J/kg/K
	  * - Ratio of Specific Heats
		- 1.4
		- 1.3
		- --
	
	
	:param body_name:
			Body for which pre-defined model settings are to be loaded. Available bodies "Earth", "Mars".
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class
	"""

def nrlmsise00(space_weather_file: str='/Users/alfonso/.tudat/resource/space_weather/sw19571001.txt') -> AtmosphereSettings:
    """Factory function for creating NRLMSISE-00 atmospheric model settings.
	
	Factory function for settings object, defining atmosphere model in accordance to the NRLMSISE-00 global reference model for Earth's atmosphere.
	
	
	:param space_weather_file:
			File to be used for space weather characteristics as a function of time (e.g. F10.7, Kp, etc.). The file is typically taken from here `celestrak <https://celestrak.org/SpaceData/sw19571001.txt>`_ (note that the file in your resources path will not be the latest version of this file; download and replace your existing file if required). Documentation on the file is given `here <https://celestrak.org/SpaceData/SpaceWx-format.php>`_
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` class
	"""

def scaled_by_constant(unscaled_atmosphere_settings: AtmosphereSettings, density_scaling: float, is_scaling_absolute: bool=False) -> AtmosphereSettings:
    """Factory function for creating scaled atmospheric model settings.
	
	Factory function for settings object, defining atmospheric model based on an scaling of an existing atmospheric settings object.
	The user can apply a scaling factor (or an absolute value) to the air densities of the existing model settings (for instance for an uncertainty analysis).
	
	
	:param unscaled_atmosphere_settings:
			Sets base settings of atmosphere model to be scaled.
	:param density_scaling:
			Constant scaling factor to be applied to the entire air density profile.
	:param is_scaling_absolute:
			Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.
	"""

def scaled_by_function(unscaled_atmosphere_settings: AtmosphereSettings, density_scaling_function: typing.Callable[[float], float], is_scaling_absolute: bool=False) -> AtmosphereSettings:
    """Factory function for creating scaled atmospheric model settings.
	
	Factory function for settings object, defining atmospheric model based on scaling an existing atmospheric settings object.
	The user can apply custom scaling factors (or absolute values) to the air densities of the existing model settings (for instance for an uncertainty analysis).
	
	
	:param unscaled_atmosphere_settings:
			Sets base settings of atmosphere model to be scaled.
	:param density_scaling_function:
			Specifies air density scaling factor as a function of time.
	:param is_scaling_absolute:
			Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.
	"""

def tabulated(atmosphere_data_file: str, dependent_variable_names: list[AtmosphereDependentVariables]=..., specific_gas_constant: float=287.0, ratio_of_specific_heats: float=1.4) -> AtmosphereSettings:
    ...

def us76() -> AtmosphereSettings:
    ...
tabulated_density: AtmosphereDependentVariables
tabulated_gas_constant: AtmosphereDependentVariables
tabulated_molar_mass: AtmosphereDependentVariables
tabulated_pressure: AtmosphereDependentVariables
tabulated_specific_heat_ratio: AtmosphereDependentVariables
tabulated_temperature: AtmosphereDependentVariables