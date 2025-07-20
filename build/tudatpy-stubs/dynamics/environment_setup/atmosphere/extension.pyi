import numpy
import pybind11_stubgen.typing_ext
from .... import data
from ....dynamics.environment_setup import aerodynamic_coefficients
import typing
__all__ = ['AtmosphereDependentVariables', 'AtmosphereSettings', 'ConstantWindModelSettings', 'CustomConstantTemperatureAtmosphereSettings', 'CustomWindModelSettings', 'ExponentialAtmosphereSettings', 'NRLMSISE00Atmosphere', 'NRLMSISE00Input', 'ScaledAtmosphereSettings', 'WindModelSettings', 'constant_wind_model', 'custom_constant_temperature', 'custom_four_dimensional_constant_temperature', 'custom_wind_model', 'exponential', 'exponential_predefined', 'nrlmsise00', 'scaled_by_constant', 'scaled_by_function', 'tabulated', 'tabulated_density', 'tabulated_gas_constant', 'tabulated_molar_mass', 'tabulated_pressure', 'tabulated_specific_heat_ratio', 'tabulated_temperature', 'us76']

class AtmosphereDependentVariables:
    """Members:
    
      tabulated_density
    
      tabulated_pressure
    
      tabulated_temperature
    
      tabulated_gas_constant
    
      tabulated_specific_heat_ratio
    
      tabulated_molar_mass"""
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
    Atmosphere model classes requiring additional information must be created using an object derived from this class."""

    @property
    def wind_settings(self) -> WindModelSettings:
        """
                 **read-only**
        
                 Wind model settings for the atmosphere model settings object.
        
                 :type: WindModelSettings
        """

    @wind_settings.setter
    def wind_settings(self, arg1: WindModelSettings) -> None:
        ...

class ConstantWindModelSettings(WindModelSettings):
    """No documentation found."""

class CustomConstantTemperatureAtmosphereSettings(AtmosphereSettings):
    """No documentation found."""

class CustomWindModelSettings(WindModelSettings):
    """No documentation found."""

class ExponentialAtmosphereSettings(AtmosphereSettings):
    """Class for providing settings for exponential atmosphere model.
    
    `AtmosphereSettings` derived class for a defining the settings of an exponential atmosphere model."""

class NRLMSISE00Atmosphere:
    """NRLMSISE00 atmosphere model.
    
                             This class uses the NRLMSISE00 model to compute the atmospheric density and temperature. The GTD7 function is used: Neutral Atmosphere Empirical Model from the surface to the lower exosphere.
    
                             Currently, the ideal gas law is used to compute the speed of sound and the specific heat ratio is assumed to be constant and equal to 1.4.
    
                             :param solar_activity_data: Solar activity data for a range of epochs as produced by tudatpy.io.read_solar_activity_data.
                             """

    def __init__(self, solar_activity_data: dict[float, data.SolarActivityData], use_ideal_gas_law: bool=True, use_storm_conditions: bool=False, use_anomalous_oxygen: bool=True) -> None:
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

    def get_use_geodetic_latitude(self) -> bool:
        ...

    def get_use_utc(self) -> bool:
        ...

    def set_use_geodetic_latitude(self, arg0: bool) -> None:
        ...

    def set_use_utc(self, arg0: bool) -> None:
        ...

class NRLMSISE00Input:
    """Input for computation of NRLMSISE00 atmospheric
                              conditions at current time and position.
    
                              Input for computation of NRLMSISE00 atmospheric
                              conditions at current time and position. The
                              computation of class may be reperformed every time
                              step, to reflect the changes in atmospheric
                              condition.
    
                              :param year: Current year
                              :param day_of_year: Day in the current year
                              :param seconds_of_day: Number of seconds into the
                              current day. :param local_solar_time: Local solar
                              time at the computation position :param f107: Current
                              daily F10.7 flux for previous day :param f107a: 81
                              day average of F10.7 flux (centered on current
                              day_of_year). :param ap_daily: Current daily magnetic
                              index :param ap_vector: Current magnetic index data
                              vector: \\sa ap_array :param switches: List of
                              NRLMSISE-specific flags: \\sa nrlmsise_flags """

    def __init__(self, year: int=0, day_of_year: int=0, seconds_of_day: float=0.0, local_solar_time: float=0.0, f107: float=0.0, f107a: float=0.0, ap_daily: float=0.0, ap_vector: list[float]=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], switches: list[int]=[]) -> None:
        ...

class ScaledAtmosphereSettings(AtmosphereSettings):
    """No documentation found."""

class WindModelSettings:
    """Class for providing settings for wind model.
    
    Functional (base) class for settings of wind models that require no information in addition to their type.
    Wind model classes requiring additional information must be created using an object derived from this class."""

def constant_wind_model(wind_velocity: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], associated_reference_frame: aerodynamic_coefficients.AerodynamicsReferenceFrames=...) -> WindModelSettings:
    """Function for creating wind model settings with constant wind velocity.
    
    Function for settings object, defining wind model entirely from constant wind velocity in a given reference frame.
    
    
    Parameters
    ----------
    wind_velocity : numpy.ndarray[numpy.float64[3, 1]]
        Constant wind velocity in the specified reference frame.
    
    associated_reference_frame : dynamics.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
        Reference frame in which constant wind velocity is defined.
    
    Returns
    -------
    ConstantWindModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ConstantWindModelSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`,
    using a constant wind-velocity vector defined in a vertical aerodynamic reference frame:
    
    .. code-block:: python
    
      # Define the wind in 3 directions in the vertical reference frame
      wind_Xv = 3     # Meridional wind of +3 m/s (pointing to the North)
      wind_Yv = 5     # Zonal wind of +5 m/s (pointing to the West)
      wind_Zv = -11   # Vertical wind of +11 m/s (pointing out of the centre of the Earth)
      # Create the constant wind settings
      constant_wind = environment_setup.atmosphere.constant_wind_model(
        [wind_Xv, wind_Yv, wind_Zv],
        environment.AerodynamicsReferenceFrames.vertical_frame)
      # Apply the constant wind settings to the Earth atmosphere settings
      body_settings.get("Earth").atmosphere_settings.wind_settings = constant_wind"""

def custom_constant_temperature(density_function: typing.Callable[[float], float], constant_temperature: float, specific_gas_constant: float=287.0, ratio_of_specific_heats: float=1.4) -> AtmosphereSettings:
    """Function for creating atmospheric model settings from custom density profile.
    
    Function for settings object, defining constant temperature atmosphere model from custom density profile.
    The user is specifying the density profile as a function of altitude.
    The value of pressure is computed by assuming hydrostatic equilibrium, temperature, gas constant and the ratio of specific heats are modelled as constants.
    
    
    Parameters
    ----------
    density_function : callable[[float], float]
        Function to retrieve the density at the current altitude.
    
    constant_temperature : float
        Constant atmospheric temperature.
    specific_gas_constant : float, default = 287.0
        Specific gas constant for (constant) atmospheric chemical composition.
    ratio_specific_heats : float, default = 1.4
        Ratio of specific heats for (constant) atmospheric chemical composition.
    Returns
    -------
    CustomConstantTemperatureAtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    with constant temperature and composition, but a density which varies with altitude according to a user-defined model:
    
    .. code-block:: python
    
      # Define the density as a function of altitude [in m]
      def density_function(h):
          # Return the density according to a modified exponential model
          return 1.15 * np.exp(-h/7300)
      # Define parameters for constant temperature and composition
      constant_temperature = 250.0
      specific_gas_constant = 300.0
      ratio_of_specific_heats = 1.4
      # Create the custom constant temperature atmosphere settings
      custom_density_settings = environment_setup.atmosphere.custom_constant_temperature(
          density_function,
          constant_temperature,
          specific_gas_constant,
          ratio_of_specific_heats)
      # Add the custom density to the body settings of "Earth"
      body_settings.get("Earth").atmosphere_settings = custom_density_settings"""

def custom_four_dimensional_constant_temperature(density_function: typing.Callable[[float, float, float, float], float], constant_temperature: float, specific_gas_constant: float=287.0, ratio_of_specific_heats: float=1.4) -> AtmosphereSettings:
    """Function for creating atmospheric model settings from custom density profile.
    
    Function for settings object, defining constant temperature atmosphere model from custom density profile.
    The user is specifying the density profile as a function of altitude, longitude, latitude and time.
    
    .. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
              The altitude is in meters, and the time is a Julian date in seconds since J2000.
    
    
    Parameters
    ----------
    density_function : callable[[float, float, float, float], float]
        Function to retrieve the density at the current altitude, longitude, latitude and time.
    
    constant_temperature : float
        Constant atmospheric temperature.
    specific_gas_constant : float, default = 287.0
        Specific gas constant for (constant) atmospheric chemical composition.
    ratio_specific_heats : float, default = 1.4
        Ratio of specific heats for (constant) atmospheric chemical composition.
    Returns
    -------
    CustomConstantTemperatureAtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    with constant temperature and composition (gas constant and ratio of specific heats), but a density which varies with altitude, longitude, latitude and time, according to a user-defined model:
    
    .. code-block:: python
    
      # Define the density as a function of altitude [m], longitude [deg], latitude [deg], and time [sec since J2000]
      def density_function(h, lon, lat, time):
          # Return the density according to an exponential model that varies with time to add noise with a sine (ignore lon/lat)
          return (1 + 0.15 * np.sin(time/10)) * np.exp(-h/7300)
      # Define the parameters for constant temperature and composition
      constant_temperature = 250.0
      specific_gas_constant = 300.0
      ratio_of_specific_heats = 1.4
      # Create the atmosphere settings and add to body settings of "Earth"
      body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature(
          density_function,
          constant_temperature,
          specific_gas_constant,
          ratio_of_specific_heats )"""

def custom_wind_model(wind_function: typing.Callable[[float, float, float, float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]], associated_reference_frame: aerodynamic_coefficients.AerodynamicsReferenceFrames=...) -> WindModelSettings:
    """Function for creating wind model settings with custom wind velocity.
    
    Function for settings object, defining wind model entirely from custom wind velocity function in a given reference frame.
    The custom wind velocity has to be given as a function of altitude, longitude, latitude and time.
    
    .. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
              The altitude is in meters, and the time is a Julian date in seconds since J2000.
    
    
    Parameters
    ----------
    wind_velocity : callable[[float, float, float, float], numpy.ndarray[numpy.float64[3, 1]]]
        Custom wind velocity function (w.r.t. altitude, longitude, latitude and time) in the specified reference frame.
    
    associated_reference_frame : dynamics.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
        Reference frame in which wind velocity is defined.
    
    Returns
    -------
    CustomWindModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.CustomWindModelSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`,
    using a user-defined wind-velocity function (of altitude, longitude, latitude and time), defined in a vertical aerodynamic reference frame:
    
    .. code-block:: python
    
      # Define the wind in 3 directions in the vertical reference frame
      def wind_function(h, lon, lat, time):
          # Meridional wind (pointing North) depends on latitude [deg] and time [sec since J2000]
          wind_Xv = lat*10/time
          # Zonal wind (pointing West) only depends on the longitude [deg]
          wind_Yv = 5/lon
          # Vertical wind (pointing out of the centre of the Earth) only depends on the altitude [m]
          wind_Zv = 1000/h
          # Return the custom wind
          return [wind_Xv, wind_Yv, wind_Zv]
      # Create the custom wind settings
      custom_wind = environment_setup.atmosphere.custom_wind_model(
          wind_function,
          environment.AerodynamicsReferenceFrames.vertical_frame)
      # Apply the custom wind settings to the Earth atmosphere settings
      body_settings.get("Earth").atmosphere_settings.wind_settings = custom_wind"""

def exponential(scale_height: float, surface_density: float, constant_temperature: float=288.15, specific_gas_constant: float=287.0, ratio_specific_heats: float=1.4) -> AtmosphereSettings:
    """Function for creating atmospheric model settings from fully parametrized exponential model.
    
    Function for settings object, defining exponential atmosphere model.
    The model is solely based on an exponentially decaying density profile with a constant temperature and composition
    (i.e. independent of time, latitude and longitude).
    
    The user has access to a fully parametrized model, meaning that in addition to the required input parameters ``scale_height`` and ``surface_density`` (ground-level air density),
    the user can specify non-standard values for constant temperature, gas constant and specific heats ratio.
    
    
    Parameters
    ----------
    scale_height : float
        Scale height for density profile of atmosphere.
    surface_density : float
        Atmospheric density at ground level.
    constant_temperature : float, default = 288.15
        Constant atmospheric temperature.
    specific_gas_constant : float, default = constants.SPECIFIC_GAS_CONSTANT_AIR
        Specific gas constant for (constant) atmospheric chemical composition.
    ratio_specific_heats : float, default = 1.4
        Ratio of specific heats for (constant) atmospheric chemical composition.
    Returns
    -------
    ExponentialAtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ExponentialAtmosphereSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    using the minimalist interface to the exponential model and taking parameters with classic values for Earth:
    
    .. code-block:: python
    
       # define parameters of an invariant exponential atmosphere model
       density_scale_height = 7.2E3
       constant_temperature = 290
       # create atmosphere settings and add to body settings of "Earth"
       body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.exponential(
           density_scale_height, density_at_zero_altitude)"""

def exponential_predefined(body_name: str) -> AtmosphereSettings:
    """Function for creating atmospheric model settings from pre-defined exponential model.
    
    Function for settings object, defining atmosphere model from pre-defined exponential model.
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
    
    
    Parameters
    ----------
    body_name : str
        Body for which pre-defined model settings are to be loaded. Available bodies "Earth", "Mars".
    
    Returns
    -------
    ExponentialAtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ExponentialAtmosphereSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Mars,
    using the interface of the predefined exponential model, using pre-encoded values:
    
    .. code-block:: python
    
       # Create atmosphere settings and add to body settings of "Mars"
       body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")"""

def nrlmsise00(space_weather_file: str='/home/dominic/.tudat/resource/space_weather/sw19571001.txt', use_storm_conditions: bool=False, use_anomalous_oxygen: bool=True) -> AtmosphereSettings:
    """Function for creating NRLMSISE-00 atmospheric model settings.
    
    Function for settings object, defining atmosphere model in accordance to the NRLMSISE-00 global reference model for Earth's atmosphere.
    The NRLMSISE-00 model implementation uses the code from `tudat-team/nrlmsise-00-cmake <https://github.com/tudat-team/nrlmsise-00-cmake>`_.
    
    
    Parameters
    ----------
    space_weather_file : str, default = :func:`~tudatpy.data.get_space_weather_path` + 'sw19571001.txt'
        File to be used for space weather characteristics as a function of time (e.g. F10.7, Kp, etc.). The file is typically taken from `celestrak <https://celestrak.org/SpaceData/sw19571001.txt>`_ (note that the file in your resources path will not be the latest version of this file; download and replace your existing file if required). Documentation on the file is given on the `celestrak website <https://celestrak.org/SpaceData/SpaceWx-format.php>`_
    use_storm_conditions : bool, default = false
        Boolean to define whether to use sub-daily Ap values when querying the NRLMSISE model, which is relevant under geomagnetic storm conditions (see `NRLMSISE code <https://github.com/tudat-team/nrlmsise-00-cmake/blob/master/nrlmsise-00.h>`_, setting this variable to true sets ``switches[9]`` to -1, with resulting details of Ap values defined in ``ap_array``).
    use_anomalous_oxygen : bool, default = true
        Boolean to define whether to use anomalous oxygen when querying the NRLMSISE model (if true, using ``gtd7d`` function, if false using ``gtd7`` function in `NRLMSISE code <https://github.com/tudat-team/nrlmsise-00-cmake/blob/master/nrlmsise-00.h>`_)
    
    Returns
    -------
    AtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    using the NRLMSISE-00 global reference model:
    
    .. code-block:: python
    
        # create atmosphere settings and add to body settings of body "Earth"
        body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.nrlmsise00()"""

def scaled_by_constant(unscaled_atmosphere_settings: AtmosphereSettings, density_scaling: float, is_scaling_absolute: bool=False) -> AtmosphereSettings:
    """Function for creating scaled atmospheric model settings.
    
    Function for settings object, defining atmospheric model based on an scaling of an existing atmospheric settings object.
    The user can apply a scaling factor (or an absolute value) to the air densities of the existing model settings (for instance for an uncertainty analysis).
    
    
    Parameters
    ----------
    unscaled_atmosphere_settings : AtmosphereSettings
        Sets base settings of atmosphere model to be scaled.
    density_scaling : float
        Constant scaling factor to be applied to the entire air density profile.
    is_scaling_absolute : bool, default=false
        Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.
    
    Returns
    -------
    ScaledAtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ScaledAtmosphereSettings` class.
    
    
    
    Notes
    -----
    At present, the scaled atmosphere model only supports scaling of the density value.
    For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
    this calculation is performed using the `unscaled` density!
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    by modifying an existing :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` object such that the resulting air density profile is scaled by a constant:
    
    .. code-block:: python
    
      # define parameter for scaling
      scaling_constant = 1.5
      # define variable containing the existing atmosphere model settings
      unscaled_atmosphere_settings = body_settings.get( "Earth" ).atmosphere_settings
      # create atmosphere settings and add to body settings of "Earth"
      body_settings.get( "Earth" ).atmosphere_settings =  environment_setup.atmosphere.scaled_by_constant(
          unscaled_atmosphere_settings,
          scaling_constant )"""

def scaled_by_function(unscaled_atmosphere_settings: AtmosphereSettings, density_scaling_function: typing.Callable[[float], float], is_scaling_absolute: bool=False) -> AtmosphereSettings:
    """Function for creating scaled atmospheric model settings.
    
    Function for settings object, defining atmospheric model based on scaling an existing atmospheric settings object.
    The user can apply custom scaling factors (or absolute values) to the air densities of the existing model settings (for instance for an uncertainty analysis).
    
    
    Parameters
    ----------
    unscaled_atmosphere_settings : AtmosphereSettings
        Sets base settings of atmosphere model to be scaled.
    density_scaling_function : Callable[[float], float]
        Specifies air density scaling factor as a function of time.
    is_scaling_absolute : bool, default=false
        Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.
    
    Returns
    -------
    ScaledAtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ScaledAtmosphereSettings` class.
    
    
    
    Notes
    -----
    At present, the scaled atmosphere model only supports scaling of the density value.
    For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
    this calculation is performed using the `unscaled` density!
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    by modifying an existing :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` object such, that the resulting air density profile is scaled with a user-defined function of time:
    
    .. code-block:: python
    
      # Define the density scaling as a function of time [sec since J2000] (to add noise with a sine)
      def scaling_function(time):
          return 1 + np.sin(time / 50) * 0.25
      # Extract the existing atmosphere model settings
      unscaled_atmosphere_settings = body_settings.get( "Earth" ).atmosphere_settings
      # Create the atmosphere settings and add to body settings of "Earth"
      body_settings.get( "Earth" ).atmosphere_settings =  environment_setup.atmosphere.scaled_by_function(
          unscaled_atmosphere_settings,
          scaling_function )"""

def tabulated(atmosphere_data_file: str, dependent_variable_names: list[AtmosphereDependentVariables]=..., specific_gas_constant: float=287.0, ratio_of_specific_heats: float=1.4) -> AtmosphereSettings:
    ...

def us76() -> AtmosphereSettings:
    """Function for creating US76 standard atmosphere model settings.
    
    Function for creating US76 standard atmosphere model settings. The model is defined using tabulated data for density, pressure and temperature,
    from an altitude of -5 km up to 1000 km. Up to 100 km, a data point is provided every 100 m. Above 100 km, a data point is provided every 1 km. The data
    are interpolated using a cubic spline interpolator. Note that this model is specific to Earth's atmosphere.
    
    Returns
    -------
    AtmosphereSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
    using the US76 standard atmosphere model:
    
    .. code-block:: python
    
       # create atmosphere settings and add to body settings of body "Earth"
       body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.us76()"""
tabulated_density: AtmosphereDependentVariables
tabulated_gas_constant: AtmosphereDependentVariables
tabulated_molar_mass: AtmosphereDependentVariables
tabulated_pressure: AtmosphereDependentVariables
tabulated_specific_heat_ratio: AtmosphereDependentVariables
tabulated_temperature: AtmosphereDependentVariables