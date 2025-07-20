import numpy
import pybind11_stubgen.typing_ext
from ....dynamics import environment
from ....math import interpolators
import typing
__all__ = ['GlobalIonosphereModelVtecCalculator', 'JakowskiVtecCalculator', 'LightTimeConvergenceCriteria', 'LightTimeCorrectionSettings', 'LightTimeFailureHandling', 'TroposphericMappingModel', 'VtecCalculator', 'WaterVaporPartialPressureModel', 'accept_without_warning', 'approximated_second_order_relativistic_light_time_correction', 'bean_and_dutton', 'dsn_tabulated_ionospheric_light_time_correction', 'dsn_tabulated_tropospheric_light_time_correction', 'first_order_relativistic_light_time_correction', 'inverse_power_series_solar_corona_light_time_correction', 'ionex_ionospheric_light_time_correction', 'jakowski_ionospheric_light_time_correction', 'light_time_convergence_settings', 'niell', 'print_warning_and_accept', 'saastamoinen_tropospheric_light_time_correction', 'set_ionosphere_model_from_ionex', 'set_vmf_troposphere_data', 'simplified_chao', 'tabulated', 'throw_exception', 'vmf3_tropospheric_light_time_correction']

class GlobalIonosphereModelVtecCalculator(VtecCalculator):

    def __init__(self, ionosphere_model: environment.IonosphereModel) -> None:
        ...

    def calculate_vtec(self, time: float, sub_ionospheric_point: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]) -> float:
        ...

class JakowskiVtecCalculator(VtecCalculator):

    def __init__(self, sun_declination_function: typing.Callable[[float], float], f10p7_function: typing.Callable[[float], float], use_utc_time_for_local_time: bool=False) -> None:
        ...

    def calculate_vtec(self, time: float, sub_ionospheric_point: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]) -> float:
        ...

class LightTimeConvergenceCriteria:
    """Base class to define criteria of light time convergence.
    
    Base class to define criteria of light time convergence.
    This class is not used for calculations of corrections, but is used for the purpose of defining the light time convergence criteria.
    Specific light time convergence criteria must be defined using an object derived from this class.
    Instances of this class are typically created via the :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings` function.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of a LightTimeConvergenceCriteria object
        from tudatpy.estimation.observable_models_setup import light_time_corrections
    
        # Create Default Light Time Convergence Settings (no args specified = setting default arguments)
        light_time_convergence_settings = light_time_corrections.light_time_convergence_settings()
    
        # Show that it is an LightTimeConvergenceCriteria object.
        print(light_time_convergence_settings)"""

class LightTimeCorrectionSettings:
    """Base class to define light time correction settings.
    
    Base class to define light time correction settings.
    This class is not used for calculations of corrections, but is used for the purpose of defining the light time correction properties.
    Specific light time correction settings must be defined using an object derived from this class.
    
    Instances of this class are typically created via the
    :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction` function
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of a LightTimeCorrectionSettings object
        from tudatpy.estimation.observable_models_setup import light_time_corrections, links
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # Case 1: perturbing body (Earth) involved in the observations
        # In this case, Earth is a receiver, so the body’s state will be evaluated at the reception time.
        perturbing_body = ['Earth']
        doppler_observation_settings = light_time_corrections.first_order_relativistic_light_time_correction(perturbing_body)
    
        # Show that it is a LightTimeCorrectionSettings object.
        print(doppler_observation_settings)
    
        # Case 2: perturbing body (Sun) not involved in the observations
        # In this case, the body's state will be evaluated at the midpoint time between the transmission and reception events.
        perturbing_body = ['Sun']
    
        # Use: light_time_corrections.first_order_relativistic_light_time_correction to create a LightTimeCorrectionSettings object
        # Note: first_order_relativistic_light_time_correction only requires the perturbing list of bodies to be passed as arguments
        doppler_observation_settings = light_time_corrections.first_order_relativistic_light_time_correction(perturbing_body)
    
        # Show that it is an LightTimeCorrectionSettings object.
        print(doppler_observation_settings.transmitter_proper_time_rate_settings)
        print(dir(doppler_observation_settings))"""

class LightTimeFailureHandling:
    """Enumeration of behaviour when failing to converge light-time with required settings.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to print all available Light Time Failure Handling Types
        from tudatpy.estimation.observable_models_setup import light_time_corrections
    
        num_LightTimeFailureHandling_types = len(light_time_corrections.LightTimeFailureHandling.__members__)
        print(f'The length of all available Tudatpy Light Time Failure Handling Types is: {num_LightTimeFailureHandling_types}')
    
        # Print all available Observation Viability Types using the "name" property
        for i in range(num_LightTimeFailureHandling_types):
            print(i, light_time_corrections.LightTimeFailureHandling(i).name)
    
    
    
          
    
    Members:
    
      accept_without_warning
    
      print_warning_and_accept
    
      throw_exception"""
    __members__: typing.ClassVar[dict[str, LightTimeFailureHandling]]
    accept_without_warning: typing.ClassVar[LightTimeFailureHandling]
    print_warning_and_accept: typing.ClassVar[LightTimeFailureHandling]
    throw_exception: typing.ClassVar[LightTimeFailureHandling]

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

class TroposphericMappingModel:
    """No documentation found.
    
    Members:
    
      simplified_chao
    
      niell"""
    __members__: typing.ClassVar[dict[str, TroposphericMappingModel]]
    niell: typing.ClassVar[TroposphericMappingModel]
    simplified_chao: typing.ClassVar[TroposphericMappingModel]

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

class VtecCalculator:
    pass

class WaterVaporPartialPressureModel:
    """No documentation found.
    
    Members:
    
      tabulated
    
      bean_and_dutton"""
    __members__: typing.ClassVar[dict[str, WaterVaporPartialPressureModel]]
    bean_and_dutton: typing.ClassVar[WaterVaporPartialPressureModel]
    tabulated: typing.ClassVar[WaterVaporPartialPressureModel]

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

@typing.overload
def approximated_second_order_relativistic_light_time_correction(perturbing_bodies: list[str]) -> LightTimeCorrectionSettings:
    """Function for creating settings for Moyer, 2000 Eq 8.55 approximated second-order relativistic light-time corrections.
    
    Function for creating settings for approximated second-order relativistic light-time corrections:  These corrections account for the delay in light travel time caused by stationary point masses, calculated up to
    :math:`c^{-2}` according to general relativity ( Moyer, 2000 Eq 8.55; correction term for Sun) and it includes the bending of light due to the perturbing body. A key consideration in the model is the time at which the states of the perturbing bodies are evaluated. This depends on their involvement in the observation link ends:
    
    * 1. **Perturbing Body as a Link End:** If the perturbing body (e.g., Earth) is directly involved in the observation (e.g., as the location of a transmitter or receiver):
    
        - The body's state is evaluated at the **transmission time** if it acts as the **transmitter**.
        - The body's state is evaluated at the **reception time** if it acts as the **receiver**.
    
    * 2. **Perturbing Body Not as a Link End:** If the perturbing body is not part of the observation link ends, its state is evaluated at the **midpoint time** between the transmission and reception events.
    
    Parameters
    ----------
    perturbing_bodies : List[str]
        A list containing the names of the bodies due to which the light-time correction is to be taken into account.
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` configured to include
        approximated second-order relativistic light-time corrections.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the first_order_relativistic_light_time_correction function
        from tudatpy.estimation.observable_models_setup import light_time_corrections, links
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # The function first_order_relativistic_light_time_correction() requires a list of strings (perturbing body/bodies) as input
        perturbing_body = ['Earth']
        doppler_observation_settings = light_time_corrections.approximated_second_order_relativistic_light_time_correction(perturbing_body)
    
        # Show that it returns a LightTimeCorrectionSettings object.
        print(doppler_observation_settings)"""

@typing.overload
def approximated_second_order_relativistic_light_time_correction(perturbing_bodies: list[str]) -> LightTimeCorrectionSettings:
    """Function for creating settings for Moyer, 2000 Eq 8.55 approximated second-order relativistic light-time corrections.
    
    Function for creating settings for approximated second-order relativistic light-time corrections:  These corrections account for the delay in light travel time caused by stationary point masses, calculated up to
    :math:`c^{-2}` according to general relativity ( Moyer, 2000 Eq 8.55) and it includes the bending of light due to the perturbing body. A key consideration in the model is the time at which the states of the perturbing bodies are evaluated. This depends on their involvement in the observation link ends:
    
    * 1. **Perturbing Body as a Link End:** If the perturbing body (e.g., Earth) is directly involved in the observation (e.g., as the location of a transmitter or receiver):
    
        - The body's state is evaluated at the **transmission time** if it acts as the **transmitter**.
        - The body's state is evaluated at the **reception time** if it acts as the **receiver**.
    
    * 2. **Perturbing Body Not as a Link End:** If the perturbing body is not part of the observation link ends, its state is evaluated at the **midpoint time** between the transmission and reception events.
    
    Parameters
    ----------
    perturbing_bodies : List[str]
        A list containing the names of the bodies due to which the light-time correction is to be taken into account.
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` configured to include
        approximated second-order relativistic light-time corrections.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the first_order_relativistic_light_time_correction function
        from tudatpy.estimation.observable_models_setup import light_time_corrections, links
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # The function first_order_relativistic_light_time_correction() requires a list of strings (perturbing body/bodies) as input
        # and a boolean value for bending (default is True).
        perturbing_body = ['Earth']
        doppler_observation_settings = light_time_corrections.approximated_second_order_relativistic_light_time_correction(perturbing_body)
    
        # Show that it returns a LightTimeCorrectionSettings object.
        print(doppler_observation_settings)"""

def dsn_tabulated_ionospheric_light_time_correction(file_names: list[str], spacecraft_name_per_id: dict[int, str]={}, quasar_name_per_id: dict[int, str]={}, reference_frequency: float=2295000000.0, body_with_atmosphere_name: str='Earth') -> LightTimeCorrectionSettings:
    """No documentation found."""

def dsn_tabulated_tropospheric_light_time_correction(file_names: list[str], body_with_atmosphere_name: str='Earth', mapping_model: TroposphericMappingModel=...) -> LightTimeCorrectionSettings:
    """No documentation found."""

def first_order_relativistic_light_time_correction(perturbing_bodies: list[str]) -> LightTimeCorrectionSettings:
    """Function for creating settings for first-order relativistic light-time corrections.
    
    Function for creating settings for first-order relativistic light-time corrections:  These corrections account for the delay in light travel time caused by stationary point masses, calculated up to
    :math:`c^{-2}` according to general relativity (e.g., Moyer, 2000 Eq 8.55). A key consideration in the model is the time at which the states of the perturbing bodies are evaluated. This depends on their involvement in the observation link ends:
    
    * 1. **Perturbing Body as a Link End:** If the perturbing body (e.g., Earth) is directly involved in the observation (e.g., as the location of a transmitter or receiver):
    
        - The body's state is evaluated at the **transmission time** if it acts as the **transmitter**.
        - The body's state is evaluated at the **reception time** if it acts as the **receiver**.
    
    * 2. **Perturbing Body Not as a Link End:** If the perturbing body is not part of the observation link ends, its state is evaluated at the **midpoint time** between the transmission and reception events.
    
    Parameters
    ----------
    perturbing_bodies : List[str]
        A list containing the names of the bodies due to which the light-time correction is to be taken into account.
    
    Returns
    -------
    :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` configured to include
        first-order relativistic light-time corrections.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the first_order_relativistic_light_time_correction function
        from tudatpy.estimation.observable_models_setup import light_time_corrections, links
    
        # Create Link Ends dictionary
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create a Link Definition Object from link_ends dictionary
        Link_Definition_Object = links.LinkDefinition(link_ends)
    
        # The function first_order_relativistic_light_time_correction() requires a list of strings (perturbing body/bodies) as input
        # and a boolean value for bending (default is True).
        perturbing_body = ['Earth']
        doppler_observation_settings = light_time_corrections.first_order_relativistic_light_time_correction(perturbing_body)
    
        # Show that it returns a LightTimeCorrectionSettings object.
        print(doppler_observation_settings)"""

def inverse_power_series_solar_corona_light_time_correction(coefficients: list[float]=[130000000000000.0, 500000000000.0], positive_exponents: list[float]=[6.0, 2.0], delay_coefficient: float=40.3, sun_body_name: str='Sun') -> LightTimeCorrectionSettings:
    """No documentation found."""

def ionex_ionospheric_light_time_correction(body_with_ionosphere_name: str, ionosphere_height: float, first_order_delay_coefficient: float=40.3) -> LightTimeCorrectionSettings:
    """Create IONEX-based ionospheric light time correction settings."""

def jakowski_ionospheric_light_time_correction(ionosphere_height: float=400000.0, first_order_delay_coefficient: float=40.3, solar_activity_data_path: str='/home/dominic/.tudat/resource/space_weather/sw19571001.txt', geomagnetic_pole_latitude: float=1.4119713648634127, geomagnetic_pole_longitude: float=-1.2671090369478832, use_utc_for_local_time_computation: bool=False, body_with_atmosphere_name: str='Earth') -> LightTimeCorrectionSettings:
    """No documentation found."""

def light_time_convergence_settings(iterate_corrections: bool=False, maximum_number_of_iterations: int=50, absolute_tolerance: float=..., failure_handling: LightTimeFailureHandling=...) -> LightTimeConvergenceCriteria:
    """Function for creating convergence settings for solving the light-time equation.
    
    Function for creating convergence settings for solving the light-time equation. Computing the light time
    :math:`s=t_{R}-t_{T}` between two receiver :math:`R` and transmitter :math:`T` requires the iterative
    solution of the following equation:
    
    .. math::
        t_{R} - t_{T} = c\\left(|\\mathbf{r}_{R}(t_{R}) - \\mathbf{r}_{T}(t_{T})| + \\Delta s(t_{R}, t_{T}, \\mathbf{r}_{R}(t_{R}), \\mathbf{r}_{T}(t_{T}))\\right)
    
    where either the reception time :math:`t_{R}` or the transmission time :math:`t_{T}` is kept fixed (reference link end time). The term :math:`\\Delta s` contains any
    deviations in the light-time from straight-line propagation at speed of light (relativistic corrections, media corrections, etc.). The algorithm starts
    at :math:`t_{R}=t_{T}`, and uses this to evaluate the right-hand side of the above equation. This leads to a new value of :math:`t_{R}` or :math:`t_{T}` (depending on which is kept fixed)
    and the right-hand side is re-evaluated in a new iteration. The input to this function defines the settings for when the iteration will terminate.
    
    Parameters
    ----------
    iterate_corrections : bool, default = False
        Boolean denoting whether the terms :math:`\\Delta s` are recomputed at each iteration or not. If false, the corrections are calculated only on the first iteration. Afterwards, the value
        is kept fixed until convergence. Once preliminarily converged, the algorithm recomputes :math:`\\Delta s`, and continues the iteration (until proper convergence) while now recomputing
        :math:`\\Delta s` each iteration. Setting this input to false is typically safe, and is computationally more efficient.
    
    maximum_number_of_iterations : int, default = 50
        Maximum number of iterations taken by the algorithm. If this number of iterations is reached without convergence (as defined by ``absolute_tolerance`` input),
        the behaviour of the algorithm is defined by the ``failure_handling`` input.
    
    absolute_tolerance : float, default = nan
        Difference in :math:`t_{R}-t_{T}` between two consecutive iterations below which the algorithm is considered to be converged. Default value is nan, which means the default value is taken.
        The default value depends on the time representation used (1 ps for float; 1 fs for Time class)
    
    failure_handling : LightTimeFailureHandling, default = accept_without_warning
        Input defines behaviour when failing to converge within the required number of iterations. NOTE: the default value should be overridden for high-accuracy applications
    
    Returns
    -------
    :class:`LightTimeConvergenceCriteria`
        Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria` with the required settings.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the light_time_convergence_settings function
        from tudatpy.estimation.observable_models_setup import light_time_corrections
    
        # The light_time_convergence_settings function can be used with default inputs as just:
        light_time_convergence_settings = light_time_corrections.light_time_convergence_settings()
        # A LightTimeConvergenceCriteria object is returned
        print(light_time_convergence_settings)
    
        # Users can also specify the following input arguments:
        # iterate_corrections, maximum_number_of_iterations, absolute_tolerance, failure_handling.
        # Let's set the failure_handling argument to LightTimeFailureHandling.print_warning_and_accept (default was LightTimeFailureHandling.accept_without_warning)
        light_time_convergence_settings = light_time_corrections.light_time_convergence_settings(
            failure_handling = light_time_corrections.LightTimeFailureHandling.print_warning_and_accept
        )
        # Again, a LightTimeConvergenceCriteria object is returned
        print(light_time_convergence_settings)"""

def saastamoinen_tropospheric_light_time_correction(body_with_atmosphere_name: str='Earth', mapping_model: TroposphericMappingModel=..., water_vapor_partial_pressure_model: WaterVaporPartialPressureModel=...) -> LightTimeCorrectionSettings:
    """No documentation found."""

def set_ionosphere_model_from_ionex(data_files: list[str], bodies: environment.SystemOfBodies, interpolator_settings: interpolators.InterpolatorSettings=None) -> None:
    ...

def set_vmf_troposphere_data(data_files: list[str], file_has_meteo: bool, file_has_gradient: bool, bodies: environment.SystemOfBodies, set_troposphere_data: bool=True, set_meteo_data: bool=True, interpolator_settings: interpolators.InterpolatorSettings=...) -> None:
    ...

def vmf3_tropospheric_light_time_correction(body_with_atmosphere_name: str='Earth', use_gradient_correction: bool=True, tropospheric_mapping_model: TroposphericMappingModel=...) -> LightTimeCorrectionSettings:
    """Create VMF3 tropospheric light time correction settings."""
accept_without_warning: LightTimeFailureHandling
bean_and_dutton: WaterVaporPartialPressureModel
niell: TroposphericMappingModel
print_warning_and_accept: LightTimeFailureHandling
simplified_chao: TroposphericMappingModel
tabulated: WaterVaporPartialPressureModel
throw_exception: LightTimeFailureHandling