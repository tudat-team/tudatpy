import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['BodyPanelReflectionLawSettings', 'CannonBallRadiationPressureInterfaceSettings', 'KnockeTypeSurfacePropertyDistributionModel', 'LuminosityModelSettings', 'PanelRadiosityModelSettings', 'RadiationPressureInterfaceSettings', 'RadiationPressureTargetModelSettings', 'RadiationPressureTargetModelType', 'RadiationPressureType', 'RadiationSourceModelSettings', 'SphericalHarmonicsSurfacePropertyDistributionModel', 'SurfacePropertyDistributionSettings', 'albedo_dlam1', 'albedo_knocke', 'cannonball', 'cannonball_radiation_pressure_interface', 'cannonball_radiation_target', 'cannonball_target', 'constant_albedo_surface_radiosity', 'constant_luminosity', 'constant_radiosity', 'constant_surface_property_distribution', 'custom', 'custom_surface_property_distribution', 'emissivity_knocke', 'irradiance_based_constant_luminosity', 'irradiance_based_time_variable_luminosity', 'isotropic_radiation_source', 'knocke_type_surface_property_distribution', 'lambertian_body_panel_reflection', 'multi_type_target', 'paneled_target', 'panelled_extended_radiation_source', 'panelled_radiation_target', 'predefined_knocke_type_surface_property_distribution', 'predefined_spherical_harmonic_surface_property_distribution', 'specular_diffuse_body_panel_reflection', 'spherical_harmonic_surface_property_distribution', 'thermal_emission_angle_based_radiosity', 'thermal_emission_blackbody_constant_emissivity', 'thermal_emission_blackbody_variable_emissivity', 'time_variable_luminosity', 'undefined_target', 'variable_albedo_surface_radiosity']

class BodyPanelReflectionLawSettings:
    """Base class for providing settings for body panel reflection law models, to be used for defining spacecraft surface properties relevant for the computation of radiation pressure acting on a macromodel."""

class CannonBallRadiationPressureInterfaceSettings(RadiationPressureInterfaceSettings):
    """No documentation found."""

class KnockeTypeSurfacePropertyDistributionModel:
    """Members:
    
      custom
    
      albedo_knocke
    
      emissivity_knocke"""
    __members__: typing.ClassVar[dict[str, KnockeTypeSurfacePropertyDistributionModel]]
    albedo_knocke: typing.ClassVar[KnockeTypeSurfacePropertyDistributionModel]
    custom: typing.ClassVar[KnockeTypeSurfacePropertyDistributionModel]
    emissivity_knocke: typing.ClassVar[KnockeTypeSurfacePropertyDistributionModel]

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

class LuminosityModelSettings:
    """Base class for providing settings for body source luminosity settings, to be used (typically but not necessarily) for defining the Sun's luminosity."""

class PanelRadiosityModelSettings:
    """Base class for providing settings for body panel radiosity models, to be used (typically but not necessarily) for defining surface radiosity of a panelled solar system body as a result of albedo and/or planetary radiation pressure"""

class RadiationPressureInterfaceSettings:
    """No documentation found."""

class RadiationPressureTargetModelSettings:
    """Base class for providing settings for properties of a radiation target (e.g. spacecraft), to be used in the context of (for instance) calculation of radiation pressure on spacecraft"""

class RadiationPressureTargetModelType:
    """No documentation found.
    
    Members:
    
      cannonball_target : No documentation found.
    
      paneled_target : No documentation found.
    
      multi_type_target : No documentation found.
    
      undefined_target : No documentation found."""
    __members__: typing.ClassVar[dict[str, RadiationPressureTargetModelType]]
    cannonball_target: typing.ClassVar[RadiationPressureTargetModelType]
    multi_type_target: typing.ClassVar[RadiationPressureTargetModelType]
    paneled_target: typing.ClassVar[RadiationPressureTargetModelType]
    undefined_target: typing.ClassVar[RadiationPressureTargetModelType]

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

class RadiationPressureType:
    """No documentation found.
    
    Members:
    
      cannonball_radiation_pressure_interface : No documentation found."""
    __members__: typing.ClassVar[dict[str, RadiationPressureType]]
    cannonball_radiation_pressure_interface: typing.ClassVar[RadiationPressureType]

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

class RadiationSourceModelSettings:
    """Base class for providing settings for properties of a radiation source (e.g. Sun), to be used in the context of (for instance) calculation of radiation pressure on spacecraft"""

class SphericalHarmonicsSurfacePropertyDistributionModel:
    """Members:
    
      albedo_dlam1"""
    __members__: typing.ClassVar[dict[str, SphericalHarmonicsSurfacePropertyDistributionModel]]
    albedo_dlam1: typing.ClassVar[SphericalHarmonicsSurfacePropertyDistributionModel]

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

class SurfacePropertyDistributionSettings:
    """Base class for providing settings for body surface property distribution settings, to be used (typically but not necessarily) for defining surface distribution of albedo and emissivity of solar system bodies for calculations of albedo and planetary radiation pressure.Note that not all albedo/emissivity models require this type of distribution model"""

def cannonball(source_body: str, reference_area: float, radiation_pressure_coefficient: float, occulting_bodies: list[str]=[]) -> RadiationPressureInterfaceSettings:
    """No documentation found."""

def cannonball_radiation_target(reference_area: float, radiation_pressure_coefficient: float, per_source_occulting_bodies: dict[str, list[str]]={}) -> RadiationPressureTargetModelSettings:
    """Function for cannonball radiation target
    
    
    Parameters
    ----------
    reference_area : float
        Cross-sectional area of cannonball [:math:`m^{2}`]
    radiation_pressure_coefficient : float
        Radiation pressure coefficient [-]
    per_source_occulting_bodies : Dict[str, List[str]]
        Names of bodies to occult the source as seen from this target
    Returns
    -------
    CannonballRadiationPressureTargetModelSettings
        Object defining settings for a cannonball radiation pressure target model"""

def constant_albedo_surface_radiosity(constant_albedo: float, original_source_name: str) -> PanelRadiosityModelSettings:
    """Function for creating settings for surface constant albedo surface radiosity of an extended source
    
    Function for creating settings for surface radiosity of an extended source, with surface radiation the result
    of albedo using a Lambertian scattering law, and a constant albedo value over the surface.
    For a surface panel normal of :math:`\\hat{\\mathbf{n}}`, a vector :math:`\\mathbf{r}` from the surface element to the target, and a
    vector :math:`\\mathbf{r}_{s}` from the surface element to the original source (typically the Sun),
    the resulting irradiance :math:`\\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):
    
    .. math::
       \\Phi=\\cos\\theta_{s}\\Phi_{s}\\frac{a}{\\pi}\\frac{A\\cos\\theta}{\\pi ||\\mathbf{r}||^{2}}
    
    with :math:`\\theta_{s}` the angle between :math:`\\hat{\\mathbf{n}}` and :math:`\\mathbf{r_{s}}`, :math:`\\Phi_{s}` the irradiance from the original source
    at the panel of the reflecting body, :math:`a` is the albedo coefficient, :math:`A` the panel area, :math:`\\theta` is the angle between :math:`\\hat{\\mathbf{n}}` and :math:`\\mathbf{r}`.
    
    
    Parameters
    ----------
    constant_albedo : float
        Constant value of the albedo coefficient :math:`a`.
    original_source_name : str
        Name of the original source from which the radiation is reflection to the target.
    Returns
    -------
    PanelRadiosityModelSettings
        Object defining settings for source panel radiosity"""

def constant_luminosity(luminosity: float) -> LuminosityModelSettings:
    """Function for creating constant radiation source luminosity settings.
    
    Function for creating constant radiation source luminosity settings, defining the total
    radiated power (in Watts) of a given source. With this function, the source luminosity is constant,
    and is assumed to emit radiation isotropically.
    
    
    Parameters
    ----------
    luminosity : float
        Constant source luminosity (in Watt)
    Returns
    -------
    LuminosityModelSettings
        Object defining settings for source luminosity"""

def constant_radiosity(radiosity: float) -> PanelRadiosityModelSettings:
    """Function for creating settings for surface constant surface radiosity of an extended source
    
    Function for creating settings for surface radiosity of an extended source, using constant Lambertian radiosity :math:`J` (in :math:`W/m^{2}`).
    For a surface panel normal of :math:`\\hat{\\mathbf{n}}` and a vector :math:`\\mathbf{r}` from the surface element to the target, the resulting
    irradiance :math:`\\Phi` (in :math:`W/m^{2}`) at the target is (if :math:`\\theta>0`, or in other words if the panel is visible from the target):
    
    .. math::
       \\Phi=J\\frac{A\\cos\\theta}{\\pi ||\\mathbf{r}||^{2}}
    
    with :math:`A` the panel area, :math:`\\theta` is the angle between :math:`\\hat{\\mathbf{n}}` and :math:`\\mathbf{r}`.
    
    
    Parameters
    ----------
    radiosity : float
        Constant Lambertian radiosity from surface in :math:`W/m^{2}`.
    Returns
    -------
    PanelRadiosityModelSettings
        Object defining settings for source panel radiosity"""

def constant_surface_property_distribution(constant_value: float) -> SurfacePropertyDistributionSettings:
    """Function for creating constant radiative surface property distribution settings.
    
    Function for creating constant radiative surface property (e.g. albedo, emissivity, etc.) distribution settings.
    
    
    Parameters
    ----------
    constant_value : float
        Constant surface property value
    Returns
    -------
    SurfacePropertyDistributionSettings
        Object defining settings for surface property distribution"""

def custom_surface_property_distribution(custom_function: typing.Callable[[float, float, float], float]) -> SurfacePropertyDistributionSettings:
    """Function for creating radiative surface property distribution settings according to a custom user-defined model.
    
    Function for creating radiative surface property (e.g. albedo, emissivity, etc.) distribution settings
    according to a custom user-defined model, as a function of latitude, longitude and time.
    
    
    Parameters
    ----------
    custom_function : Callable[[float, float, float], float]
        Function providing surface property as a function of latitude, longitude and time (in that order).
    Returns
    -------
    SurfacePropertyDistributionSettings
        Object defining settings for surface property distribution"""

def irradiance_based_constant_luminosity(constant_irradiance: float, reference_distance: float) -> LuminosityModelSettings:
    """Function for creating source luminosity settings based on the irradiance at a reference distance.
    
    Function for creating source luminosity based on the irradiance at a reference distance. For instance,
    one can provide the solar irradiance at 1 AU, and this will be translated to the Sun's luminosity. With this function,
    the source luminosity is constant, and is assumed to emit radiation isotropically.
    
    
    Parameters
    ----------
    constant_irradiance : float
        Irradiance at reference distance from center of source (in :math:`W/m^{2}`)
    reference_distance : float
        Distance from center of source at which the irradiance is defined
    Returns
    -------
    LuminosityModelSettings
        Object defining settings for source luminosity
    
    Examples
    --------
    In this example, we create a constant luminosity model for the Sun, based on the solar constant, which defines the irradiance at 1 AU. Different values exist in literature for the solar constant, in this case we assume 1367 W/m^2 at 1 AU (Vallado 2013).
    Assuming an isotropic radiation source, we then create the radiation source settings for the Sun.
    
    .. code-block:: python
    
        ...
    
        from tudatpy import constants
    
        irradiance_at_1AU = 1367.0  # W/m^2, Vallado 2013
    
        luminosity_model_settings = (
            environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
                irradiance_at_1AU, constants.ASTRONOMICAL_UNIT
            )
        )
        radiation_source_settings_sun = (
            environment_setup.radiation_pressure.isotropic_radiation_source(
                luminosity_model_settings
            )
        )
    
        body_settings.get("Sun").radiation_source_settings = radiation_source_settings_sun"""

def irradiance_based_time_variable_luminosity(irradiance_function: typing.Callable[[float], float], reference_distance: float) -> LuminosityModelSettings:
    """Function for creating time-variable source luminosity settings based on the irradiance at a reference distance.
    
    Function for creating source time-variable luminosity based on the irradiance at a reference distance. For instance,
    one can provide the solar irradiance at 1 AU as a function of time, and this will be translated to the Sun's luminosity.
    With this function, the source is assumed to emit radiation isotropically.
    
    
    Parameters
    ----------
    irradiance_function : Callable[[float], float]
        Function returning irradiance at reference distance from center of source (in :math:`W/m^{2}`) as a function fo time
    reference_distance : float
        Distance from center of source at which the irradiance is defined
    Returns
    -------
    LuminosityModelSettings
        Object defining settings for source luminosity"""

def isotropic_radiation_source(luminosity_model: LuminosityModelSettings) -> RadiationSourceModelSettings:
    """Function for creating settings for an isotropic radiation source
    
    Function for creating settings for a radiation source that emits isotropically. The source is provided
    with a luminosity model :math:`L(t)` as a (possible) function of time :math:`t`. The irradiance :math:`\\Phi` at a relative position
    :math:`\\mathbf{r}` from the source's center is then computed from:
    
    .. math::
       \\Phi=\\frac{L}{4\\pi||\\mathbf{r}||^{2}}
    
    
    Parameters
    ----------
    luminosity_model : LuminosityModelSettings
        Settings for the luminosity model.
    Returns
    -------
    RadiationSourceModelSettings
        Object defining settings for source model irradiance
    Examples
    --------
    In this example, we create a constant luminosity model for the Sun, based on the solar constant, which defines the irradiance at 1 AU. Different values exist in literature for the solar constant, in this case we assume 1367 W/m^2 at 1 AU (Vallado 2013).
    Assuming an isotropic radiation source, we then create the radiation source settings for the Sun.
    
    .. code-block:: python
    
        ...
    
        from tudatpy import constants
    
        irradiance_at_1AU = 1367.0  # W/m^2, Vallado 2013
    
        luminosity_model_settings = (
            environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
                irradiance_at_1AU, constants.ASTRONOMICAL_UNIT
            )
        )
        radiation_source_settings_sun = (
            environment_setup.radiation_pressure.isotropic_radiation_source(
                luminosity_model_settings
            )
        )
    
        body_settings.get("Sun").radiation_source_settings = radiation_source_settings_sun"""

def knocke_type_surface_property_distribution(constant_contribution: float, constant_degree_one_contribution: float, cosine_periodic_degree_one_contribution: float, sine_periodic_degree_one_contribution: float, constant_degree_two_contribution: float, reference_epoch: float, period: float) -> ...:
    """Function for creating radiative surface property distribution settings according to 'Knocke-type' model
    
    Function for creating radiative surface property (e.g. albedo, emissivity, etc.) distribution settings
    according to a model such as the one used by Knocke (1988). This model uses a degree two zonal spherical harmonic model, with
    a sinusoidal variation in the degree one coefficient. The surface property :math:`k` is computed from:
    
    .. math::
       k(\\phi,t)=a_{0}+a_{1}P_{1}(\\sin\\phi)+a_{2}P_{2}(\\sin\\phi)
    .. math::
       a_{1}=c_{0}+c_{1}\\cos\\left(\\frac{2\\pi(t-t_{0})}{T}\\right)+c_{2}\\sin\\left(\\frac{2\\pi(t-t_{0})}{T}\\right)
    
    with the angle :math:`\\phi` denotes the body-fixed latitude of the evaluation point, and :math:`t`, :math:`t_{0}` and :math:`T` define the current time,
    reference time and period of the variation, respectively. The coefficients :math:`a_{0}, a_{2}, c_{0}, c_{1}, c_{2}` are provided by the user.
    
    
    Parameters
    ----------
    constant_contribution : float
        Value of :math:`a_{0}` in above formulation.
    constant_degree_one_contribution : float
        Value of :math:`c_{0}` in above formulation.
    cosine_periodic_degree_one_contribution : float
        Value of :math:`c_{1}` in above formulation.
    sine_periodic_degree_one_contribution : float
        Value of :math:`c_{2}` in above formulation.
    constant_degree_two_contribution : float
        Value of :math:`a_{2}` in above formulation.
    reference_epoch : float
        Reference epoch :math:`t_{0}` of the periodic variation.
    period : float
        Period :math:`T` of the periodic variation.
    Returns
    -------
    SurfacePropertyDistributionSettings
        Object defining settings for surface property distribution"""

def lambertian_body_panel_reflection(reflectivity: float) -> BodyPanelReflectionLawSettings:
    """Function for creating settings for target panel reflection law using a Lambertian model
    
    Function for creating settings for target panel reflection law used for a radiation pressure target, with a
    purely Lambertian model. The implementation is as :func:`specular_diffuse_body_panel_reflection`, with
    :math:`\\rho=0` and no instantaneous reradiation. The only free parameter is the reflectivity :math:`\\delta`, such that
    :math:`\\alpha=1-\\delta`.
    
    
    Parameters
    ----------
    reflectivity : float
        Reflectivity :math:`\\delta`
    Returns
    -------
    BodyPanelReflectionLawSettings
        Object defining settings for target panel reflection law"""

def panelled_extended_radiation_source(panel_radiosity_settings: list[PanelRadiosityModelSettings], number_of_panels_per_ring: list[int], original_source_occulting_bodies: dict[str, list[str]]={}) -> RadiationSourceModelSettings:
    """Function for creating settings for a dynamically panelled extended radiation source
    
    Function for creating settings for a radiation source that is modelled as an anisotropic extended source,
    such as a source due to albedo or planetary IR. The model can combined any number of superimposed surface panel radiosity models
    (e.g. albedo, direct radiation), each of which may or may not involve an 'original source' (e.g. the Sun for albedo).
    Each time the radiation at a given target is computed, the surface of the body is re-panelled, using the algorithm described by
    Knocke (1989). In short, a single panel is placed at the zenith of the evaluation point, with any number of rings around it, each of
    which has any number of (equispaced) panels on it. The width of each ring is such that all panels have the same projected, attenuated area.
    The panelling settings are defined by the user to this function. The
    The irradiance :math:`\\Phi` at a relative position :math:`\\mathbf{r}` from the source's center is then computed from all
    :math:`N` panels, each of which evaluated :math:`M` panel radiosity models
    
    .. math::
       \\Phi=\\sum_{i=1}^{N}\\sum_{j=1}\\Phi_{i,j}
    
    where :math:`\\Phi_{i,j}` denotes the contribution to the total irradiance of panel radiosity model :math:`j` on panel :math:`i`.
    
    
    Parameters
    ----------
    luminosity_model : LuminosityModelSettings
        Settings for the luminosity model.
    Returns
    -------
    RadiationSourceModelSettings
        Object defining settings for source model irradiance"""

def panelled_radiation_target(source_to_target_occulting_bodies: dict[str, list[str]]={}, maximum_number_of_pixels_per_source: dict[str, int]={}) -> RadiationPressureTargetModelSettings:
    """Function for creating settings for a paneled radiation pressure target model
    
    Function for creating settings for a paneled radiation pressure target model. Each source can have
    its own set of occulting bodies.
    This model requires the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.vehicle_shape_settings` of type :class:`~tudatpy.dynamics.environment_setup.vehicle_systems.FullPanelledBodySettings` to be defined.
    The functions to define the panelled body settings are available in the :ref:`vehicle_systems` module.
    
    Parameters
    ----------
    source_to_target_occulting_bodies : dict[str, list[str]]
        Dictionary (source name -> list of occulting body names) of bodies to occult sources as seen from this target.
    maximum_number_of_pixels : Dict[str, int]
        Maximum number of pixels used in the self-shadowing algorithm per source, omitting a the value or setting it to zero equals to not considering self-shadowing for a given source.
    Returns
    -------
    RadiationPressureTargetModelSettings
        Object defining settings for a radiation pressure target model"""

def predefined_knocke_type_surface_property_distribution(predefined_model: KnockeTypeSurfacePropertyDistributionModel) -> SurfacePropertyDistributionSettings:
    """Function for creating radiative surface property distribution settings according to a predefined 'Knocke-type` model.
    
    As :func:`spherical_harmonic_surface_property_distribution`, but with a predefined spherical harmonic distribution.
    
    
    Parameters
    ----------
    predefined_model : KnockeTypeSurfacePropertyDistributionModel
        Identifier for predefined Knocke-type surface property model.
    Returns
    -------
    SurfacePropertyDistributionSettings
        Object defining settings for surface property distribution"""

def predefined_spherical_harmonic_surface_property_distribution(predefined_model: SphericalHarmonicsSurfacePropertyDistributionModel) -> SurfacePropertyDistributionSettings:
    """Function for creating radiative surface property distribution settings according to a predefined spherical harmonic model.
    
    As :func:`spherical_harmonic_surface_property_distribution`, but with a predefined spherical harmonic distribution.
    
    
    Parameters
    ----------
    predefined_model : SphericalHarmonicsSurfacePropertyDistributionModel
        Identifier for predefined spherical harmonic surface property model.
    Returns
    -------
    SurfacePropertyDistributionSettings
        Object defining settings for surface property distribution"""

def specular_diffuse_body_panel_reflection(specular_reflectivity: float, diffuse_reflectivity: float, with_instantaneous_reradiation: bool) -> BodyPanelReflectionLawSettings:
    """Function for creating settings for target panel reflection law using a specular-diffuse model
    
    Function for creating settings for target panel reflection law used for a radiation pressure target, with a
    specular diffuse model. The details of the implementation are given by Montenbruck et al. (2015). The reflection
    law is defined by the absorption coefficient :math:`\\alpha`, diffuse reflectivity :math:`\\delta` and specular reflectivity
    :math:`\\rho`, which must meet the condition :math:`\\alpha+\\delta+\\rho=1`. For the model definition, the user provides
    :math:`\\alpha` and :math:`\\delta` (and :math:`\\rho` is calculated). The reaction vector :math:`\\hat{\\mathbf{f}}` for a panel
    with surface normal :math:`\\hat{\\mathbf{n}}`, and unit vector from panel surface to source :math:`\\hat{\\mathbf{r}}` then becomes:
    
    .. math::
       \\hat{\\mathbf{f}}=\\cos\\theta\\left((\\alpha+\\delta)\\hat{\\mathbf{r}}+(\\frac{2}{3}\\delta+2\\rho\\cos\\theta)\\hat{\\mathbf{n}} \\right)
    
    In addition, it can be specified whether the absorbed radiation is also instantaneously retransmitted (according to Lambert's law), in
    which case the above is modified to become:
    
    .. math::
       \\hat{\\mathbf{f}}=\\cos\\theta\\left((\\alpha+\\delta)\\left(\\hat{\\mathbf{r}}+\\frac{2}{3}\\hat{\\mathbf{n}}\\right)+2\\rho\\cos\\theta\\hat{\\mathbf{n}} \\right)
    
    
    Parameters
    ----------
    specular_reflectivity : float
        Specular reflectivity :math:`\\rho`.
    diffuse_reflectivity : float
        Diffuse reflectivity :math:`\\delta`.
    with_instantaneous_reradiation : bool
        Boolean denoting whether absorbed radiation is instantaneously retransmitted (yes, if true).
    Returns
    -------
    BodyPanelReflectionLawSettings
        Object defining settings for target panel reflection law"""

def spherical_harmonic_surface_property_distribution(cosine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], sine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> SurfacePropertyDistributionSettings:
    """Function for creating radiative surface property distribution settings according to a spherical harmonic model.
    
    Function for creating radiative surface property (e.g. albedo, emissivity, etc.) distribution settings
    according to a spherical harmonic model. The user provides unnormalized cosine and sine coefficients :math:`C_{lm}` and :math:`S_{lm}`,
    from which the surface property :math:`k` is computed from:
    
    .. math::
       k(\\phi,\\theta)=\\sum_{l=0}^{l_{max}}\\sum_{m=0}^{l}{P}_{lm}(\\sin\\phi)\\left({C}_{lm}\\cos m\\theta+{S}_{lm}\\sin m\\theta\\right)
    
    with the angles :math:`\\phi` and :math:`\\theta` the body-fixed latitude and longitude of the evaluation point.
    
    
    Parameters
    ----------
    cosine_coefficients : numpy.ndarray
        Cosine coefficients of surface distribution. Entry (i,j) denotes coefficient :math:`{C}_{ij}` at degree i and order j.
    sine_coefficients : numpy.ndarray
        Sine coefficients of surface distribution. Entry (i,j) denotes coefficient :math:`{C}_{ij}` at degree i and order j.
    Returns
    -------
    SurfacePropertyDistributionSettings
        Object defining settings for surface property distribution"""

def thermal_emission_angle_based_radiosity(minimum_temperature: float, maximum_temperature: float, constant_emissivity: float, original_source_name: str) -> PanelRadiosityModelSettings:
    """Function for creating settings for surface radiosity of an extended source with surface temperature from Lemoine (2013)
    
    Function for creating settings for surface radiosity of an extended source from an isotropically heated body (e.g. IR radiation)
    with surface temperature :math:`T` computed from the angle of the surface normal and the original source as follows:
    
    .. math::
       T=\\max\\left(T_{max}(\\cos\\phi_{s})^{1/4},T_{min} \\right)
    
    with :math:`phi_{s}` the angle along a great circle arc from the panel to the subsolar (for the Sun as original source) point; for
    a circular body equivalent to the angle of the vector to the original source and the surface normal. The minimum and
    maximum temperatures are user parameters.
    
    For a surface panel normal of :math:`\\hat{\\mathbf{n}}`, a vector :math:`\\mathbf{r}` from the surface element to the target,
    the resulting irradiance :math:`\\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):
    
    .. math::
       \\Phi=\\epsilon kT^{4}\\frac{A\\cos\\theta}{\\pi ||\\mathbf{r}||^{2}}
    
    with :math:`\\epsilon` the emissivity, :math:`k` the Stefan-Boltzmann constant, :math:`A` the panel area, :math:`\\theta` is the angle between
    :math:`\\hat{\\mathbf{n}}` and :math:`\\mathbf{r}`.
    
    
    Parameters
    ----------
    minimum_temperature : float
        Minimum surface temperature :math:`T_{min}`.
    maximum_temperature : float
        Maximum surface temperature :math:`T_{min}`.
    constant_emissivity : float
        Constant emissivity of the surface :math:`\\epsilon`.
    original_source_name : str
        Name of the original source from which the radiation is reflection to the target.
    Returns
    -------
    PanelRadiosityModelSettings
        Object defining settings for source panel radiosity"""

def thermal_emission_blackbody_constant_emissivity(constant_emissivity: float, original_source_name: str) -> PanelRadiosityModelSettings:
    """Function for creating settings for surface radiosity of an extended source from an isotropically heated body with constant emissivity
    
    Function for creating settings for surface radiosity of an extended source from an isotropically heated body (e.g. IR radiation) with constant surface
    emissivity,
    where the emitted power of the body is computed from the assumption that all heat absorbed from an original source is
    emitted isotropically by the body. For instance, for Earth with Sun as original source, this model is equivalent to
    assuming that a given fraction of all heat incident of the Sun on the Earth is absorbed and causes the full Earth surface to
    heat to a constant temperature, which then results in the body emitting infrared radiation from its surface.
    
    For a surface panel normal of :math:`\\hat{\\mathbf{n}}`, a vector :math:`\\mathbf{r}` from the surface element to the target,
    the resulting irradiance :math:`\\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):
    
    .. math::
       \\Phi=\\frac{\\epsilon\\Phi_{s}}{4}\\frac{A\\cos\\theta}{\\pi ||\\mathbf{r}||^{2}}
    
    with :math:`\\epsilon` the emissivity, :math:`\\Phi_{s}` the irradiance from the original source,  :math:`A` the panel area, :math:`\\theta` is the angle between
    :math:`\\hat{\\mathbf{n}}` and :math:`\\mathbf{r}`.
    
    
    Parameters
    ----------
    constant_emissivity : float
        Constant emissivity of the surface :math:`\\epsilon`.
    original_source_name : str
        Name of the original source from which the radiation is reflection to the target.
    Returns
    -------
    PanelRadiosityModelSettings
        Object defining settings for source panel radiosity"""

def thermal_emission_blackbody_variable_emissivity(emissivity_distribution_model: SurfacePropertyDistributionSettings, original_source_name: str) -> PanelRadiosityModelSettings:
    """Function for creating settings for surface radiosity of an extended source from an isotropically heated body with variable emissivity
    
    As :func:`thermal_emission_blackbody_constant_emissivity`, but with the surface emissivity :math:`\\epsilon` defined by a surface distribution model.
    
    
    Parameters
    ----------
    emissivity_distribution_model : SurfacePropertyDistributionSettings
        Model for the surface distribution of the emissivity :math:`\\epsilon`.
    original_source_name : str
        Name of the original source from which the radiation is reflection to the target.
    Returns
    -------
    PanelRadiosityModelSettings
        Object defining settings for source panel radiosity"""

def time_variable_luminosity(luminosity_function: typing.Callable[[float], float]) -> LuminosityModelSettings:
    """Function for creating time-variable radiation source luminosity settings.
    
    Function for creating time-variable radiation source luminosity settings, defining the total
    radiated power (in Watts) of a given source as a function of time. With this function, the source
    is assumed to emit radiation isotropically.
    
    
    Parameters
    ----------
    luminosity_function : Callable[[float], float]
        Function returning source luminosity (in Watt) as a function of time
    Returns
    -------
    LuminosityModelSettings
        Object defining settings for source luminosity"""

def variable_albedo_surface_radiosity(albedo_distribution_settings: SurfacePropertyDistributionSettings, original_source_name: str) -> PanelRadiosityModelSettings:
    """Function for creating settings for surface variable albedo surface radiosity of an extended source
    
    As :func:`constant_albedo_surface_radiosity`, but with the surface albedo :math:`a` defined by a surface distribution model.
    
    
    Parameters
    ----------
    albedo_distribution_settings : SurfacePropertyDistributionSettings
        Model for the surface distribution of the albedo :math:`a`.
    original_source_name : str
        Name of the original source from which the radiation is reflection to the target.
    Returns
    -------
    PanelRadiosityModelSettings
        Object defining settings for source panel radiosity"""
albedo_dlam1: SphericalHarmonicsSurfacePropertyDistributionModel
albedo_knocke: KnockeTypeSurfacePropertyDistributionModel
cannonball_radiation_pressure_interface: RadiationPressureType
cannonball_target: RadiationPressureTargetModelType
custom: KnockeTypeSurfacePropertyDistributionModel
emissivity_knocke: KnockeTypeSurfacePropertyDistributionModel
multi_type_target: RadiationPressureTargetModelType
paneled_target: RadiationPressureTargetModelType
undefined_target: RadiationPressureTargetModelType