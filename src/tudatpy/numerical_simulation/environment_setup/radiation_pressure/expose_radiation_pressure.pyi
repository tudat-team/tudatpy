import numpy
import typing
__all__ = ['BodyPanelReflectionLawSettings', 'CannonBallRadiationPressureInterfaceSettings', 'KnockeTypeSurfacePropertyDistributionModel', 'LuminosityModelSettings', 'PanelRadiosityModelSettings', 'RadiationPressureInterfaceSettings', 'RadiationPressureTargetModelSettings', 'RadiationPressureTargetModelType', 'RadiationPressureType', 'RadiationSourceModelSettings', 'SphericalHarmonicsSurfacePropertyDistributionModel', 'SurfacePropertyDistributionSettings', 'albedo_dlam1', 'albedo_knocke', 'cannonball', 'cannonball_radiation_pressure_interface', 'cannonball_radiation_target', 'cannonball_target', 'constant_albedo_surface_radiosity', 'constant_luminosity', 'constant_radiosity', 'constant_surface_property_distribution', 'custom', 'custom_surface_property_distribution', 'emissivity_knocke', 'irradiance_based_constant_luminosity', 'irradiance_based_time_variable_luminosity', 'isotropic_radiation_source', 'knocke_type_surface_property_distribution', 'lambertian_body_panel_reflection', 'multi_type_target', 'paneled_target', 'panelled_extended_radiation_source', 'panelled_radiation_target', 'predefined_knocke_type_surface_property_distribution', 'predefined_spherical_harmonic_surface_property_distribution', 'specular_diffuse_body_panel_reflection', 'spherical_harmonic_surface_property_distribution', 'thermal_emission_angle_based_radiosity', 'thermal_emission_blackbody_constant_emissivity', 'thermal_emission_blackbody_variable_emissivity', 'time_variable_luminosity', 'undefined_target', 'variable_albedo_surface_radiosity']

class BodyPanelReflectionLawSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class CannonBallRadiationPressureInterfaceSettings(RadiationPressureInterfaceSettings):
    """Class for defining model settings of a cannonball radiation pressure interface.
	
	`RadiationPressureInterfaceSettings` derived class for cannonball radiation pressure interface model settings.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class KnockeTypeSurfacePropertyDistributionModel:
    """Members:
	
	custom
	
	albedo_knocke
	
	emissivity_knocke
	"""
    __members__: typing.ClassVar[dict[str, KnockeTypeSurfacePropertyDistributionModel]]
    albedo_knocke: typing.ClassVar[KnockeTypeSurfacePropertyDistributionModel]
    custom: typing.ClassVar[KnockeTypeSurfacePropertyDistributionModel]
    emissivity_knocke: typing.ClassVar[KnockeTypeSurfacePropertyDistributionModel]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class PanelRadiosityModelSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RadiationPressureInterfaceSettings:
    """Base class for providing settings for radiation pressure interface models.
	
	Functional (base) class for settings of radiation pressure interface models that require no information in addition to their type.
	Radiation pressure interface model settings requiring additional information must be defined using an object derived from this class.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RadiationPressureTargetModelSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RadiationPressureTargetModelType:
    """Members:
	
	cannonball_target :
	
	paneled_target :
	
	multi_type_target :
	
	undefined_target :
	"""
    __members__: typing.ClassVar[dict[str, RadiationPressureTargetModelType]]
    cannonball_target: typing.ClassVar[RadiationPressureTargetModelType]
    multi_type_target: typing.ClassVar[RadiationPressureTargetModelType]
    paneled_target: typing.ClassVar[RadiationPressureTargetModelType]
    undefined_target: typing.ClassVar[RadiationPressureTargetModelType]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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
    """Enumeration of available radiation pressure types.
	
	
	:member cannonball_radiation_pressure_interface:
	:member panelled_radiation_pressure_interface:
	:member solar_sailing_radiation_pressure_interface:
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, RadiationPressureType]]
    cannonball_radiation_pressure_interface: typing.ClassVar[RadiationPressureType]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class SphericalHarmonicsSurfacePropertyDistributionModel:
    """Members:
	
	albedo_dlam1
	"""
    __members__: typing.ClassVar[dict[str, SphericalHarmonicsSurfacePropertyDistributionModel]]
    albedo_dlam1: typing.ClassVar[SphericalHarmonicsSurfacePropertyDistributionModel]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def cannonball(source_body: str, reference_area: float, radiation_pressure_coefficient: float, occulting_bodies: list[str]=[]) -> RadiationPressureInterfaceSettings:
    """Factory function for creating cannonball radiation pressure interface model settings.
	
	Factory function for settings object, defining a cannonball radiation pressure interface model,
	In this model the effective force is co-linear with the vector from radiation source to the body experiencing the force.
	
	
	:param source_body:
			Name of body emitting the radiation.
	:param reference_area:
			Surface area that undergoes radiation pressure.
	:param radiation_pressure_coefficient:
			Radiation pressure coefficient.
	:param occulting_bodies:
			List of bodies causing (partial) occultation.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.RadiationPressureInterfaceSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.CannonBallRadiationPressureInterfaceSettings` class
	"""

def cannonball_radiation_target(reference_area: float, radiation_pressure_coefficient: float, per_source_occulting_bodies: dict[str, list[str]]={}) -> RadiationPressureTargetModelSettings:
    ...

def constant_albedo_surface_radiosity(constant_albedo: float, original_source_name: str) -> PanelRadiosityModelSettings:
    ...

def constant_luminosity(luminosity: float) -> LuminosityModelSettings:
    ...

def constant_radiosity(radiosity: float) -> PanelRadiosityModelSettings:
    ...

def constant_surface_property_distribution(constant_value: float) -> SurfacePropertyDistributionSettings:
    ...

def custom_surface_property_distribution(custom_function: typing.Callable[[float, float, float], float]) -> SurfacePropertyDistributionSettings:
    ...

def irradiance_based_constant_luminosity(constant_irradiance: float, reference_distance: float) -> LuminosityModelSettings:
    ...

def irradiance_based_time_variable_luminosity(irradiance_function: typing.Callable[[float], float], reference_distance: float) -> LuminosityModelSettings:
    ...

def isotropic_radiation_source(luminosity_model: LuminosityModelSettings) -> RadiationSourceModelSettings:
    ...

def knocke_type_surface_property_distribution(constant_contribution: float, constant_degree_one_contribution: float, cosine_periodic_degree_one_contribution: float, sine_periodic_degree_one_contribution: float, constant_degree_two_contribution: float, reference_epoch: float, period: float) -> ...:
    ...

def lambertian_body_panel_reflection(reflectivity: float) -> BodyPanelReflectionLawSettings:
    ...

def panelled_extended_radiation_source(panel_radiosity_settings: list[PanelRadiosityModelSettings], number_of_panels_per_ring: list[int], original_source_occulting_bodies: dict[str, list[str]]={}) -> RadiationSourceModelSettings:
    ...

def panelled_radiation_target(source_to_target_occulting_bodies: dict[str, list[str]]={}) -> RadiationPressureTargetModelSettings:
    ...

def predefined_knocke_type_surface_property_distribution(predefined_model: KnockeTypeSurfacePropertyDistributionModel) -> SurfacePropertyDistributionSettings:
    ...

def predefined_spherical_harmonic_surface_property_distribution(predefined_model: SphericalHarmonicsSurfacePropertyDistributionModel) -> SurfacePropertyDistributionSettings:
    ...

def specular_diffuse_body_panel_reflection(specular_reflectivity: float, diffuse_reflectivity: float, with_instantaneous_reradiation: bool) -> BodyPanelReflectionLawSettings:
    ...

def spherical_harmonic_surface_property_distribution(cosine_coefficients: numpy.ndarray, sine_coefficients: numpy.ndarray) -> SurfacePropertyDistributionSettings:
    ...

def thermal_emission_angle_based_radiosity(minimum_temperature: float, maximum_temperature: float, constant_emissivity: float, original_source_name: str) -> PanelRadiosityModelSettings:
    ...

def thermal_emission_blackbody_constant_emissivity(constant_emissivity: float, original_source_name: str) -> PanelRadiosityModelSettings:
    ...

def thermal_emission_blackbody_variable_emissivity(emissivity_distribution_model: SurfacePropertyDistributionSettings, original_source_name: str) -> PanelRadiosityModelSettings:
    ...

def time_variable_luminosity(luminosity_function: typing.Callable[[float], float]) -> LuminosityModelSettings:
    ...

def variable_albedo_surface_radiosity(albedo_distribution_settings: SurfacePropertyDistributionSettings, original_source_name: str) -> PanelRadiosityModelSettings:
    ...
albedo_dlam1: SphericalHarmonicsSurfacePropertyDistributionModel
albedo_knocke: KnockeTypeSurfacePropertyDistributionModel
cannonball_radiation_pressure_interface: RadiationPressureType
cannonball_target: RadiationPressureTargetModelType
custom: KnockeTypeSurfacePropertyDistributionModel
emissivity_knocke: KnockeTypeSurfacePropertyDistributionModel
multi_type_target: RadiationPressureTargetModelType
paneled_target: RadiationPressureTargetModelType
undefined_target: RadiationPressureTargetModelType