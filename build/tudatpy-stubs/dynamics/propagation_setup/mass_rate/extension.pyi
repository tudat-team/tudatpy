import typing
__all__ = ['AvailableMassRateModels', 'CustomMassRateSettings', 'FromThrustMassRateSettings', 'MassRateModelSettings', 'custom', 'custom_mass_rate', 'custom_mass_rate_type', 'from_thrust', 'from_thrust_mass_rate_type', 'undefined_mass_rate_type']

class AvailableMassRateModels:
    """Enumeration of available mass rate models.
    
             Enumeration of mass rate models supported by tudat.
    
    
    
    
    
          
    
    Members:
    
      undefined_mass_rate_type : 
          
    
      custom_mass_rate_type : 
          
    
      from_thrust_mass_rate_type : 
          """
    __members__: typing.ClassVar[dict[str, AvailableMassRateModels]]
    custom_mass_rate_type: typing.ClassVar[AvailableMassRateModels]
    from_thrust_mass_rate_type: typing.ClassVar[AvailableMassRateModels]
    undefined_mass_rate_type: typing.ClassVar[AvailableMassRateModels]

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

class CustomMassRateSettings(MassRateModelSettings):
    """`MassRateModelSettings`-derived class to define settings for a custom mass rate model.
    
    `MassRateModelSettings`-derived class to define settings for a custom mass rate model."""

class FromThrustMassRateSettings(MassRateModelSettings):
    """`MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model.
    
    `MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model."""

class MassRateModelSettings:
    """Functional base class to define settings for mass rates.
    
    Base class for providing settings for a mass rate model, that defines the models to be used to numerically propagate the
    mass of a body during a simulation. If any additional information (besides the type of the mass rate model) is required,
    these must be implemented in a derived class (dedicated for each mass rate model type)."""

def custom(mass_rate_function: typing.Callable[[float], float]) -> MassRateModelSettings:
    ...

def custom_mass_rate(mass_rate_function: typing.Callable[[float], float]) -> MassRateModelSettings:
    """Creates the settings for a mass rate model defined from a thrust model.
    
    Creates the settings for a custom mass rate model defined through a mass rate function. The function must have
    time as an independent variable.
    
    
    Parameters
    ----------
    mass_rate_function : callable[[float], float]
        Function of time defining the custom mass rate.
    Returns
    -------
    CustomMassRateSettings
        Custom mass rate settings object."""

def from_thrust(use_all_thrust_models: bool=1, associated_thrust_source: str='') -> MassRateModelSettings:
    """Creates the settings for a mass rate model defined from a thrust model.
    
    Creates the settings for a mass rate model defined from a thrust model. The mass rate model is derived from
    the associated body's engine model. It is possible to consider only a specific engine or all engines.
    
    
    Parameters
    ----------
    use_all_thrust_models : bool, default=true
        Denotes whether all engines of the associated body are to be combined into a single thrust model.
    associated_thrust_source : str, default=""
        Name of engine model from which thrust is to be derived (must be empty if the first argument is set to true).
    Returns
    -------
    FromThrustMassRateSettings
        From thrust mass rate settings object."""
custom_mass_rate_type: AvailableMassRateModels
from_thrust_mass_rate_type: AvailableMassRateModels
undefined_mass_rate_type: AvailableMassRateModels