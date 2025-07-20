import numpy
import pybind11_stubgen.typing_ext
from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
from ....estimation.observations_setup import ancillary_settings
from ....estimation.observations_setup import observations_dependent_variables
import typing
__all__ = ['ObservationCollectionParser', 'ObservationFilterBase', 'ObservationFilterType', 'ObservationParserType', 'ObservationSetSplitterBase', 'ObservationSetSplitterType', 'absolute_value_filtering', 'ancillary_settings_parser', 'dependent_variable_filtering', 'empty_parser', 'epochs_filtering', 'link_end_id_parser', 'link_end_str_parser', 'link_end_type_parser', 'link_ends_parser', 'multi_type_parser', 'nb_observations_splitter', 'observable_type_parser', 'observation_filter', 'observation_parser', 'observation_set_splitter', 'residual_filtering', 'single_link_end_parser', 'time_bounds_filtering', 'time_bounds_parser', 'time_interval_splitter', 'time_span_splitter', 'time_tags_splitter']

class ObservationCollectionParser:
    """No documentation found."""

class ObservationFilterBase:
    """No documentation found."""

class ObservationFilterType:
    """No documentation found.
    
    Members:
    
      residual_filtering
    
      absolute_value_filtering
    
      epochs_filtering
    
      time_bounds_filtering
    
      dependent_variable_filtering"""
    __members__: typing.ClassVar[dict[str, ObservationFilterType]]
    absolute_value_filtering: typing.ClassVar[ObservationFilterType]
    dependent_variable_filtering: typing.ClassVar[ObservationFilterType]
    epochs_filtering: typing.ClassVar[ObservationFilterType]
    residual_filtering: typing.ClassVar[ObservationFilterType]
    time_bounds_filtering: typing.ClassVar[ObservationFilterType]

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

class ObservationParserType:
    """No documentation found.
    
    Members:
    
      empty_parser
    
      observable_type_parser
    
      link_ends_parser
    
      link_end_str_parser
    
      link_end_id_parser
    
      link_end_type_parser
    
      single_link_end_parser
    
      time_bounds_parser
    
      ancillary_settings_parser
    
      multi_type_parser"""
    __members__: typing.ClassVar[dict[str, ObservationParserType]]
    ancillary_settings_parser: typing.ClassVar[ObservationParserType]
    empty_parser: typing.ClassVar[ObservationParserType]
    link_end_id_parser: typing.ClassVar[ObservationParserType]
    link_end_str_parser: typing.ClassVar[ObservationParserType]
    link_end_type_parser: typing.ClassVar[ObservationParserType]
    link_ends_parser: typing.ClassVar[ObservationParserType]
    multi_type_parser: typing.ClassVar[ObservationParserType]
    observable_type_parser: typing.ClassVar[ObservationParserType]
    single_link_end_parser: typing.ClassVar[ObservationParserType]
    time_bounds_parser: typing.ClassVar[ObservationParserType]

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

class ObservationSetSplitterBase:
    """No documentation found."""

class ObservationSetSplitterType:
    """No documentation found.
    
    Members:
    
      time_tags_splitter
    
      time_interval_splitter
    
      time_span_splitter
    
      nb_observations_splitter"""
    __members__: typing.ClassVar[dict[str, ObservationSetSplitterType]]
    nb_observations_splitter: typing.ClassVar[ObservationSetSplitterType]
    time_interval_splitter: typing.ClassVar[ObservationSetSplitterType]
    time_span_splitter: typing.ClassVar[ObservationSetSplitterType]
    time_tags_splitter: typing.ClassVar[ObservationSetSplitterType]

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
def observation_filter(filter_type: ObservationFilterType, filter_value: float, filter_out: bool=True, use_opposite_condition: bool=False) -> ObservationFilterBase:
    """No documentation found."""

@typing.overload
def observation_filter(filter_type: ObservationFilterType, filter_value: list[float], filter_out: bool=True, use_opposite_condition: bool=False) -> ObservationFilterBase:
    """No documentation found."""

@typing.overload
def observation_filter(filter_type: ObservationFilterType, first_filter_value: float, second_filter_value: float, filter_out: bool=True, use_opposite_condition: bool=False) -> ObservationFilterBase:
    """No documentation found."""

@typing.overload
def observation_filter(filter_type: ObservationFilterType, filter_value: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], filter_out: bool=True, use_opposite_condition: bool=False) -> ObservationFilterBase:
    """No documentation found."""

@typing.overload
def observation_filter(dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, filter_value: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], filter_out: bool=True, use_opposite_condition: bool=False) -> ObservationFilterBase:
    """No documentation found."""

@typing.overload
def observation_parser() -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(observable_type: model_settings.ObservableType, use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(observable_type_vector: list[model_settings.ObservableType], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_ends: dict[links.LinkEndType, links.LinkEndId], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_ends_vector: list[dict[links.LinkEndType, links.LinkEndId]], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_ends_str: str, is_reference_point: bool=False, use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_ends_str_vector: list[str], is_reference_point: bool=False, use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_end_id: tuple[str, str], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_end_ids_vector: list[tuple[str, str]], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_end_type: links.LinkEndType, use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(link_end_types_vector: list[links.LinkEndType], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(single_link_end: tuple[links.LinkEndType, links.LinkEndId], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(single_link_ends_vector: list[tuple[links.LinkEndType, links.LinkEndId]], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(time_bounds: tuple[float, float], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(time_bounds_vector: list[tuple[float, float]], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(ancillary_settings: ancillary_settings.ObservationAncilliarySimulationSettings, use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(ancillary_settings_vector: list[ancillary_settings.ObservationAncilliarySimulationSettings], use_opposite_condition: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_parser(observation_parsers: list[ObservationCollectionParser], combine_conditions: bool=False) -> ObservationCollectionParser:
    """No documentation found."""

@typing.overload
def observation_set_splitter(splitter_type: ObservationSetSplitterType, splitter_value: list[float], min_number_observations: int=0) -> ObservationSetSplitterBase:
    """No documentation found."""

@typing.overload
def observation_set_splitter(splitter_type: ObservationSetSplitterType, splitter_value: float, min_number_observations: int=0) -> ObservationSetSplitterBase:
    """No documentation found."""

@typing.overload
def observation_set_splitter(splitter_type: ObservationSetSplitterType, splitter_value: int, min_number_observations: int=0) -> ObservationSetSplitterBase:
    """No documentation found."""
absolute_value_filtering: ObservationFilterType
ancillary_settings_parser: ObservationParserType
dependent_variable_filtering: ObservationFilterType
empty_parser: ObservationParserType
epochs_filtering: ObservationFilterType
link_end_id_parser: ObservationParserType
link_end_str_parser: ObservationParserType
link_end_type_parser: ObservationParserType
link_ends_parser: ObservationParserType
multi_type_parser: ObservationParserType
nb_observations_splitter: ObservationSetSplitterType
observable_type_parser: ObservationParserType
residual_filtering: ObservationFilterType
single_link_end_parser: ObservationParserType
time_bounds_filtering: ObservationFilterType
time_bounds_parser: ObservationParserType
time_interval_splitter: ObservationSetSplitterType
time_span_splitter: ObservationSetSplitterType
time_tags_splitter: ObservationSetSplitterType