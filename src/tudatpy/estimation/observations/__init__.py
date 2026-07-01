from tudatpy._deprecation import property_deprecation
from tudatpy.estimation.observable_models_setup.links import (
    LinkDefinition,
    LinkEndId,
    LinkEndType,
    observed_body,
    observer,
    receiver,
    transmitter,
)
from tudatpy.estimation.observable_models_setup.model_settings import (
    ObservableType,
    angular_position_type as angular_position,
    angular_position_type,
    azimuth_elevation_type as azimuth_elevation,
    azimuth_elevation_type,
    one_way_instantaneous_doppler_type as one_way_doppler,
    one_way_instantaneous_doppler_type,
    one_way_range_type as one_way_range,
    one_way_range_type,
)
from tudatpy.kernel.estimation.observations import *

from ._query import observation_query

for _name, _object in list(globals().items()):
    if getattr(_object, "__module__", None) == "tudatpy.kernel.estimation.observations":
        _object.__module__ = "tudatpy.estimation.observations"

SingleObservationSet.ancilliary_settings = property_deprecation(
    "SingleObservationSet.ancilliary_settings", "SingleObservationSet.ancillary_settings"
)(SingleObservationSet.ancillary_settings)
