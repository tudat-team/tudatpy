from tudatpy._deprecation import property_deprecation
from tudatpy.kernel.estimation.observations import *
from tudatpy.estimation.observations._observation_collection_helpers import *

SingleObservationSet.ancilliary_settings = property_deprecation(
    "SingleObservationSet.ancilliary_settings", "SingleObservationSet.ancillary_settings"
)(SingleObservationSet.ancillary_settings)
