from tudatpy._deprecation import deprecated_property_alias
from tudatpy.kernel.estimation.observations import *

SingleObservationSet.ancilliary_settings = deprecated_property_alias("SingleObservationSet.ancilliary_settings", "SingleObservationSet.ancillary_settings")(SingleObservationSet.ancillary_settings)
