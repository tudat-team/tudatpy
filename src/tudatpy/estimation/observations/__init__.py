from tudatpy.kernel.estimation.observations import *
import warnings


def with_deprecation_warning(old_name, new_name, function):

    warnings.warn(
        (f"Function {old_name} is deprecated" f". Use {new_name} instead"),
    )

    return function


SingleObservationSet.ancilliary_settings = with_deprecation_warning(
    "SingleObservationSet.ancilliary_settings",
    "SingleObservationSet.ancillary_settings",
    SingleObservationSet.ancillary_settings,
)
