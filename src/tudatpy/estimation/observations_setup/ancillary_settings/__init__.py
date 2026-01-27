from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *
import sys
import warnings


def with_deprecation_warning(old_name, new_name, function):

    warnings.warn(f"{old_name} is deprecated. Use {new_name} instead.")

    return function


mod = sys.modules[__name__]
for argument in dir(mod):

    if "ancillary" in argument:

        old_name = argument.replace("ancillary", "ancilliary")
        interface = with_deprecation_warning(
            old_name, argument, getattr(mod, argument)
        )
        mod.__setattr__(old_name, interface)

# dsn_n_way_doppler_ancillary_settings = dsn_n_way_doppler_ancilliary_settings
# dsn_n_way_range_ancillary_settings = dsn_n_way_range_ancilliary_settings
