import warnings
warnings.warn(
    "tudatpy.astro.time_conversion is deprecated. Use tudatpy.astro.time_representation instead.",
    FutureWarning,
    stacklevel=1
)


from tudatpy.kernel.astro.time_representation import *
