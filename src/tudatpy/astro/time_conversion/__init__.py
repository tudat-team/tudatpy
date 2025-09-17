import warnings
warnings.warn(
    "tudatpy.astro.time_conversion is deprecated as of v1.0 (see https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html). Use tudatpy.astro.time_representation instead.",
    FutureWarning,
    stacklevel=1
)


from tudatpy.kernel.astro.time_representation import *
