import warnings

warnings.warn(
    "`tudatpy.io`is deprecated since version 0.7 and will be removed two minor "
    "versions hence: Import output utilities from `tudatpy.util` and input-data "
    "helpers from `tudatpy.data_input` instead.",
    DeprecationWarning,
    2,
)
