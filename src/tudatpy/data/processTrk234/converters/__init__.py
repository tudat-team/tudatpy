from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "Converter": "tudatpy.data_input.tracking_data.tnf._converters.Converter",
    "RadioBase": "tudatpy.data.processTrk234._legacy_converters.RadioBase",
    "RampConverter": "tudatpy.data_input.tracking_data.tnf._converters.RampConverter",
    "OpenRampHandling": "tudatpy.data_input.tracking_data.tnf._converters.OpenRampHandling",
    "DerivedDopplerConverter": "tudatpy.data.processTrk234._legacy_converters.DerivedDopplerConverter",
    "DerivedSraRangeConverter": "tudatpy.data.processTrk234._legacy_converters.DerivedSraRangeConverter",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
