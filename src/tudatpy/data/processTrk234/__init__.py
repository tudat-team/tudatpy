from tudatpy.data._compat import deprecated_dir, deprecated_getattr
from .processor import Trk234Processor

_ALIASES = {
    "OpenRampHandling": "tudatpy.data_input.tracking_data.tnf.OpenRampHandling",
}

__all__ = sorted(set(_ALIASES) | {"Trk234Processor"})


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
