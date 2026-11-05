from tudatpy.data.sbdb import _ALIASES, __all__

from tudatpy.data._compat import deprecated_dir, deprecated_getattr


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), __all__)
