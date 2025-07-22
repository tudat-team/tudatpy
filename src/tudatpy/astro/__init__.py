from tudatpy.kernel.astro import *

import sys
import importlib
def __getattr__(name):
    if name == "time_conversion":
        mod = importlib.import_module("tudatpy.astro.time_conversion")
        sys.modules["tudatpy.astro.time_conversion"] = mod
        return mod
    raise AttributeError(f"module {__name__} has no attribute {name}")