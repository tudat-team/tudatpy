from . import (
    element_conversion,
    ephemerides,
    frame_conversion,
    fundamentals,
    gravitation,
    polyhedron_utilities,
    time_representation,
    two_body_dynamics,
)


def __getattr__(name):

    if name == "time_conversion":
        import sys
        import importlib

        mod = importlib.import_module("tudatpy.astro.time_conversion")
        sys.modules["tudatpy.astro.time_conversion"] = mod
        return mod
    raise AttributeError(f"module {__name__} has no attribute {name}")
