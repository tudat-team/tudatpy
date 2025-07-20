from .extension import element_conversion, frame_conversion, fundamentals, gravitation, polyhedron_utilities, time_representation, two_body_dynamics
import sys
import sys
import importlib
import importlib

def __getattr__(name):
    if name == 'time_conversion':
        mod = importlib.import_module('tudatpy.astro.time_conversion')
        sys.modules['tudatpy.astro.time_conversion'] = mod
        return mod
    raise AttributeError(f'module {__name__} has no attribute {name}')
__all__ = ['element_conversion', 'frame_conversion', 'fundamentals', 'gravitation', 'polyhedron_utilities', 'time_representation', 'two_body_dynamics']