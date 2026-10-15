"""
This module contains functions and classes to read tables of transformation parameters between terrestrial reference
frame realizations (for instance between ITRF2020 and ITRF2014), as published by IGN.

The parameters are read without evaluating the transformation. They are used by the
:class:`~tudatpy.astro.frame_conversion.TerrestrialFrameTransformationModel` class and the
:func:`~tudatpy.astro.frame_conversion.evaluate_linearized_helmert` function.
"""

from tudatpy.kernel.data_input.environment_data.terrestrial_frames import *
