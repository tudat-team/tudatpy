"""Neutral atmospheric S1/S2 loading records and strict SI coefficient readers.

Records preserve the provider axes, both longitude seam endpoints, original
units and source metadata. Reading does not choose a geographic convention,
interpolate, synthesize displacement or change the source reference origin.
"""

from tudatpy.kernel.data_input.environment_data.atmospheric_loading import *
