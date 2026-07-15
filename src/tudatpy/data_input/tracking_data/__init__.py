"""Tracking-data readers and containers.

This module exposes Tudat tracking-data containers together with readers and
helpers for radiometric, optical, SLR, TNF/TRK-2-34, ODF, IFMS, Fdets, PSF, and
generic tracking text-file inputs.
"""

from tudatpy.kernel.data_input.tracking_data import *

import importlib as _importlib

from . import (
    fdets,
    generic_text_file,
    ifms,
    mpc,
    odf,
    optical_utilities,
    psf,
    slr,
    tnf,
)

globals()["80_column"] = _importlib.import_module(__name__ + ".80_column")
