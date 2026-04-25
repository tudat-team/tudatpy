from pathlib import Path

# In editable installs, this file may be symlinked from source while the package
# directory in site-packages is stale. Extend __path__ with the resolved source
# directory so newly added subpackages are discoverable without requiring a full
# reinstall after every package-structure change.
_current_package_dir = Path(__file__).parent
_resolved_package_dir = Path(__file__).resolve().parent
if _resolved_package_dir != _current_package_dir and str(_resolved_package_dir) not in __path__:
    __path__.append(str(_resolved_package_dir))

from tudatpy.kernel.data import *
from ._support import save2txt, save_time_history_to_file
from .mission_data_downloader import LoadPDS, DownloadAtmosphericData
from .processTrk234 import Trk234Processor
from .ancillary import (
    IonexProduct,
    IonexResolution,
    VmfTechnique,
    VmfProcessing,
    DownloadResult,
    AncillaryDownloadError,
    AuthenticationError,
    download_ionex,
    download_vmf,
    download_ancillary,
)

# This would generate a circular import (SBDB in environment_setup)
# from . import horizons, mpc, sbdb, spacetrack, discos, mission_data_downloader
