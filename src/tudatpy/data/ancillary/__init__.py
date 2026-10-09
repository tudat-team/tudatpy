from tudatpy.data._compat import deprecated_dir, deprecated_getattr

_ALIASES = {
    "IonexProduct": "tudatpy.data_input.data_retrieval.media_corrections.IonexProduct",
    "IonexResolution": "tudatpy.data_input.data_retrieval.media_corrections.IonexResolution",
    "VmfTechnique": "tudatpy.data_input.data_retrieval.media_corrections.VmfTechnique",
    "VmfProcessing": "tudatpy.data_input.data_retrieval.media_corrections.VmfProcessing",
    "DownloadResult": "tudatpy.data_input.data_retrieval.media_corrections.DownloadResult",
    "AncillaryDownloadError": "tudatpy.data_input.data_retrieval.media_corrections.AncillaryDownloadError",
    "AuthenticationError": "tudatpy.data_input.data_retrieval.media_corrections.AuthenticationError",
    "download_ionex": "tudatpy.data_input.data_retrieval.media_corrections.download_ionex",
    "download_vmf": "tudatpy.data_input.data_retrieval.media_corrections.download_vmf",
    "download_ancillary": "tudatpy.data_input.data_retrieval.media_corrections.download_ancillary",
}

__all__ = sorted(_ALIASES)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
