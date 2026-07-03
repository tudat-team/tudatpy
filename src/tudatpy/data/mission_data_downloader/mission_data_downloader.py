#!/usr/bin/env python
# coding: utf-8

"""Backwards-compatible re-export shim.

The implementation used to live entirely in this file; it has since been
split into ``_patterns``, ``_download``, ``_parsing``, ``_meta_kernel`` and
the per-mission modules under ``missions/``, composed together in ``_core``.
This module keeps the historical import path
(``tudatpy.data.mission_data_downloader.mission_data_downloader``) working.

Atmospheric (IONEX/VMF) downloads used to live here as ``DownloadAtmosphericData``;
that has been superseded by ``tudatpy.data.ancillary`` and removed.
"""

from ._core import LoadPDS

__all__ = ["LoadPDS"]
