"""MPC 80-column optical astrometry loading."""

from .parsers import parse_80cols_file


def read_80_column_data(*args, **kwargs):
    from tudatpy.data_access.tracking.mpc import read_80_column_data as _read_80_column_data

    return _read_80_column_data(*args, **kwargs)


__all__ = ["parse_80cols_file", "read_80_column_data"]
