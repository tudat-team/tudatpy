from ._io import *
from ._path import *
from ._format import *


def recurse_dict(dict, keys):
    """Retrieve nested dict value with list of keys.

    Parameters
    ----------
    dict : Dict[Dict[...]]
        Recursively nested dictionary structure.
    keys : List[str]


    Returns
    -------

    """
    _dict = dict
    for key in keys:
        _dict = _dict[key]
    return _dict
