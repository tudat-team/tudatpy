import typing

def save2txt(solution, filename, directory: str='./') -> None:
    """Save a vector or matrix history to a file.
    
    This function can be used to save a dictionary that maps epochs to a vector or matrix at the given epoch.
    
    Parameters
    ----------
    solution: Dict[float, numpy.ndarray]
        Dictionary mapping floats (e.g. the simulation time steps) to arrays (e.g. the propagated
        state time series).
    filename: str
        Name of the text file that is to be saved.
    directory: str, optional, default="./"
        Directory in which to save the text file."""

def save_time_history_to_file(solution, filename, directory: str='./') -> None:
    """Save a propagated time history to a file.
    
    This function can be used to save a propagated state history to a text file.
    It can also be used for instance to save a dependent variable history, or a sensitivity matrix history.
    
    .. note::
        This function is essentially the same method as :func:`save2txt`, offering the same functionality under a different name.
    
    Parameters
    ----------
    solution: Dict[float, numpy.ndarray]
        Dictionary mapping the simulation time steps to the propagated
        state time series.
    filename: str
        Name of the text file that is to be saved.
    directory: str, optional, default="./"
        Directory in which to save the text file."""