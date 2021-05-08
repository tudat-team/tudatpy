import pandas as pd
import numpy as np
import os
from .kernel.io import *


def save2txt(solution, filename, directory="./"):
    """

    Parameters
    ----------
    solution
    filename
    directory
    column_names

    Returns
    -------

    """
    if os.path.exists(directory):
        pass
    else:
        os.makedirs(directory)
    df = pd.DataFrame(index=solution.keys(),
                      data=np.vstack(list(solution.values())))
    if len(filename.split('.')) > 1:
        _filename = filename
    else:
        _filename = filename + ".txt"
    df.index.name = "time"
    df.to_csv(os.path.join(directory, filename),header=False,sep='\t')
