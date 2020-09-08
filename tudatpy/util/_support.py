import numpy as np


def result2array(result):
    """Initial prototype function to convert dict result from DynamicsSimulator

    return
    ------

    numpy.ndarray: Array of converted results. First column is time.

    """
    # Convert dict_values into list before stacking them.
    dict_values_list = list(result.values())

    # Stack to form m x 6 array of independent variables.
    independent_array = np.vstack(dict_values_list)

    # Get time list from dict.
    dict_keys_list = result.keys()

    # Convert time list to array and result to m X 1
    time_array = np.array(list(dict_keys_list)).reshape(-1, 1)

    # Stack horizontally and return.
    return np.hstack((time_array, independent_array))
