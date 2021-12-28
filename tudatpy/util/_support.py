import numpy as np
from ..kernel.math import interpolators


def result2array(result):
    """Initial prototype function to convert dict result from DynamicsSimulator

    The `state_history` and `dependent_history` retrieved from classes
    deriving from the :class:`~tudatpy.numerical_simulation.SingleArcSimulator`
    return these time series as a mapping illustrated by:

    >>> state_history = {
    ...  t[0]: np.array([pos_x[0], pos_y[0], pos_z[0], vel_x[0], vel_y[0], vel_z[0]]),
    ...  t[1]: np.array([pos_x[1], pos_y[1], pos_z[1], vel_x[1], vel_y[1], vel_z[1]]),
    ...  t[2]: np.array([pos_x[2], pos_y[2], pos_z[2], vel_x[2], vel_y[2], vel_z[2]]),
    ...  # ... clipped
    ...  t[-1]: np.array([pos_x[-1], pos_y[-1], pos_z[-1], vel_x[-1], vel_y[-1], vel_z[-1]]),
    ...  }

    `result2array` is a utility function to convert the above into the
    more conventional :class:`numpy.ndarray` as illustrated by:

    >>> states_array = np.array([
    ...  [t[0], pos_x[0], pos_y[0], pos_z[0], vel_x[0], vel_y[0], vel_z[0]],
    ...  [t[1], pos_x[1], pos_y[1], pos_z[1], vel_x[1], vel_y[1], vel_z[1]],
    ...  [t[2], pos_x[2], pos_y[2], pos_z[2], vel_x[2], vel_y[2], vel_z[2]],
    ...  # ... clipped
    ...  [t[-1], pos_x[-1], pos_y[-1], pos_z[-1], vel_x[-1], vel_y[-1], vel_z[-1]],
    ...  ])

    Parameters
    ----------
    result : Dict[float, numpy.ndarray]
        Dictionary mapping the simulation time steps to the propagated
        state time series.

    Returns
    -------
    array : numpy.ndarray
        Array of converted results. First column is time.

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

def compare_results(baseline_results, new_results, difference_epochs):
    """Compare the results of a baseline simulation with the results of a new different simulation.

    This uses a 8th-order Lagrange interpolator to compute the difference in state of the two simulations at specified epochs.
    Alternatively, the dependent variables can be input instead of the states, to compute the difference in dependent variables at these epochs.

    Parameters
    ----------
    baseline_results : Dict[float, numpy.ndarray]
        Dictionary mapping the simulation time steps to the propagated state time series of the baseline simulation.
    new_results : Dict[float, numpy.ndarray]
        Dictionary mapping the simulation time steps to the propagated state time series of a new simulation different from the baseline.
    difference_epochs : numpy.array or list
        Array containing the epochs at which to compute the difference between the results from the distinct simulations.

    Returns
    -------
    results_comparison :  Dict[float, numpy.ndarray]
        Dictionary with difference_epochs as keys, and values corresponding to the difference between the baseline results and the new results.
    """
    # Setup an 8th-order Lagrange interpolator
    interpolator_settings = interpolators.lagrange_interpolation(8, boundary_interpolation=interpolators.use_boundary_value)

    # Setup the interpolator for the baseline and the new simulations
    baseline_results_interpolator = interpolators.create_one_dimensional_vector_interpolator(baseline_results, interpolator_settings)
    new_results_interpolator = interpolators.create_one_dimensional_vector_interpolator(new_results, interpolator_settings)

    # Compute the different between the baseline and the new results
    results_comparison = {
        epoch: new_results_interpolator.interpolate(epoch) - baseline_results_interpolator.interpolate(epoch)
        for epoch in difference_epochs }

    # Return the difference between the results
    return results_comparison