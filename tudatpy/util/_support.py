import numpy as np
from ..kernel.math import interpolators
import os
from typing import List, Dict


def result2array(result: Dict[float, np.array]):
    """Initial prototype function to convert dict result from DynamicsSimulator

    The `state_history` and `dependent_history` retrieved from classes
    deriving from the :class:`~tudatpy.numerical_simulation.SingleArcSimulator`
    return these time series as a mapping illustrated by:

    .. code-block:: python

        state_history = {
            t[0]: np.array([pos_x[0], pos_y[0], pos_z[0], vel_x[0], vel_y[0], vel_z[0]]),
            t[1]: np.array([pos_x[1], pos_y[1], pos_z[1], vel_x[1], vel_y[1], vel_z[1]]),
            t[2]: np.array([pos_x[2], pos_y[2], pos_z[2], vel_x[2], vel_y[2], vel_z[2]]),
            # ... clipped
            t[-1]: np.array([pos_x[-1], pos_y[-1], pos_z[-1], vel_x[-1], vel_y[-1], vel_z[-1]]),
        }

    `result2array` is a utility function to convert the above into the
    more conventional :class:`numpy.ndarray` as illustrated by:

    .. code-block:: python 

        states_array = np.array([
            [t[0], pos_x[0], pos_y[0], pos_z[0], vel_x[0], vel_y[0], vel_z[0]],
            [t[1], pos_x[1], pos_y[1], pos_z[1], vel_x[1], vel_y[1], vel_z[1]],
            [t[2], pos_x[2], pos_y[2], pos_z[2], vel_x[2], vel_y[2], vel_z[2]],
            # ... clipped
            [t[-1], pos_x[-1], pos_y[-1], pos_z[-1], vel_x[-1], vel_y[-1], vel_z[-1]],
        ])

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

def compare_results(baseline_results: Dict[float, np.array], new_results: Dict[float, np.array], difference_epochs: List[float]):
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

    Examples
    --------
    .. code-block:: python

        # Create a first orbital simulation around Earth
        simulation_baseline = SingleArcSimulator(..., integrator_settings, ...)
        # Extract the simulation epochs for the baseline
        epochs_baseline = list(simulation_baseline.state_history.keys())

        # Create the same orbital simulation with a lower tolerance for the integrator
        simulation_faster = SingleArcSimulator(..., integrator_settings_lower_tolerance, ...)

        # Setup a list of epochs starting at the beginning of the baseline simulation and spanning 3 hours, with a timestep of 30 seconds
        compare_times = np.arange(epochs_baseline[0], epochs_baseline[0]+3*24*3600, 30)

        # Compute the difference between the two simulations
        simulations_difference = util.compare_results(simulation_baseline, simulation_faster, compare_times)
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

class redirect_std():
    """Redirect any print that is sent by a noisy function by encapsulating it with this class.

    The print will successfully be redirected even if they are sent by a C++ function (or from other language).

    Exceptions that are raised will still show in the terminal as excepted.


    Parameters
    ----------
    redirect_file_path : None or string, default=None
        If None, the prints are redirected to Dev Null, and are thus suppressed.
        Else, redirect_file specifies in what file the messages are saved to.
    redirect_out : Boolean, optional, default=True
        If True, redirect anything that is sent to the terminal trough the STD OUT method.
    redirect_err : Boolean, optional, default=True
        If True, redirect anything that is sent to the terminal trough the STD ERR method.

    Examples
    --------
    The following code will for instance run a single arc simulation with no print to the console.
    Any outputs that would normally be printed on the terminal are save in the file `C:/log/single_arc_log.txt`.

    .. code-block:: python

        with util.redirect_std(redirect_file_path="C:/log/single_arc_log.txt"):
            simulation_baseline = SingleArcSimulator(...)
    """
    # Class adapted from https://stackoverflow.com/q/11130156/11356694
    def __init__(self, redirect_file_path=None, redirect_out=True, redirect_err=True):
        # Create the file in Dev Null to dispose of the unwanted prints
        if redirect_file_path is None:
            self.null_files = os.open(os.devnull,os.O_RDWR)
        # Create the file at the specified path to redirect the print messages
        elif redirect_out or redirect_err:
            f = open(redirect_file_path, "w")
            f.close()
            self.null_files = os.open(redirect_file_path,os.O_RDWR)
        # Save the streams of the real STD OUT and STD ERR
        self.std_streams = [os.dup(1), os.dup(2)]
        # Save what STD method should be muted
        self.redirect_out, self.redirect_err = redirect_out, redirect_err

    def __enter__(self):
        # Link any print trough STD OUT to the Null pointer
        if self.redirect_out:
            os.dup2(self.null_files,1)
        # Link any print trough STD ERR to the Null pointer
        if self.redirect_err:
            os.dup2(self.null_files,2)

    def __exit__(self, *_):
        # Link STD OUT and ERR back to their real stream
        os.dup2(self.std_streams[0],1), os.dup2(self.std_streams[1],2)
        # Close all links
        for link in [self.null_files] + self.std_streams:
            os.close(link)