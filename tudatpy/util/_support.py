import numpy as np
from ..kernel.math import interpolators
from ..kernel import numerical_simulation
from ..kernel.numerical_simulation import propagation_setup
from ..kernel.numerical_simulation.propagation_setup.dependent_variable import SingleDependentVariableSaveSettings, get_dependent_variable_id
import os
from typing import List, Dict, Union
from tudatpy.util import result2array

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
    redirect_file_path : None or string, optional, default=None
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
    def __init__(self, redirect_file_path:Union[None,str]=None, redirect_out: bool=True, redirect_err: bool=True):
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

def pareto_optimums(points: list, operator:Union[None,List[Union[min,max]]]=None):
    """Compute Pareto optimums from a set of points.

    These points are all individually optimums, meaning that to be better in one dimension, they have to be worse in another one.

    This function also allows to specify wether each dimension (which can be seen as an optimisation objective), should be minimised or maximised.

    The set of points can indicate 2D or 3D coordinates, but can also extend to as many dimensions (and objectives) as one wants.

    Parameters
    ----------
    points : list or numpy.ndarray
        Multi-dimensional list that contains the set of points to compute Pareto optimums from.
        If the points are spread in 3D, this list should have 3 columns, and as many rows as there are points.
    operator : None or list[min or max], optional, default=None
        If None, it will be considered that the optimums are the minimums of each axis (dimension).
        Otherwise, a list of `min` or `max` functions can be passed in this input, to specify whether a point along a given dimension should be minimum or maximum to be considered optimum.

    Examples
    --------
    The following code defines a set of equispaced points in 3 dimensions, with spacing of 0.5, ranging from -0.5 to 1.
    
    Wether each of these points is a Pareto optimum is then computed, taking into account that an optimum is a minimum for x and y, and a maximum for z.

    Then, a 3D plot is made, showing the Pareto optimum points in green, and the other ones in red.
    The plot generated by this example code can be seen in the image below.

    .. image:: _static/pareto_optimums_test.png
       :width: 400
       :align: center

    .. code-block:: python

        # Define coordinates at which we want points
        coordinates = [-0.5, 0, 0.5, 1]
        # Assemble every possible combination of the coordinates to make the set of points
        points = np.asarray([[i, j, k] for i in coordinates for j in coordinates for k in coordinates])
        # Compute which points are Pareto Optimums, when taking the minimum for x and y, and maximum for z
        pareto_optimums = TU.pareto_optimums(points, operator=[min, min, max])
        # Create a matplotlib.pyplot figure
        fig = plt.figure()
        # Add axis with a 3D projection to the figure
        ax = fig.add_subplot(projection='3d')
        # Define a color green for the Pareto Optimums (computed as 'True'). Others are red
        color = ["green" if opti else "red" for opti in pareto_optimums]
        # Plot all points with their color
        ax.scatter(points[:,0], points[:,1], points[:,2], c=color)
        # Add axis labels
        ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
        # Show the plot
        plt.show()
    """
    points = np.asarray(points)
    if operator is None:
        sign = np.ones(points.shape[1])
    else:
        if len(operator) != points.shape[1]:
            raise IndexError("The length of the sign argument does not correspond with the number of points.")
        sign = np.asarray([1 if o == min else -1 for o in operator])
    points = np.asarray([p*sign for p in points])
    pareto_optimal = np.ones(points.shape[0], dtype=bool)
    for i, c in enumerate(points):
        if pareto_optimal[i]:
            pareto_optimal[pareto_optimal] = np.any(points[pareto_optimal]<=c, axis=1)
    return pareto_optimal

def split_history(state_history: Dict[float, np.array], propagator_settings: propagation_setup.propagator.PropagatorSettings):
    """Split the state history into a distinct state histories for each body.

    Creates a dictionnary of state histories based on the unified `state_history`
    from the propagation of multiple bodies. Each dictionnary key contains the name of a propagated body,
    and the value is the state history for the given propagated body.

    Parameters
    -----------
    state_history : Dict[float, numpy.ndarray]
        Dictionary mapping the simulation time steps to the propagated
        state time series.

    propagator_settings : tudatpy.kernel.numerical_simulation.propagation_setup.propagator.PropagatorSettings
        Settings used for the propagation.

    Returns
    -----------
    state_history_book : Dict[str,[Dict[float, numpy.ndarray]]]
        Dictionnary containing the name of the propagated body as key, and the state history as value.
    """    
    # Get the propagated state types and names of the propagated bodies from the integrator settings.
    integrated_type_and_body_list = numerical_simulation.get_integrated_type_and_body_list(propagator_settings)

    # Extract the states and epochs from the state history.
    states_list = np.asarray(list(state_history.values()))
    epochs = list(state_history.keys())
    
    # Loop trough the state types and bodies name to save them beforehand.
    n_bodies, body_names = None, None
    propagated_states_sizes = []
    for state_type, body_list in integrated_type_and_body_list.items():
        if n_bodies is None:
            n_bodies = len(body_list)
            body_names = [body_list[i][0] for i in range(n_bodies)]
        # Get the state size for the current state type.
        state_size = numerical_simulation.get_single_integration_size(state_type)
        propagated_states_sizes.append(state_size)

    # Create the empty state history book.
    state_history_book = {body_name: {epoch: [] for epoch in epochs} for body_name in body_names}

    # Loop through the epochs and states to fill the state history book.
    for epoch, state in zip(epochs, states_list):
        state_idx = 0
        # Loop trough the state types by their propagated vector size.
        for propagated_state_size in propagated_states_sizes:
            # Loop trough the propagated bodies names.
            for body_name in body_names:
                # Add the section of the state related to the current state type and body to the state history book.
                state_history_book[body_name][epoch].extend(state[state_idx:state_idx+propagated_state_size])
                # Update the state index cursor.
                state_idx += propagated_state_size

    # Return the state history book.
    return state_history_book

def vector2matrix(flat_matrix: np.ndarray) :
    """Convert a flattened matrix into a matrix.

    Following Tudat standards, a rotation matrix is returned as a nine-entry vector in the dependent variable output,
    where entry (i,j) of the matrix is stored in entry (3i+j) of the vector with i,j = 0,1,2.
    This is detailled in the :func:`~tudatpy.numerical_simulation.propagation_setup.dependent_variable.inertial_to_body_fixed_rotation_frame` docs.

    Parameters
    -----------
    flat_matrix : numpy.ndarray
        Vector containing a flattened rotation matrix.

    Returns
    -----------
    rotation_matrix: numpy.ndarray
        Rotation matrix (3x3 orthogonal matrix).
    """
    return flat_matrix.reshape(3,3)