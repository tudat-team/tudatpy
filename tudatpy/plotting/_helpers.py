import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Union
from ..util import result2array, pareto_optimums

def dual_y_axis(
    x_data: list,
    y_data_1: list,
    y_data_2: list, 
    x_label: str, 
    y_label_1: str,
    y_label_2: str,
    title:str = "", 
    c1:str = "tab:blue", 
    c2:str = "tab:red", 
    grid:str = ""   ):
    """Plot two y-axis that share a common x-axis.

    Parameters
    ----------
    x_data : list
        List of common values that constitute the x-axis.
    y_data_1 : list
        List of values that constitute the first (left) y-axis.
    y_data_2 : list
        List of values that constitute the second (right) y-axis.
    x_label : str
        Label of the common x-axis.
    y_label_1 : str
        Label of the first (left) y-axis.
    y_label_1 : str
        Label of the second (right) y-axis.
    title : str, optional
        Title of the plot.
    c1 : str, optional, default = "tab:blue"
        Color of the first (left) y-axis.
    c2 : str, optional, default = "tab:red"
        Color of the second (right) y-axis.
    grid : str, optional
        Indication on where to use a grid. This variable can take the value "left", "right", or "both"
        to indicate relative to which y-axis the grid should be made.        
        Using any of these options will put a grid relative to the common x-axis.
    

    Examples
    --------
    .. code-block:: python

        # Define data to be plotted
        x = np.arange(0, 20, 0.01)
        y1 = x**5*np.sin(x)
        y2 = 5**x*(np.cos(x)+1)

        # Plot the results from the two functions of x
        fig, ax1, ax2 = dual_y_axis(
            x,
            y1,
            y2,
            "x",
            "$x^5 \cdot \sin x$",
            "$5^x \cdot (\cos x + 1)$",
            title="Some basic functions",
            grid="both")
        # Set the second (right) axis scale as a log
        ax2.set_yscale("log")
        # Show the plot
        plt.show()

    .. image:: _static/dual_y_plot.png
       :width: 700
       :align: center

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the dual y-axis plot.
    ax1 : matplotlib.axes.Axes
        First axis, corresponding to the left y-axis (and containing the shared x-axis).
    ax2 : matplotlib.axes.Axes
        Second axis, corresponding to the right y-axis (and containing the shared x-axis).

    """
    # Setup a new figure
    fig, ax1 = plt.subplots(figsize=(11,6))
    # Add the common label on the x axis
    ax1.set_xlabel(x_label)
    # Add the label on the left y axis (in the correct color)
    ax1.set_ylabel(y_label_1, color=c1)
    # Plot the data for the left axis
    ax1.plot(x_data, y_data_1, color=c1)
    # Re-color the left y axis
    ax1.tick_params(axis="y", labelcolor=c1)
    if grid in ["right", "left", "both"]:
        # Add the grid based on the x-axis (common to both)
        ax1.grid(axis="x")
        if grid in ["left", "both"]:
            # Add the grid based on the left y-axis (dashed and in the correct color)
            ax1.grid(linestyle="dashed", color=c1, axis="y")
    # Create a second y axis that shares the same x axis
    ax2 = ax1.twinx()
    # Add the label on the right y axis (in the correct color)
    ax2.set_ylabel(y_label_2, color=c2)
    # Plot the data on the right axis
    ax2.plot(x_data, y_data_2, color=c2)
    # Re-color the right y axis
    ax2.tick_params(axis="y", labelcolor=c2)
    # Add the grid based on the right y-axis (dotted and in the correct color)
    if grid in ["right", "both"]:
        ax2.grid(linestyle="dotted", color=c2, axis="y")
    # Add a title
    if len(title) > 0:
        ax2.set_title(title)
    # Save space by using a tight layout
    fig.tight_layout()
    # Return the figure and axis
    return fig, ax1, ax2

def trajectory_3d(
    vehicles_states: Dict[float, np.ndarray],
    vehicles_names:List[str],
    central_body_name:str,
    spice_bodies:List[str] = [],
    frame_orientation:str = "J2000",
    center_plot:bool = False,
    colors:List[str] = [],
    linestyles:List[str] = [] ):
    """Plot the trajectory specified bodies in 3D.

    This function allows to plot the 3D trajectory of vehicles of which the state has been propagated, as well as
    the trajectory of bodies taken directly from the SPICE interface.

    Parameters
    ----------
    vehicles_states : Dict[float, numpy.ndarray]
        Dictionary mapping the simulation time steps to the propagated state time series of the simulation bodies.
    vehicles_names : List[str]
        List of the name of the simulated body for which the trajectory must be plotted. Use an empty string in the list to skip a specific body.
    central_body_name : str
        Name of the central body in the simulation
    spice_bodies : List[str], optional
        List of the name of bodies for which the trajectory has to be taken from SPICE and plotted.
    frame_orientation : str, optional, default="J2000"
        Orientation of the reference frame used in the simulation.
    center_plot : bool, optional, default=False
        If True, the central body will be centered on the plot.
    colors : List[str], optional
        List of colors to use for the trajectories. In order, the colors will first be used for the vehicles and then for the SPICE bodies.
    linestyles : List[str], optional
        List of linestyles to use for the trajectories. In order, the linestyles will first be used for the vehicles and then for the SPICE bodies.
    

    Examples
    --------
    After the propagation of two CubeSats on which thrust is applied, we can for instance plot both of their trajectories, as well as the trajectory of the Moon,
    using the following code snippet:

    .. code-block:: python

        # Plot the trajectory of two satellite and the Moon around the Earth
        fig, ax = plotting.trajectory_3d(
            vehicles_states=dynamics_simulator.state_history,
            vehicles_names=["Lunar CubeSat A", "Lunar CubeSat B"],
            central_body_name="Earth",
            spice_bodies=["Moon"],
            linestyles=["dotted", "dashed", "solid"],
            colors=["blue", "green", "grey"]
        )
        # Change the size of the figure
        fig.set_size_inches(8, 8)
        # Show the plot
        plt.show()

    .. image:: _static/trajectory_3D.png
       :width: 450
       :align: center

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the 3D plot.
    ax : matplotlib.axes.Axes
        3D axis system used to plot the trajectory.

    """
    # Import SPICE
    from ..kernel.interface import spice

    # Save color and linestyle index
    i_c, i_ls = 0, 0

    # Set color and linestyles if nothing specified
    if len(colors) == 0:
        colors = ["C%i" % i for i in range(10)]
    if len(linestyles) == 0:
        linestyles = ["-"]*(len(vehicles_names)+len(spice_bodies))

    # If the SPICE kernels are not loaded, load THEM
    if spice.get_total_count_of_kernels_loaded() < 4:
        spice.load_standard_kernels()

    # Make sure inputs are the correct format
    if type(vehicles_names) == str:
        vehicles_names = [vehicles_names]
    if type(spice_bodies) == str:
        spice_bodies = [spice_bodies]

    # Convert the states to a ndarray
    vehicles_states_array = result2array(vehicles_states)
    sim_epochs = vehicles_states_array[:,0]

    # Make a list of positions per vehicle
    vehicles_positions = []
    for i in range(len(vehicles_names)):
        vehicles_positions.append(vehicles_states_array[:, 1+i*6:4+i*6])

    # Create a figure with a 3D projection for the Moon and vehicle trajectory around Earth
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Save the minimum and maximum positions
    min_pos, max_pos = 0, 0

    # Plot the trajectory of each vehicle
    for i, vehicle_name in enumerate(vehicles_names):
        if len(vehicle_name) != 0:
            # Update the minimum and maximum positions
            min_pos, max_pos = min(min_pos, np.min(vehicles_positions[i])), max(max_pos, np.max(vehicles_positions[i]))
            # Select appropriate color and linestyle
            _color = colors[i_c]
            _linestyle = linestyles[i_ls]
            # Increment color and linestyle indexes
            i_c = i_c + 1 if i_c != len(colors) else 0
            i_ls = i_ls + 1 if i_ls != len(linestyles) else 0
            # Plot the trajectory of the vehicle
            ax.plot(vehicles_positions[i][:,0], vehicles_positions[i][:,1], vehicles_positions[i][:,2], label=vehicle_name, color=_color, linestyle=_linestyle)

    for i, spice_body in enumerate(spice_bodies):
        # Get the position of the body from SPICE
        body_state_array = np.array([
            spice.get_body_cartesian_position_at_epoch(spice_body, central_body_name, frame_orientation, "None", epoch)
            for epoch in sim_epochs
        ])
        # Update the minimum and maximum positions
        min_pos, max_pos = min(min_pos, np.min(body_state_array)), max(max_pos, np.max(body_state_array))
        # Select appropriate color and linestyle
        _color = colors[i_c]
        _linestyle = linestyles[i_ls]
        # Increment color and linestyle indexes
        i_c = i_c + 1 if i_c != len(colors) else 0
        i_ls = i_ls + 1 if i_ls != len(linestyles) else 0
        # Plot the trajectory of the body
        ax.plot(body_state_array[:,0], body_state_array[:,1], body_state_array[:,2], label=spice_body, color=_color, linestyle=_linestyle)

    # Plot the central body position
    ax.scatter(0.0, 0.0, 0.0, label=central_body_name, marker="o", color="k")

    # Make the plot centered if needed
    if center_plot:
        min_pos, max_pos = -max(np.fabs(min_pos), max_pos), max(np.fabs(min_pos), max_pos)

    # Add a legend, set the plot limits, and add axis labels
    ax.legend()
    ax.set_xlim([min_pos*1.2, max_pos*1.2]), ax.set_ylim([min_pos*1.2, max_pos*1.2]), ax.set_zlim([min_pos*1.2, max_pos*1.2])
    ax.set_xlabel("x [m]"), ax.set_ylabel("y [m]"), ax.set_zlabel("z [m]")

    # Use a tight layout for the figure (do last to avoid trimming axis)
    fig.tight_layout()

    # Return the figure and the axis system
    return fig, ax

def pareto_front(
    x_objective: list,
    y_objective: list,
    x_label: str,
    y_label: str,
    title:str = "",
    c_parameter:Union[list,None] = None,
    c_label:str = "",
    cmap:Union[matplotlib.colors.Colormap,str] = "viridis",
    operator:List[Union[min,max]] = [min, min],
    alpha:float = 0.85):
    """Plot dual-objective optimisation solutions with their associated Pareto front.

    This function allows to plot the cloud of points that results from an optimisation that contains two objectives, with the x and y axis both representing one of the optimisation objectives.
    A Pareto front is also automatically added to the plot, to represent the suite of optimum solutions.
    Optionally, information on a third design or optimisation parameter can also be incorporated to the plot in the form of a colormap.

    Parameters
    ----------
    x_objective : List[float]
        List of values representing the score of each solution according to a first optimisation objective.
    y_objective : List[float]
        List of values representing the score of each solution according to a second optimisation objective.
    x_label : str
        Label of the x-axis.
    y_label : str
        Label of the y-axis.
    title : str, optional
        Title of the plot.
    c_parameter : List[float], optional
        List of values representing an extra design or optimisation parameter, coloring the cloud of points.
    c_label : str, optional
        Label of the colorbar that is used to represent the optional third objective.
    cmap : str or matplotlib.colors.Colormap, optional, default="viridis"
        Colormap to use to represent the optional third objective.
    operator : List[min,max], optional, default=[min,min]
        List used to indicated if each of the two main objectives should be minimised or maximised.
    alpha : float, optional, default=0.85

    Examples
    -------
    The following example is based on the optimisation of the orbit of a CubeSat flying at very low altitude.
    The two main optmisation objectives were to minimise the mean orbital altitude, and minimise the loss in periapsis over a few years.
    The `pareto_front` function is then used to plot the optimisation results related to both of these objectives.
    In addition, a colormap is used to represent the mean power that the CubeSat was receiving in each of the computed orbits.

    .. code-block:: python

        # Generate fake result from an optimisation 
        mean_altitudes = np.random.gamma(100, 2, 1000)
        periapsis_decays = np.random.gamma(10, 5, 1000)
        mean_powers = np.random.normal(8, 3, 1000)

        # Plot the two main objectives, using colors to represent the mean power
        fig, ax = plotting.pareto_front(
            x_objective=mean_altitudes,
            y_objective=periapsis_decays,
            x_label="Mean altitude [km]",
            y_label="Periapsis decay [km]",
            title="CubeSat altitude vs periapsis decay objectives, as a function of power",
            c_parameter=mean_powers,
            c_label="Mean power [W]",
            cmap="plasma",
            alpha=0.65
        )

        # Show the plot
        plt.show()

    Running this code snippet results in a figure similar to the one below.
    Do note that the code snippet uses randomly generated data, while the figure shown below is the real optimisation result.

    .. image:: _static/pareto_front.png
       :width: 700
       :align: center

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the plot.
    ax : matplotlib.axes.Axes
        Axis used to make the plot.

    """
    if c_parameter is None:
        # Create the figure and its axis
        fig, ax = plt.subplots(figsize=(10,6))
    else:
        # Create the figure and its axis, with a 2.5% space reserved on its right for the colormap
        fig, (ax, cax) = plt.subplots(
            figsize=(10,6), 
            ncols=2,
            gridspec_kw={"width_ratios":[1, 0.025]}
        )
    # Plot the two main objectives, with the appropriate color if the color-associated objective is specified
    objectives_scatter = ax.scatter(x_objective, y_objective, alpha=alpha, c=c_parameter, cmap=cmap)
    # If a color-associated objective was specified, add a colormap on the side of the figure
    if c_parameter is not None:
        plt.colorbar(objectives_scatter, cax=cax, label=c_label)
        # Print a warning if the colorbar is left with no label (but accept this)
        if c_label == "":
            print("Warning: a colorbar has been added to the Pareto front plot, but no label was specified for it.")
    # Compute which of the points from the two main objectives are the Pareto optimums
    optimum_mask = pareto_optimums(np.array([x_objective, y_objective]).T, operator=operator)
    # Plot the pareto front (use the step function to not falsely indicate that a region is feasible without proof of it)
    ax.step(sorted(x_objective[optimum_mask], reverse=operator[0]==min), sorted(y_objective[optimum_mask], reverse=operator[1]==max), color="#418F3E", linewidth=2)
    # Add a label to both axis (this is not optional because a plot without axis label does not mean anything)
    ax.set_xlabel(x_label), ax.set_ylabel(y_label)
    # Add a grid
    ax.grid()
    # If a title is specified, add it
    if title != "":
        ax.set_title(title)
    # Use a tight layout to save space
    plt.tight_layout()
    # Return the figure and axis
    return fig, ax