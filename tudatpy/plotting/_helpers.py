import matplotlib.pyplot as plt


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
    fig : matplotlib.pyplot.figure
        Figure containing the dual y-axis plot.
    ax1 : matplotlib.pyplot.axis
        First axis, corresponding to the left y-axis (and containing the shared x-axis).
    ax2 : matplotlib.pyplot.axis
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