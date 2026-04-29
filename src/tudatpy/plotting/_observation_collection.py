import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime
from tudatpy.estimation import observations

def plot_sky_map(
        collection: observations.ObservationCollection,
        object_designation: str | None = None,
        projection: str = None,
        figsize: tuple = (14, 7)
):
    """
    Generates a sky map from a Tudatpy ObservationCollection.

    Parameters:
    - collection: The Tudatpy ObservationCollection (C++ object).
    - designator: String name of the object for labeling.
    - projection: Matplotlib projection string (e.g., 'mollweide', 'aitoff', None).
    - figsize: Tuple for figure dimensions.
    """

    # 1. Extract data from the C++ Collection
    # Times are in TDB seconds since J2000
    times_tdb = np.array(collection.get_concatenated_observation_times())
    # Observations are a flat 1D array [RA1, Dec1, RA2, Dec2...]
    obs_values = np.array(collection.get_concatenated_observations())

    if len(times_tdb) == 0:
        print(f"Warning: No observations found for {object_designation}.")
        return None, None

    # 2. Reshape 1D observation array to Nx2 (RA, Dec)
    obs_values_2d = obs_values.reshape(-1, 2)
    ra_rad = obs_values_2d[:, 0]
    dec_rad = obs_values_2d[:, 1]

    # 3. Create the Figure
    fig, ax = plt.subplots(
        1, 1, subplot_kw={"projection": projection}, figsize=figsize
    )

    # 4. Prepare Time Scale for the Colorbar (TDB -> UTC)
    time_scale_converter = time_representation.default_time_scale_converter()
    utc_epochs = [
        time_scale_converter.convert_time(
            input_scale=time_representation.tdb_scale,
            output_scale=time_representation.utc_scale,
            input_value=t
        ) for t in times_tdb
    ]

    vmin, vmax = min(utc_epochs), max(utc_epochs)

    # 5. Handle Projections and Units
    if projection is not None:
        # Matplotlib projections (Mollweide/Aitoff) require Radians
        # and usually a domain of [-pi, pi]
        x_plot = (ra_rad + np.pi) % (2 * np.pi) - np.pi
        y_plot = dec_rad
    else:
        # Standard Cartesian plots use Degrees
        x_plot = np.degrees(ra_rad)
        y_plot = np.degrees(dec_rad)

    # 6. Scatter Plot
    sc = ax.scatter(
        x_plot,
        y_plot,
        s=30,
        marker="o",
        cmap=cm.plasma,
        label=object_designation,
        c=utc_epochs,
        vmin=vmin,
        vmax=vmax,
        edgecolors='black',
        linewidths=0.2
    )

    # 7. Axes Formatting
    if projection is None:
        ax.set_xlabel("Right Ascension [deg]")
        ax.set_ylabel("Declination [deg]")
    else:
        # For celestial projections, clean up labels
        ax.grid(True, linestyle=':', alpha=0.6)
        _format_projection_ticks(ax, projection)

    # 8. Colorbar with Date Formatting
    _add_date_colorbar(fig, ax, sc, vmin, vmax)

    # 9. Final Polish
    start_str = DateTime.from_epoch(vmin).to_iso_string().split(' ')[0]
    end_str = DateTime.from_epoch(vmax).to_iso_string().split(' ')[0]

    fig.suptitle(f"Object: {object_designation}\n{start_str} to {end_str}")
    ax.legend(loc='upper right', frameon=True, markerscale=1.5)
    fig.tight_layout()

    return fig, ax

def _add_date_colorbar(fig, ax, mappable, vmin, vmax):
    """Internal helper to add a colorbar with human-readable dates."""
    ticks = np.linspace(vmin, vmax, 6)
    cbar = fig.colorbar(mappable, ax=ax, label="Observation Date (UTC)", ticks=ticks, pad=0.05)

    # Convert numeric UTC seconds back to ISO Date strings (YYYY-MM-DD)
    date_labels = [
        DateTime.from_epoch(t).to_iso_string().split(' ')[0]
        for t in ticks
    ]
    cbar.ax.set_yticklabels(date_labels)

def _format_projection_ticks(ax, projection):
    """Internal helper to format ticks for Aitoff/Mollweide projections."""
    try:
        # Extract tick positions and convert to integer degrees
        # Note: we add 180 to RA to shift from [-180, 180] to [0, 360]
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()

        ax.set_xticklabels([f"{int(np.degrees(x) + 180)}°" for x in xticks])
        ax.set_yticklabels([f"{int(np.degrees(y))}°" for y in yticks])
    except:
        # Fallback if the projection doesn't support standard label setting
        pass