''' 
Copyright (c) 2010-2020, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

# General imports
import numpy as np
from scipy import ndimage

# Plotting imports
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Tudat imports
from tudatpy.kernel import constants
from tudatpy.kernel.astro import time_conversion


def plot_porkchop_of_single_field(
    departure_body: str,
    target_body: str,
    departure_epochs: np.ndarray, 
    arrival_epochs: np.ndarray, 
    DV: np.ndarray,
    C3: bool = False,
    threshold: float = 10,
    upscale: bool = True,
    # Plot arguments
    filled_contours: bool = True,
    plot_out_of_range_field: bool = True,
    plot_isochrones: bool = True,
    number_of_levels: int = 20,
    line_width: float = 0.5,
    line_color: str = 'black',
    font_size: float = 7,
    # Figure arguments
    colorbar: bool = False,
    percent_margin: float = 5,
    fig: matplotlib.figure.Figure = None,
    ax:  matplotlib.axes._axes.Axes = None,
    figsize: tuple[int, int] = (8, 8),
    show: bool = True,
    save: bool = False,
    filename: str = 'porkchop.png',
    ):

    #--------------------------------------------------------------------
    #%% DATA PREPARATION
    #--------------------------------------------------------------------

    # Transpose ΔV array such as to show departure on the horizontal axis and arrival on the vertical axis
    DV = DV.T
    
    # Divide ΔV array to display results in km/s
    DV = DV / 1000

    # Interpolate data to improve plot resolution
    if upscale:
        DV = ndimage.zoom(DV, 10)
        departure_epochs = np.linspace(departure_epochs.min(), departure_epochs.max(), DV.shape[1])
        arrival_epochs   = np.linspace(arrival_epochs.min(), arrival_epochs.max(), DV.shape[0])

    # Transform ΔV to characteristic energy (C3) if needed
    field = DV**2 if C3 else DV

    # Mask field array to discard excessive values
    field_within_range = np.ma.array(field, mask=field >= threshold)
    
    # Calculate departure and arrival time spans
    departure_epoch_span = (departure_epochs.max() - departure_epochs.min())
    arrival_epoch_span   = (arrival_epochs.max() - arrival_epochs.min())

    #--------------------------------------------------------------------
    #%% PLOT
    #--------------------------------------------------------------------

    # Use monospaced font for all text
    plt.rcParams["font.family"] = plt.rcParams["font.monospace"][0]
    plt.rcParams["mathtext.default"] = 'tt'

    # Create axis and figure
    if fig is None and ax is None:

        fig, ax = plt.subplots(figsize=figsize)

        # Title
        fig.suptitle(
            f'{departure_body}-{target_body} Transfer {"C3-Launch" if C3 else "ΔV-Launch"}',
            y=0.97
        )

    # Layout
    plt.subplots_adjust(
        top=0.92,
        bottom=0.155,
        left=0.15,
        right=0.875,
        hspace=0.2,
        wspace=0.2
    )

    # Filled contour
    if filled_contours:
        # Define custom colormap
        cmap = LinearSegmentedColormap.from_list(
            'TudatMap',
            ['#d1edff', '#0076C2'],
            N=256)
        contour = plt.contourf(
            departure_epochs,
            arrival_epochs,
            field_within_range,
            cmap=cmap,
            levels=np.logspace(np.log10(field_within_range.min()), np.log10(field_within_range.max()), number_of_levels),
            zorder=1.5)
    # Levels
    contour_lines = plt.contour(
        departure_epochs,
        arrival_epochs,
        field,
        colors=line_color,
        linewidths=line_width,
        levels=np.logspace(np.log10(field_within_range.min()), np.log10(field_within_range.max()), number_of_levels),
        zorder=1.5)
    plt.clabel(
        contour_lines,
        fontsize=font_size,
    )

    # Colorbar
    if colorbar:
        cax = fig.add_axes([
            ax.get_position().x1 + 0.01,
            ax.get_position().y0,
            0.02,
            ax.get_position().height]
        )
        cbar = fig.colorbar(
            contour,
            cax=cax)
        cbar.ax.set_title(
            '$m^2/s^2$' if C3 else 'm/s',
             pad=12.5,
             x=2.5
        )
        plt.sca(ax)
    
    # Out-of-range field plot
    if plot_out_of_range_field:
        if filled_contours:
            plt.contourf(
                departure_epochs,
                arrival_epochs,
                np.log(field),
                cmap='Reds',
                zorder=1.25)
        else:
            plt.contour(
                departure_epochs,
                arrival_epochs,
                np.log(field),
                cmap='Reds',
                levels=np.logspace(np.log10(field_within_range.min()), np.log10(field_within_range.max()), number_of_levels),
                linewidths=line_width,
                zorder=1.25)

    # Plot isochrones
    if plot_isochrones:
        x, y = np.meshgrid(departure_epochs, arrival_epochs)
        isoc = (y - x)
        isochrones = plt.contour(
            x, y, (isoc - isoc.min()) / constants.JULIAN_DAY,
            colors='black', alpha=0.5,
            levels=10,
            linewidths=2,
            linestyles='--',
            zorder=3)
        plt.clabel(
            isochrones,
            colors='black',
            fontsize=font_size + 3,
        )

    # Axis labels
    ax.set_xlabel(
        'Departure date',
        fontsize=font_size+3,
        weight='bold'
    )
    ax.set_ylabel(
        'Arrival date',
        fontsize=font_size+3,
        weight='bold'
    )

    # Axis limits
    ax.set_xlim([
        departure_epochs.min() - departure_epoch_span * percent_margin / 2 / 100,
        departure_epochs.max() + departure_epoch_span * percent_margin / 2 / 100,
    ])
    ax.set_ylim([
        arrival_epochs.min() - arrival_epoch_span * percent_margin / 2 / 100,
        arrival_epochs.max() + arrival_epoch_span * percent_margin / 2 / 100,
    ])

    # Tick labels
    dt = 40
    # X axis
    nx = int(np.floor(
        departure_epoch_span / (dt * constants.JULIAN_DAY))
    )
    x_ticks = np.linspace(departure_epochs.min(), departure_epochs.max(), nx)
    ax.xaxis.set_ticks(
        x_ticks,
        [f'{time_conversion.date_time_from_epoch(t).iso_string()[:10]}' for t in x_ticks],
        rotation=90
    )
    # Y axis
    ny = int(np.floor(
        arrival_epoch_span / (dt * constants.JULIAN_DAY))
    )
    y_ticks = np.linspace(arrival_epochs.min(), arrival_epochs.max(), ny)
    ax.yaxis.set_ticks(
        y_ticks,
        [f'{time_conversion.date_time_from_epoch(t).iso_string()[:10]}' for t in y_ticks]
    )

    # Grid
    plt.grid(True, linewidth=0.5, color='black')

    # Save
    if save:
        plt.savefig(filename, dpi=300)

    # Show
    if show:
        plt.show()


def plot_porkchop(
    departure_body: str,
    target_body: str,
    departure_epochs: np.ndarray, 
    arrival_epochs: np.ndarray, 
    DV: np.ndarray,
    C3: bool = False,
    total: bool = False,
    threshold: float = 10,
    upscale: bool = True,
    number_of_levels: int = 10,
    # Figure arguments
    percent_margin: float = 5,
    figsize: tuple[int, int] = (8, 8),
    show: bool = True,
    save: bool = False,
    filename: str = 'porkchop.png',
    ):
    
    if total:

        # Plot departure ΔV or C3
        plot_porkchop_of_single_field(
            departure_body           = departure_body,
            target_body              = target_body,
            departure_epochs         = departure_epochs,
            arrival_epochs           = arrival_epochs,
            DV                       = DV.sum(axis=2),
            C3                       = C3,
            threshold                = threshold,
            upscale                 = upscale,
            # Plot arguments
            line_width               = 0.5,
            line_color               = 'black',
            filled_contours          = True,
            plot_out_of_range_field  = True,
            number_of_levels         = number_of_levels,
            # Figure arguments
            colorbar                 = True,
            percent_margin           = percent_margin,
            figsize                  = figsize,
            show                     = False,
        )

    else:

        # Plot departure ΔV or C3
        plot_porkchop_of_single_field(
            departure_body           = departure_body,
            target_body              = target_body,
            departure_epochs         = departure_epochs,
            arrival_epochs           = arrival_epochs,
            DV                       = DV[:, :, 0],
            C3                       = C3,
            threshold                = threshold,
            upscale                 = upscale,
            # Plot arguments
            plot_out_of_range_field  = False,
            number_of_levels         = number_of_levels,
            # Figure arguments
            colorbar                 = True,
            percent_margin           = percent_margin,
            figsize                  = figsize,
            show                     = False,
        )

        # Plot departure ΔV or C3
        plot_porkchop_of_single_field(
            departure_body           = departure_body,
            target_body              = target_body,
            departure_epochs         = departure_epochs,
            arrival_epochs           = arrival_epochs,
            DV                       = DV[:, :, 1],
            C3                       = C3,
            threshold                = threshold,
            upscale                 = upscale,
            # Plot arguments
            line_width               = 0.5,
            line_color               = '#d90028',
            filled_contours          = False,
            plot_out_of_range_field  = False,
            plot_isochrones          = False,
            number_of_levels         = number_of_levels,
            # Figure arguments
            percent_margin           = percent_margin,
            fig                      = plt.gcf(),
            ax                       = plt.gca(),
            show                     = False,
        )

    # Save
    if save:
        plt.savefig(filename, dpi=300)

    # Show
    if show:
        plt.show()