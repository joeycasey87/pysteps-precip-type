# -*- coding: utf-8 -*-
"""
pysteps.visualization.precipTypefields
==================================

Methods for plotting precipitation types.

.. autosummary::
    :toctree: ../generated/

    plot_precip_field
    get_colormap
"""
import copy
import warnings

import matplotlib.pylab as plt
import numpy as np
from matplotlib import cm, colors

from pysteps.visualization.utils import get_geogrid, get_basemap_axis


def plot_precipType_field(
        precipType,
        ax=None,
        geodata=None,
        bbox=None,
        colorscale="pysteps",
        title=None,
        colorbar=True,
        cBarLabel="",
        categoryNr=4,
        axis="on",
        cax=None,
        map_kwargs=None,
):
    """
    Function to plot a precipitation types field with a colorbar.

    .. _Axes: https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes

    .. _SubplotSpec: https://matplotlib.org/api/_as_gen/matplotlib.gridspec.SubplotSpec.html

    Parameters
    ----------
    precipType: array-like
        Two-dimensional array containing the input precipitation types.
    geodata: dictionary or None, optional
        Optional dictionary containing geographical information about
        the field. Required is map is not None.

        If geodata is not None, it must contain the following key-value pairs:

        .. tabularcolumns:: |p{1.5cm}|L|

        +-----------------+---------------------------------------------------+
        |        Key      |                  Value                            |
        +=================+===================================================+
        |    projection   | PROJ.4-compatible projection definition           |
        +-----------------+---------------------------------------------------+
        |    x1           | x-coordinate of the lower-left corner of the data |
        |                 | raster                                            |
        +-----------------+---------------------------------------------------+
        |    y1           | y-coordinate of the lower-left corner of the data |
        |                 | raster                                            |
        +-----------------+---------------------------------------------------+
        |    x2           | x-coordinate of the upper-right corner of the     |
        |                 | data raster                                       |
        +-----------------+---------------------------------------------------+
        |    y2           | y-coordinate of the upper-right corner of the     |
        |                 | data raster                                       |
        +-----------------+---------------------------------------------------+
        |    yorigin      | a string specifying the location of the first     |
        |                 | element in the data raster w.r.t. y-axis:         |
        |                 | 'upper' = upper border, 'lower' = lower border    |
        +-----------------+---------------------------------------------------+
    bbox : tuple, optional
        Four-element tuple specifying the coordinates of the bounding box. Use
        this for plotting a subdomain inside the input grid. The coordinates are
        of the form (lower left x, lower left y ,upper right x, upper right y).
        If 'geodata' is not None, the bbox is in map coordinates, otherwise
        it represents image pixels.
    colorscale : {'pysteps', 'STEPS-BE', 'STEPS-NL', 'BOM-RF3'}, optional
        Which colorscale to use. TO BE DEFINED
    title : str, optional
        If not None, print the title on top of the plot.
    colorbar : bool, optional
        If set to True, add a colorbar on the right side of the plot.
    cBarLabel :
        Set color bar label.
    categoryNr :
        Number of categories to be plotted (2 to 6)
    axis : {'off','on'}, optional
        Whether to turn off or on the x and y axis.
    cax : Axes_ object, optional
        Axes into which the colorbar will be drawn. If no axes is provided
        the colorbar axes are created next to the plot.

    Other parameters
    ----------------
    map_kwargs: dict
        Optional parameters that need to be passed to
        :py:func:`pysteps.visualization.basemaps.plot_geography`.

    Returns
    -------
    ax : fig Axes_
        Figure axes. Needed if one wants to add e.g. text inside the plot.
    """

    if map_kwargs is None:
        map_kwargs = {}

    if len(precipType.shape) != 2:
        raise ValueError("The input is not two-dimensional array")

    # Assumes the input dimensions are lat/lon
    nlat, nlon = precipType.shape

    x_grid, y_grid, extent, regular_grid, origin = get_geogrid(
        nlat, nlon, geodata=geodata
    )

    ax = get_basemap_axis(extent, ax=ax, geodata=geodata, map_kwargs=map_kwargs)

    precipType = np.ma.masked_invalid(precipType)
    # plot rainfield
    if regular_grid:
        im = _plot_field(precipType, ax, colorscale, categoryNr, extent, origin=origin)
    else:
        im = _plot_field(
            precipType, ax, colorscale, categoryNr, extent, x_grid=x_grid, y_grid=y_grid
        )

    plt.title(title, loc='center', fontsize=25)

    # add colorbar
    if colorbar:
        # get colormap and color levels
        _, _, clevs, clevs_str = get_colormap(colorscale, categoryNr)
        cbar = plt.colorbar(
            im, ticks=clevs, spacing="uniform", extend="neither", shrink=0.8, cax=cax, drawedges=False
        )
        if clevs_str is not None:
            cbar.ax.set_yticklabels('')
            cbar.ax.tick_params(size=0)
            cbar.ax.set_yticks([i + .5 for i in clevs][:-1], minor=True)
            cbar.ax.set_yticklabels(clevs_str[:-1], minor=True, fontsize=15)
    cbar.set_label(cBarLabel)

    if geodata is None or axis == "off":
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])

    if bbox is not None:
        ax.set_xlim(bbox[0], bbox[2])
        ax.set_ylim(bbox[1], bbox[3])

    return ax


def _plot_field(precipType, ax, colorscale, categoryNr, extent, origin=None, x_grid=None, y_grid=None):
    precipType = precipType.copy()

    # Get colormap and color levels
    cmap, norm, _, _ = get_colormap(colorscale, categoryNr)

    if (x_grid is None) or (y_grid is None):
        im = ax.imshow(
            precipType,
            cmap=cmap,
            norm=norm,
            extent=extent,
            interpolation="nearest",
            origin=origin,
            zorder=10,
        )
    else:
        im = ax.pcolormesh(
            x_grid,
            y_grid,
            precipType,
            cmap=cmap,
            norm=norm,
            zorder=10,
        )

    return im


def get_colormap(colorscale="pysteps", categoryNr=4):
    """
    Function to generate a colormap (cmap) and norm.

    Parameters
    ----------
    colorscale : {'pysteps', 'STEPS-BE', 'STEPS-NL', 'BOM-RF3'}, optional
        Which colorscale to use. Applicable if units is 'mm/h', 'mm' or 'dBZ'.

    Returns
    -------
    cmap : Colormap instance
        colormap
    norm : colors.Normalize object
        Colors norm
    clevs: list(float)
        List of precipitation values defining the color limits.
    clevs_str: list(str)
        List of precipitation values defining the color limits (with correct
        number of decimals).
        :param categoryNr:
    """
    # Get list of colors
    color_list, clevs, clevs_str = _get_colorlist(colorscale, categoryNr)
    cmap = colors.LinearSegmentedColormap.from_list(
        "cmap", color_list, len(clevs) - 1
    )
    cmap.set_over("darkred", 1)
    cmap.set_bad("gray", alpha=0.5)
    cmap.set_under("none")
    norm = colors.BoundaryNorm(clevs, cmap.N)

    return cmap, norm, clevs, clevs_str


def _get_colorlist(colorscale="pysteps", categoryNr=4):
    """
    Function to get a list of colors to generate the colormap.

    Parameters
    ----------
    colorscale : str
        Which colorscale to use (BOM-RF3, pysteps, STEPS-BE, STEPS-NL)
    categoryNr  :
        How many categories should be plotted

    Returns
    -------
    color_list : list(str)
        List of color strings.

    clevs : list(float)
        List of precipitation values defining the color limits.

    clevs_str : list(str)
        List of precipitation type names
    """

    if categoryNr < 1 or categoryNr > 6:
        raise ValueError("Invalid category index [1 to 6] " + str(categoryNr))

    if colorscale == "pysteps":
        color_list = ["#ffe38f", "#ceda86", "#009489", "#3897ed", "#b0a0dc", "#ec623b"]
    # elif colorscale == 'other color scale': ... [6 colors]
    else:
        print("Invalid colorscale", colorscale)
        raise ValueError("Invalid colorscale " + colorscale)

    # Ticks and labels
    clevs = [1, 2, 3, 4, 5, 6, 7]
    clevs_str = ['Rain', 'Wet Snow', 'Snow', 'Freezing Rain', 'Hail', 'Severe Hail']

    # filter by category number
    color_list = color_list[0:categoryNr]
    clevs = clevs[0:(categoryNr + 1)]
    clevs_str = clevs_str[0:categoryNr] + ['']


    return color_list, clevs, clevs_str
