"""
Spectrum plotting functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from matplotlib import pyplot
import numpy as np

def plot_spectrum(wavelength, flux, filename="spectrum.png",
                  lambda_limits=None, flux_limits=None,
                  title=None, label=None,
                  stagger=0.2):
    """
    Plot a spectrum or a collection of spectra and save to disk

    Parameters

    wavelength : array or list of arrays
        wavelength vals in angstroms

    flux : array or list of arrays
        relative flux (from 0 to 1)

    filename : string, optional

    title : string, optional
        title for plot

    label : string or list of strings, optional
        label for each spectrum to be plotted

    stagger : float, optional
        if plotting multiple spectra on top of each other, do we stagger them?
        If None, no.  If set to a float, it is the value in relative flux to
        stagger each spectrum
    """

    # number of rows and columns
    n_rows = 1
    n_columns = 1

    # blank space between edge of figure and active plot area
    top_buffer = 0.07
    bottom_buffer = 0.15
    left_buffer = 0.05
    right_buffer = 0.03

    # blank space between plots
    hor_buffer = 0.05
    vert_buffer = 0.05

    # calculate the height and width of each panel
    panel_width = ((1.0 - left_buffer - right_buffer -
                    ((n_columns-1)*hor_buffer)) / n_columns)
    panel_height = ((1.0 - top_buffer - bottom_buffer -
                     ((n_rows-1)*vert_buffer)) / n_rows)

    # create a figure (figsize is in inches)
    pyplot.figure(figsize=(12, 4))

    # get the row and column number
    my_row = 0
    my_column = 0

    # calculate the position of the bottom, left corner of this plot
    left_side = left_buffer + (my_column * panel_width) + \
                my_column * hor_buffer
    top_side = 1.0 - (top_buffer + (my_row * panel_height) + \
               my_row * vert_buffer)
    bottom_side = top_side - panel_height

    # create an axes object on which we will make the plot
    my_axes = pyplot.axes((left_side, bottom_side, panel_width, panel_height))

    # Are we overplotting several spectra?  or just one?
    if not (isinstance(wavelength, list) and isinstance(flux, list)):
        wavelengths = [wavelength]
        fluxs = [flux]
        if label is not None: labels = [label]
    else:
        wavelengths = wavelength
        fluxs = flux
        if label is not None: labels = label

    # A running maximum of flux for use in ylim scaling in final plot
    max_flux = 0.

    for i, (wavelength, flux) in enumerate(zip(wavelengths, fluxs)):

        # Do we stagger the fluxes?
        if stagger is not None:
            flux -= stagger * i

        # Do we include labels and a legend?
        if label is not None:
            my_axes.plot(wavelength, flux, label=labels[i])
        else:
            my_axes.plot(wavelength, flux)

        # Do we include a title?
        if title is not None:
            pyplot.title(title)

        new_max_flux = np.max(flux)
        if new_max_flux > max_flux:
            max_flux = new_max_flux

    if lambda_limits is None:
        lambda_limits = (wavelength.min(), wavelength.max())
    my_axes.set_xlim(lambda_limits[0], lambda_limits[1])

    if flux_limits is None:
        flux_limits = (0, 1.1*max_flux)
    my_axes.set_ylim(flux_limits[0], flux_limits[1])
    my_axes.xaxis.set_label_text("$\\lambda$ [$\\AA$]")
    my_axes.yaxis.set_label_text("Relative Flux")

    # Don't let the x-axis switch to offset values for tick labels
    my_axes.get_xaxis().get_major_formatter().set_useOffset(False)

    if label is not None: pyplot.legend()
    pyplot.savefig(filename)
