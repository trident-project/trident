"""
Spectrum plotting functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import \
    mylog
from matplotlib import \
    pyplot
import numpy as np

def plot_spectrum(wavelength, flux, filename="spectrum.png",
                  lambda_limits=None, flux_limits=None,
                  title=None, label=None, figsize=None, step=False,
                  stagger=0.2, features=None, axis_labels=None):
    """
    Plot a spectrum or a collection of spectra and save to disk.

    This function wraps some Matplotlib plotting functionality for
    plotting spectra generated with the :class:`~trident.SpectrumGenerator`.
    In its simplest form, it accepts a wavelength array consisting of 
    wavelength values and a corresponding flux array consisting of relative
    flux values, and it plots them and saves to disk.

    In addition, it can plot several spectra on the same axes simultaneously
    by passing a list of arrays to the ``wavelength``, ``flux`` arguments 
    (and optionally to the ``label`` and ``step`` keywords)..

    **Parameters**

    :wavelength: array of floats or list of arrays of floats

        Wavelength values in angstroms.  Either as an array of floats in the
        case of plotting a single spectrum, or as a list of arrays of floats
        in the case of plotting several spectra on the same axes.

    :flux: array of floats or list of arrays of floats

        Relative flux values (from 0 to 1) corresponding to wavelength array.
        Either as an array of floats in the case of plotting a single 
        spectrum, or as a list of arrays of floats in the case of plotting 
        several spectra on the same axes.

    :filename: string, optional

        Output filename of the plotted spectrum.  Will be a png file.
        Default: 'spectrum.png'

    :lambda_limits: tuple or list of floats, optional

        The minimum and maximum of the wavelength range (x-axis) for the plot
        in angstroms.  If specified as None, will use whole lambda range
        of spectrum. Example: (1200, 1400) for 1200-1400 Angstroms
        Default: None

    :flux_limits: tuple or list of floats, optional

        The minimum and maximum of the flux range (y-axis) for the plot.
        If specified as None, limits are automatically from
        [0, 1.1*max(flux)]. Example: (0, 1) for normal flux range before
        postprocessing.
        Default: None

    :step: boolean or list of booleans, optional

        Plot the spectrum as a series of step functions.  Appropriate for 
        plotting processed and noisy data.  Use a list of booleans when
        plotting multiple spectra, where each boolean corresponds to the entry
        in the ``wavelength`` and ``flux`` lists.

    :title: string, optional

        Optional title for plot
        Default: None

    :label: string or list of strings, optional

        Label for each spectrum to be plotted. Useful if plotting multiple
        spectra simultaneously.  Will automatically trigger a legend to be
        generated.
        Default: None

    :stagger: float, optional

        If plotting multiple spectra on the same axes, do we offset them in
        the y direction?  If set to None, no.  If set to a float, stagger them 
        by the flux value specified by this parameter.

    :features: dict, optional

        Include vertical lines with labels to represent certain spectral
        features.  Each entry in the dictionary consists of a key string to
        be overplot and the value float as to where in wavelength space it 
        will be plot as a vertical line with the corresponding label.

        Example: features={'Ly a' : 1216, 'Ly b' : 1026}

        Default: None

    :axis_labels: tuple of strings, optional

        Optionally set the axis labels directly.  If set to None, defaults to
        ('Wavelength [$\\AA$]', 'Relative Flux').
        Default: None

    **Example**

    Plot a flat spectrum

    >>> import numpy as np
    >>> import trident
    >>> wavelength = np.arange(1200, 1400)
    >>> flux = np.ones(len(wavelength))
    >>> trident.plot_spectrum(wavelength, flux)

    Generate a one-zone ray, create a Lyman alpha spectrum from it, and add
    gaussian noise to it.  Plot both the raw spectrum and the noisy spectrum
    on top of each other.

    >>> import trident
    >>> ray = trident.make_onezone_ray(column_densities={'H_p0_number_density':1e21})
    >>> sg = trident.SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    >>> sg.make_spectrum(ray, lines=['Ly a'])
    >>> sg.save_spectrum('spec.h5')
    >>> sg.add_gaussian_noise(10)
    >>> sg.save_spectrum('noise.h5')
    >>> sg1 = trident.SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    >>> sg2 = trident.SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    >>> sg1.load_spectrum('spec.h5')
    >>> sg2.load_spectrum('noise.h5')
    >>> trident.plot_spectrum([sg1.lambda_field, sg2.lambda_field], 
    ... [sg1.flux_field, sg2.flux_field], stagger=0, step=[False, True],
    ... filename='raw_and_noise.png')
    """

    # number of rows and columns
    n_rows = 1
    n_columns = 1

    # blank space between edge of figure and active plot area
    top_buffer = 0.07
    bottom_buffer = 0.15
    left_buffer = 0.06
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
    if figsize is None:
        figsize = (12, 4)
    pyplot.figure(figsize=figsize)

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
        labels = [label]
        steps = [step]
    else:
        wavelengths = wavelength
        fluxs = flux
        if label is not None: 
            labels = label
        else: 
            labels = [None]*len(fluxs)
        if step is not None: 
            steps = step
        else:
            steps = [None]*len(fluxs)

    # A running maximum of flux for use in ylim scaling in final plot
    max_flux = 0.

    for i, (wavelength, flux) in enumerate(zip(wavelengths, fluxs)):

        # Do we stagger the fluxes?
        if stagger is not None:
            flux -= stagger * i

        # Do we include labels and a legend?
        if steps[i]:
            my_axes.step(wavelength, flux, label=labels[i])
        else:
            my_axes.plot(wavelength, flux, label=labels[i])

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
    if axis_labels is None:
        axis_labels = ('Wavelength [$\\AA$]', 'Relative Flux')
    my_axes.xaxis.set_label_text(axis_labels[0])
    my_axes.yaxis.set_label_text(axis_labels[1])

    # Don't let the x-axis switch to offset values for tick labels
    my_axes.get_xaxis().get_major_formatter().set_useOffset(False)

    if label is not None: pyplot.legend()

    # Overplot the relevant features on the plot
    if features is not None:
        for feature in features:
            label = feature
            wavelength = features[feature]
            # Draw line
            my_axes.plot([wavelength, wavelength], flux_limits, '--', color='k')
            # Write text
            text_location = flux_limits[1] - 0.05*(flux_limits[1] - flux_limits[0])
            my_axes.text(wavelength, text_location, label, 
                    horizontalalignment='right', 
                    verticalalignment='top', rotation='vertical')
                    #transform=ax.transAxes) 

    #pyplot.tight_layout()
    mylog.info("Writing spectrum plot to png file: %s" % filename)
    pyplot.savefig(filename)
    pyplot.close()
