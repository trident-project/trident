import h5py
import numpy as np
import os
import roman
import sys
from yt.convenience import load
from yt.analysis_modules.absorption_spectrum.api import AbsorptionSpectrum
from yt.funcs import mylog, YTArray
from matplotlib import pyplot
from ion_balance import \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field, \
    atomic_mass
from line_database import LineDatabase

class SpectrumGenerator(AbsorptionSpectrum):
    def __init__(self, instrument=None, lambda_min=None, lambda_max=None,
                 n_lambda=None, dlambda=None, lsf_kernel=None, 
                 line_database='lines.txt', ionization_table=None):
        """
        SpectrumGenerator is a subclass of yt's AbsorptionSpectrum class
        with additional functionality like line lists, adding spectral
        templates, and plotting.

        Parameters
        ----------

        lambda_min, lambda_max : float
        The wavelength extrema in angstroms

        n_lambda : int
        The number of wavelength bins in the spectrum

        line_database : string, optional
        A text file listing the various lines to insert into the line database.
        The line database provides a list of all possible lines that could
        be added to the spectrum. The file should 4 tab-delimited columns of
        name (e.g. MgII), wavelength in angstroms, gamma of transition, and
        f-value of transition.  See example datasets in trident/data/line_lists
        for examples.

        ionization_table: hdf5 file, optional
        An HDF5 file used for computing the ionization fraction of the gas
        based on its density, temperature, metallicity, and redshift.
        The format of this file should be... <THIS NEEDS TO BE FINISHED>

        < IN GENERAL THESE DOCS ARE INCOMPLETE >

        """
        if instrument is None and lambda_min is None:
            instrument = 'COS'
            mylog.info("No parameters specified, defaulting to COS instrument.")
        elif instrument is None:
            instrument = Instrument(lambda_min=lambda_min,
                                    lambda_max=lambda_max,
                                    n_lambda=n_lambda,
                                    dlambda=dlambda,
                                    lsf_kernel=lsf_kernel, name="Custom")
        self.set_instrument(instrument)
        mylog.info("Setting instrument to %s" % self.instrument.name)

        AbsorptionSpectrum.__init__(self,
                                    self.instrument.lambda_min,
                                    self.instrument.lambda_max,
                                    self.instrument.n_lambda)
        
        # instantiate the LineDatabase
        self.line_database = LineDatabase(line_database)

        # store the ionization table in the SpectrumGenerator object
        if ionization_table is not None:
            # figure out where the user-specified files lives
            ionization_table = os.path.join(os.path.dirname(__file__), "..",
                                            "data", "ion_balance", filename)
            if not os.path.isfile(ionization_table):
                ionization_table = filename
            if not os.path.isfile(ionization_table):
                raise RuntimeError("ionization_table %s is not found in local "
                                   "directory or in trident/data/ion_balance "
                                   % (filename.split('/')[-1]))
            self.ionization_table = ionization_table
        else:
            table_dir = os.path.join(os.path.dirname(__file__), '../data/ion_balance')
            filelist = os.listdir(table_dir)
            ion_files = [i for i in filelist if i.endswith('.h5')]
            if 'hm2012_hr.h5' in ion_files: ionization_table = 'hm2012_hr.h5'
            elif 'hm2012_lr.h5' in ion_files: ionization_table = 'hm2012_lr.h5'
            else:
                mylog.info("No ionization file specified, using %s" %ion_files[0])
                ionization_table = ion_files[0]
            self.ionization_table = os.path.join(os.path.dirname(__file__), "..",
                                                 "data", "ion_balance", ionization_table)

    def make_spectrum(self, input_ds, lines=None,
                      output_file="spectrum.h5",
                      use_peculiar_velocity=True, njobs="auto"):
        """
        Make spectrum from ray data using the line list.

        Parameters
        ----------

        input_ds : string or dataset
           path to input ray data or a loaded ray dataset
        lines: FILL THIS IN
        output_file : optional, string
           path for output file.  File formats are chosen based on the
           filename extension.  ``.h5`` for hdf5, ``.fits`` for fits,
           and everything else is ASCII.
           Default: "spectrum.h5"
        use_peculiar_velocity : optional, bool
           if True, include line of sight velocity for shifting lines.
           Default: True
        njobs : optional, int or "auto"
           the number of process groups into which the loop over
           absorption lines will be divided.  If set to -1, each
           absorption line will be deposited by exactly one processor.
           If njobs is set to a value less than the total number of
           available processors (N), then the deposition of an
           individual line will be parallelized over (N / njobs)
           processors.  If set to "auto", it will first try to
           parallelize over the list of lines and only parallelize
           the line deposition if there are more processors than
           lines.  This is the optimal strategy for parallelizing
           spectrum generation.
           Default: "auto"
        """


        if isinstance(input_ds, str):
            input_ds = load(input_ds)
        ad = input_ds.all_data()

        active_lines = self.line_database.parse_subset(lines)
        
        # Make sure we've produced all the necessary
        # derived fields if they aren't native to the data
        for line in active_lines:
            try:
                disk_field = ad._determine_fields(line.field)[0]
            except:
                if line.field not in input_ds.derived_field_list:
                    my_ion = \
                      line.field[:line.field.find("number_density")]
                    on_ion = my_ion.split("_")
                    if on_ion[1]:
                        my_lev = int(on_ion[1][1:]) + 1
                    else:
                        my_lev = 1
                add_ion_number_density_field(on_ion[0], my_lev,
                                             self.ionization_table,
                                             input_ds)
            self.add_line(line.identifier, line.field,
                          float(line.wavelength),
                          float(line.f_value),
                          float(line.gamma),
                          atomic_mass[line.element],
                          label_threshold=1e3)

        AbsorptionSpectrum.make_spectrum(self, input_ds,
                                         output_file=output_file,
                                         line_list_file=None,
                                         use_peculiar_velocity=use_peculiar_velocity,
                                         njobs=njobs)

    def _get_qso_spectrum(self, redshift=0.0, filename=None):
        """
        Read in the composite QSO spectrum and return an interpolated version
        to fit the desired wavelength interval and binning.
        """

        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "..", "data",
                                    "spectral_templates",
                                    "qso_background_COS_HST.txt")

        data = np.loadtxt(filename)
        qso_lambda = YTArray(data[:, 0], 'angstrom')
        qso_lambda += qso_lambda * redshift
        qso_flux = data[:, 1]

        index = np.digitize(self.lambda_bins, qso_lambda)
        np.clip(index, 1, qso_lambda.size - 1, out=index)
        slope = (qso_flux[index] - qso_flux[index - 1]) / \
          (qso_lambda[index] - qso_lambda[index - 1])
        my_flux = slope * (self.lambda_bins - qso_lambda[index]) + qso_flux[index]
        return my_flux

    def _get_milky_way_foreground(self, filename=None):
        """
        Read in the composite QSO spectrum and return an interpolated version
        to fit the desired wavelength interval and binning.
        """

        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "..", "data",
                                    "spectral_templates",
                                    "mw_foreground_COS.txt")

        data = np.loadtxt(filename)
        MW_lambda = YTArray(data[:, 0], 'angstrom')
        MW_flux = data[:, 1]

        index = np.digitize(self.lambda_bins, MW_lambda)
        np.clip(index, 1, MW_lambda.size - 1, out=index)
        slope = (MW_flux[index] - MW_flux[index - 1]) / \
          (MW_lambda[index] - MW_lambda[index - 1])
        my_flux = slope * (self.lambda_bins - MW_lambda[index]) + MW_flux[index]
        # just set values that go beyond the data to 1
        my_flux[self.lambda_bins > 1799.9444] = 1.0
        return my_flux

    def add_milky_way_foreground(self, flux_field=None,
                                 filename=None):
        if flux_field is None:
            flux_field = self.flux_field
        MW_spectrum = self._get_milky_way_foreground(filename=filename)
        flux_field *= MW_spectrum

    def add_qso_spectrum(self, flux_field=None,
                         redshift=0.0, filename=None):
        if flux_field is None:
            flux_field = self.flux_field
        qso_spectrum = self._get_qso_spectrum(redshift=redshift,
                                              filename=filename)
        flux_field *= qso_spectrum

    def add_gaussian_noise(self, snr, n_bins=None, out=None, seed=None):
        np.random.seed(seed)
        if n_bins is None:
            n_bins = self.lambda_bins.size
        if out is None:
            out = self.flux_field
        np.add(out, np.random.normal(loc=0.0, scale=1/float(snr), size=n_bins),
               out=out)
        return out

    def apply_lsf(self, function=None, width=None, filename=None):
        """
        Apply the LSF to the flux_field of the spectrum.
        If an instrument already supplies a valid filename and no keywords
        are supplied, it is used by default.  Otherwise, the user can
        specify a filename of a user-defined kernel or a function+width
        for a kernel.  Valid functions are: "boxcar" and "gaussian".
        """
        # if nothing is specified, then use the Instrument-defined kernel
        if function is None and width is None and filename is None:
            if self.instrument.lsf_kernel is None:
                raise RuntimeError("To apply a line spread function, you "
                                   "must specify one or use an instrument "
                                   "where one is defined.")
            else:
                mylog.info("Applying default line spread function for %s." % \
                       self.instrument.name)
                lsf = LSF(filename=self.instrument.lsf_kernel)
        else:
            mylog.info("Applying specified line spread function.")
            lsf = LSF(function=function, width=width, filename=filename)
        self.flux_field = np.convolve(lsf.kernel,self.flux_field,'same')

    def load_spectrum(self, filename=None):
        if not filename.endswith(".h5"):
            raise RuntimeError("Only hdf5 format supported for loading spectra.")
        in_file = h5py.File(filename, "r")
        self.lambda_bins = in_file['wavelength'].value
        self.flux_field = in_file['flux'].value
        in_file.close()

    def make_flat_spectrum(self):
        """
        Makes a flat spectrum devoid of any lines.
        """
        self.flux_field = np.ones(self.lambda_bins.size)
        return (self.lambda_bins, self.flux_field)

    def set_instrument(self, instrument):
        """
        Sets the appropriate range of wavelengths and binsize for the
        output spectrum as well as the line spread function.

        set_instrument accepts either the name of a valid instrument or
        a fully specified Instrument object.

        Valid instruments are: %s
        """ % valid_instruments.keys()

        if isinstance(instrument, str):
            if instrument not in valid_instruments:
                raise RuntimeError("set_instrument accepts only Instrument "
                                   "objects or the names of valid "
                                   "instruments: ", valid_instruments.keys())
            self.instrument = valid_instruments[instrument]
        elif isinstance(instrument, Instrument):
            self.instrument = instrument
        else:
            raise RuntimeError("set_instrument accepts only Instrument "
                               "objects or the names of valid instruments: ",
                               valid_instruments.keys())

    def add_line_to_database(self, element, ion_state, wavelength, gamma, 
                             f_value, field=None, identifier=None):
        self.line_database.add_line(element, ion_state, wavelength,
                                    gamma, f_value, field=field,
                                    identifier=identifier)

class Instrument():
    """
    An instrument template for specifying a spectrograph/telescope pair
    """
    def __init__(self, lambda_min, lambda_max, n_lambda=None,
                 dlambda=None, lsf_kernel=None, name=None):
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.lsf_kernel = lsf_kernel
        if n_lambda is None and dlambda is None:
            raise RuntimeError("Either n_lambda or dlambda must be set to "
                               "specify the binsize")
        elif dlambda is not None:
            n_lambda = (lambda_max - lambda_min) / dlambda
        self.n_lambda = n_lambda
        if name is not None:
            self.name = name

class LSF():
    """
    Line Spread Function class

    The user must define either a filename or a function and a width

    Parameters
    ----------

    function : string, optional
        the function defining the LSF kernel.
        valid functions are "boxcar" or "gaussian"

    width : int, optional
        the width of the LSF kernel

    filename : string , optional
        the filename of a textfile for a user-specified kernel. each line
        in the textfile is the non-normalized flux value of the kernel
    """
    def __init__(self, function=None, width=None, filename=None):
        self.kernel = []
        # if filename is defined, use it
        if filename is not None:
            # Check to see if the file is in the local dir
            if os.path.isfile(filename):
                lsf_file = open(filename, 'r')
            # otherwise use the file in the lsf_kernels dir
            else:
                filename2 = os.path.join(os.path.dirname(__file__), "..",
                                         "data", "lsf_kernels", filename)
                if os.path.isfile(filename2):
                    lsf_file = open(filename2, 'r')
                else:
                    sys.exit("filename must be in local directory or in",
                             "trident/data/lsf_kernels directory")
            for line in lsf_file:
                self.kernel.append(float(line.split()[1]))
            lsf_file.close()
        elif function is not None and width is not None:
            if function == 'boxcar':
               self.kernel = np.ones(width)
            #XXX Define more functional forms
            #elif function == 'gaussian':
            #   self.kernel = np.gaussian(width)
        else:
            sys.exit("Either filename OR function+width must be specified.")

def plot_spectrum(wavelength, flux, filename="spectrum.png",
                  lambda_limits=None, flux_limits=None, 
                  title=None, label=None,
                  stagger=0.2):
    """
    Plot a spectrum or a collection of spectra and save to disk

    Parameters
    ----------
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

# Valid instruments
valid_instruments = \
    {'COS' :
       Instrument(1150, 1450, dlambda=0.01, lsf_kernel='avg_COS.txt', name='COS'),
     'HIRES' :
       Instrument(1200, 1400, dlambda=0.01, name='HIRES'),
     'UVES' :
       Instrument(1200, 1400, dlambda=0.01, name='UVES'),
     'MODS' :
       Instrument(1200, 1400, dlambda=0.01, name='MODS'),
     'SDSS' :
       Instrument(1200, 1400, dlambda=0.01, name='SDSS')}
