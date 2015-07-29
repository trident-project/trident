import h5py
import numpy as np
import os
from yt.analysis_modules.absorption_spectrum.api import \
      AbsorptionSpectrum
from yt.funcs import mylog, YTArray
from matplotlib import pyplot
import sys
      
atomic_mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941,
               'Be': 9.012182, 'B': 10.811, 'C': 12.0107,
               'N': 14.0067, 'O': 15.9994, 'F': 18.9984032,
               'Ne': 20.1797, 'Na': 22.989770, 'Mg': 24.3050,
               'Al': 26.981538, 'Si': 28.0855, 'P': 30.973761,
               'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
               'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955910,
               'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961,
               'Mn': 54.938049, 'Fe': 55.845, 'Co': 58.933200,
               'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409}


class SpectrumGenerator(AbsorptionSpectrum):
    def __init__(self, instrument=None, lambda_min=None, lambda_max=None, 
                 n_lambda=None, dlambda=None, lsf_kernel=None, line_list=None):
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

        line_list : string, optional
        A text file listing the various lines to deposit in the spectrum.
        File is 4 tab-delimited columns of name (e.g. MgII), wavelength in
        angstroms, gamma of transition, and f-value of transition.  See
        example datasets in trident/data/line_lists for examples.

        """
        if instrument is None and lambda_min is None:
            instrument = 'COS'
            print "No parameters specified, defaulting to COS instrument"
        elif instrument is None:
            instrument = Instrument(lambda_min=lambda_min, 
                                    lambda_max=lambda_max,
                                    n_lambda=n_lambda,
                                    dlambda=dlambda,
                                    lsf_kernel=lsf_kernel, name="Custom")
        self.set_instrument(instrument)
        print "Setting instrument to %s" % self.instrument.name

        AbsorptionSpectrum.__init__(self, 
                                    self.instrument.lambda_min,
                                    self.instrument.lambda_max,
                                    self.instrument.n_lambda)
        # Load a line list by default if one is provided
        self.load_line_list(line_list)

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
                print "Applying default line spread function for %s." % \
                       self.instrument.name
                lsf = LSF(filename=self.instrument.lsf_kernel)
        else:
            print "Applying specified line spread function."
            lsf = LSF(function=function, width=width, filename=filename)
        self.flux_field = np.convolve(lsf.kernel,self.flux_field,'same')

    def load_spectrum(self, filename=None):
        if not filename.endswith(".h5"):
            raise RuntimeError("Only hdf5 format supported for loading spectra.")
        in_file = h5py.File(filename, "r")
        self.lambda_bins = in_file['wavelength'].value
        self.flux_field = in_file['flux'].value
        in_file.close()
    
    def load_line_list(self, filename=None):
        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "..", "data",
                                    "line_lists",
                                    "Nist_elem_list.txt")
        else:
            # check to see if file exists, if not, check in 
            # trident/data/line_lists
            if not os.path.isfile(filename):
                filename = os.path.join(os.path.dirname(__file__), "..", 
                                        "data", "line_lists", filename)
            if not os.path.isfile(filename):
                raise RuntimeError("line_list %s is not found in local "
                                   "directory or in trident/data/line_lists " 
                                   % (filename.split('/')[-1]))

        for line in file(filename).readlines():
            online = line.split()
            if line.startswith("#") or "--" in line or len(online) != 4: continue
            list_ion, wavelength, gamma, f_value = online
            label = list_ion
            ion = list_ion
            if "Ly" in list_ion:
                ion = "HI"
                label = "HI %s" % list_ion
            if "*" in ion:
                ion = ion[:ion.find("*")]
            element = ion[:2]
            if element[1].isupper():
                element = element[:1]        

            # # This is how things should work if there are H_number_density
            # # fields, but I cam commenting it out for now and we'll just
            # # use the post-processed fields
            if "Ly" in list_ion:
                field = "H_number_density" 
            else:
                field = "%s_Cloudy_eq_NumberDensity_post" % ion
                
            self.add_line(label, field, float(wavelength),
                          float(f_value), float(gamma),
                          atomic_mass[element], label_threshold=1e3)
        mylog.info("Load %d lines from %s." % (len(self.line_list), filename))

    def _write_spectrum_ascii(self, filename):
        print "Writing spectrum to ascii file: %s." % filename
        f = open(filename, 'w')
        f.write("# wavelength[A] flux rootflux\n")
        for i in xrange(self.lambda_bins.size):
            f.write("%e %e %e\n" % (self.lambda_bins[i],
                                    self.flux_field[i], 
                                    self.flux_field[i]**0.5))
        f.close()

    def _write_spectrum_fits(self, filename):
        """
        Write spectrum to a fits file.
        """
        try:
            import pyfits
        except:
            print "Could not import the pyfits module.  Please install pyfits."
            return

        print "Writing spectrum to fits file: %s." % filename
        col1 = pyfits.Column(name='wavelength', format='E', array=self.lambda_bins)
        col2 = pyfits.Column(name='flux', format='E', array=self.flux_field)
        col3 = pyfits.Column(name='rootflux', format='E', array=self.flux_field**0.5)
        cols = pyfits.ColDefs([col1, col2, col3])
        tbhdu = pyfits.new_table(cols)
        tbhdu.writeto(filename, clobber=True)

    def _write_spectrum_hdf5(self, filename):
        """
        Write spectrum to an hdf5 file.

        """
        print "Writing spectrum to hdf5 file: %s." % filename
        output = h5py.File(filename, 'w')
        output.create_dataset('wavelength', data=self.lambda_bins)
        output.create_dataset('flux', data=self.flux_field)
        output.create_dataset('rootflux', data=self.flux_field**0.5)
        output.close()

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
                  lambda_min=None, lambda_max=None, title=None, label=None, 
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
            
    if lambda_min is None and lambda_max is None:
        my_axes.set_xlim(wavelength.min(), wavelength.max())
    else:
        my_axes.set_xlim(lambda_min, lambda_max)
    my_axes.set_ylim(0, 2)
    my_axes.xaxis.set_label_text("$\\lambda$ [$\\AA$]")
    my_axes.yaxis.set_label_text("Relative Flux")
    if label is not None: pyplot.legend()
    pyplot.savefig(filename)

# Valid instruments
valid_instruments = \
    {'COS' : 
       Instrument(1150, 1450, dlambda=0.005, lsf_kernel='avg_COS.txt', name='COS'),
     'HIRES' :
       Instrument(1200, 1400, dlambda=0.01, name='HIRES'),
     'UVES' :
       Instrument(1200, 1400, dlambda=0.01, name='UVES'),
     'MODS' :
       Instrument(1200, 1400, dlambda=0.01, name='MODS'),
     'SDSS' :
       Instrument(1200, 1400, dlambda=0.01, name='SDSS')}

