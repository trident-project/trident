"""
SpectrumGenerator class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import os

from yt.analysis_modules.absorption_spectrum.api import \
    AbsorptionSpectrum
from yt.convenience import \
    load
from yt.funcs import \
    mylog, \
    YTArray

from trident.instrument import \
    Instrument
from trident.ion_balance import \
    add_ion_number_density_field, \
    atomic_mass
from trident.line_database import \
    LineDatabase
from trident.lsf import \
    LSF
from trident.plotting import \
    plot_spectrum
from trident.utilities import \
    trident_path

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

class SpectrumGenerator(AbsorptionSpectrum):
    """
    SpectrumGenerator is a subclass of yt's AbsorptionSpectrum class
    with additional functionality like line lists, adding spectral
    templates, and plotting.

    **Parameters**

    lambda_min, lambda_max : int
        The wavelength extrema in angstroms
        Defaults: None

    n_lambda : int
        The number of wavelength bins in the spectrum
        Default: None

    dlambda : float
        The desired wavelength bin width of the spectrum (in angstroms)

    lsf_kernel : string, optional
        The filename for the LSF kernel

    line_database : string, optional
        A text file listing the various lines to insert into the line database.
        The line database provides a list of all possible lines that could
        be added to the spectrum. The file should 4 tab-delimited columns of
        name (e.g. MgII), wavelength in angstroms, gamma of transition, and
        f-value of transition.  See example datasets in trident/data/line_lists
        for examples.
        Default: lines.txt

    ionization_table: hdf5 file, optional
        An HDF5 file used for computing the ionization fraction of the gas
        based on its density, temperature, metallicity, and redshift.  
    """
    def __init__(self, instrument=None, lambda_min=None, lambda_max=None,
                 n_lambda=None, dlambda=None, lsf_kernel=None,
                 line_database='lines.txt', ionization_table=None):
        if instrument is None and lambda_min is None:
            instrument = 'COS'
            mylog.info("No parameters specified, defaulting to COS instrument.")
        elif instrument is None:
            instrument = Instrument(lambda_min=lambda_min,
                                    lambda_max=lambda_max,
                                    n_lambda=n_lambda,
                                    dlambda=dlambda,
                                    lsf_kernel=lsf_kernel, name="Custom")
        self.observing_redshift = 0.
        self.set_instrument(instrument)
        mylog.info("Setting instrument to %s" % self.instrument.name)
        self.dlambda = self.instrument.dlambda

        AbsorptionSpectrum.__init__(self,
                                    self.instrument.lambda_min,
                                    self.instrument.lambda_max,
                                    self.instrument.n_lambda)

        # instantiate the LineDatabase
        self.line_database = LineDatabase(line_database)

        # Instantiate the spectrum to be zeros and ones for tau_field and
        # flux_field respectively.
        self.clear_spectrum()

        if ionization_table is not None:
            # figure out where the user-specified files lives
            ion_table_dir, ion_table_file = parse_config()
            ion_table_filepath = os.path.join(ion_table_dir, ion_table_file)
            if os.path.isfile(ion_table_file):
                self.ionization_table = ion_table_file
            elif os.path.isfile(ion_table_filepath):
                self.ionization_table = ion_table_filepath
            else:
                raise RuntimeError("ionization_table %s is not found in local "
                                   "directory or in %s" % 
                                   (filename.split('/')[-1], 
                                    ion_table_dir))
        else:
            self.ionization_table = None

    def make_spectrum(self, input_ds, lines=None,
                      output_file=None,
                      use_peculiar_velocity=True, 
                      observing_redshift=0.0,
                      njobs="auto"):
        """
        Make spectrum from ray data using the line list.

        **Parameters**

        input_ds : string or dataset
            Path to input ray data or a loaded ray dataset
        lines: list of strings
            List of strings that determine which lines will be added
            to the spectrum.  List can include things like "C", "O VI",
            or "Mg II ####", where #### would be the integer wavelength
            value of the desired line.
        output_file : optional, string
            Path for output file.  File formats are chosen based on the
            filename extension.  ".h5" for HDF5, ".fits" for FITS,
            and everything else is ASCII.
            Default: None
        use_peculiar_velocity : optional, bool
            If True, include line of sight velocity for shifting lines.
            Default: True
        observing_redshift : optional, float
            This is the value of the redshift at which the observer of this
            spectrum exists.  In most cases, this will be a redshift of 0.
            Default: 0.
        njobs : optional, int or "auto"
            The number of process groups into which the loop over
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
        self.observing_redshift = observing_redshift

        if isinstance(input_ds, str):
            input_ds = load(input_ds)
        ad = input_ds.all_data()

        # Clear out any previous spectrum that existed first
        self.clear_spectrum()

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
                add_ion_number_density_field(on_ion[0], my_lev, input_ds, 
                                             self.ionization_table)
            self.add_line(line.identifier, line.field,
                          float(line.wavelength),
                          float(line.f_value),
                          float(line.gamma),
                          atomic_mass[line.element],
                          label_threshold=1e3)

        AbsorptionSpectrum.make_spectrum(self, input_ds,
                                         output_file=None,
                                         line_list_file=None,
                                         use_peculiar_velocity=use_peculiar_velocity,
                                         observing_redshift=observing_redshift,
                                         njobs=njobs)

    def _get_qso_spectrum(self, emitting_redshift, observing_redshift, 
                          filename=None):
        """
        Read in the composite QSO spectrum and return an interpolated version
        to fit the desired wavelength interval and binning.
        """
        if observing_redshift is None:
            observing_redshift = self.observing_redshift
        if emitting_redshift is None:
            emitting_redshift = 0.
        # Following Hogg (2000) eq. 13 for the effective redshift z12 of 
        # observing at z1 redshift light emitted at z2:
        # 1 + z12 = (1 + z2) / (1 + z1)
        redshift_eff = (1 + emitting_redshift) / (1 + observing_redshift) - 1

        if filename is None:
            filename = os.path.join(trident_path(),  "data",
                                    "spectral_templates",
                                    "qso_background_COS_HST.txt")

        data = np.loadtxt(filename)
        qso_lambda = YTArray(data[:, 0], 'angstrom')
        qso_lambda += qso_lambda * redshift_eff
        qso_flux = data[:, 1]

        index = np.digitize(self.lambda_field, qso_lambda)
        np.clip(index, 1, qso_lambda.size - 1, out=index)
        slope = (qso_flux[index] - qso_flux[index - 1]) / \
          (qso_lambda[index] - qso_lambda[index - 1])
        my_flux = slope * (self.lambda_field - qso_lambda[index]) + qso_flux[index]
        return my_flux

    def _get_milky_way_foreground(self, filename=None):
        """
        Read in the Milky Way foreground spectrum and return an interpolated 
        version to fit the desired wavelength interval and binning.
        """

        if filename is None:
            filename = os.path.join(trident_path(), "data",
                                    "spectral_templates",
                                    "mw_foreground_COS.txt")

        data = np.loadtxt(filename)
        MW_lambda = YTArray(data[:, 0], 'angstrom')
        MW_flux = data[:, 1]

        index = np.digitize(self.lambda_field, MW_lambda)
        np.clip(index, 1, MW_lambda.size - 1, out=index)
        slope = (MW_flux[index] - MW_flux[index - 1]) / \
          (MW_lambda[index] - MW_lambda[index - 1])
        my_flux = slope * (self.lambda_field - MW_lambda[index]) + MW_flux[index]
        # just set values that go beyond the data to 1
        my_flux[self.lambda_field > 1799.9444] = 1.0
        return my_flux

    def add_milky_way_foreground(self, flux_field=None,
                                 filename=None):
        """
        Add a Milky Way foreground flux to the spectrum.  Data
        from Charles Danforth. Median-filter of 92 normalized 
        COS/G130M+G160M AGN spectra spanning the wavelength range of 
        1130 to 1800 Angstroms in 0.07 Angstrom bin size.

        **Parameters**

        flux_field : optional, array
            Array of flux values to which the Milky Way foreground is applied.
            Default: None
        filename : string
            Filename where the Milky Way foreground values used to modify
            the flux are stored.
            Default: None
        """
        if flux_field is None:
            flux_field = self.flux_field
        MW_spectrum = self._get_milky_way_foreground(filename=filename)
        flux_field *= MW_spectrum

    def add_qso_spectrum(self, flux_field=None,
                         emitting_redshift=None,
                         observing_redshift=None,
                         filename=None):
        """
        Add a composite QSO spectrum to the spectrum.  Uses data from 
        Telfer et al., ApJ, 565, 773 
        "The Rest-Frame Extreme Ultraviolet Spectral Properties of QSO"
        HST Radio Quiet composite for < 1275 Ang, SDSS composite > 2000 Ang,
        mean in between 8251 0 

        **Parameters**

        flux_field : array, optional
            Array of flux values to which the quasar background is applied.
            Default: None
        emitting_redshift: float, optional
            Redshift value at which the QSO emitted its light.  If specified
            as None, use 0.
            Default: None
        observing_redshift: float, optional
            Redshift value at which the quasar is observed.  If specified as
            None, use the observing_redshift value specified in make_spectrum()
            which defaults to 0.
            Default: None
        filename : string
            Filename where the Milky Way foreground values used to modify
            the flux are stored.
            Default: None
        """
        if flux_field is None:
            flux_field = self.flux_field
        qso_spectrum = self._get_qso_spectrum(emitting_redshift=emitting_redshift,
                                              observing_redshift=observing_redshift,
                                              filename=filename)
        flux_field *= qso_spectrum

    def add_gaussian_noise(self, snr, n_bins=None, out=None, seed=None):
        """
        Add random gaussian noise to the spectrum.

        **Parameters**

        snr : int
            The desired signal-to-noise ratio for adding the gaussian noise
        n_bins: int
            <I'm not entirely sure what functionality this has>
            Default: None
        out : array
            Array of flux values to which the noise will be added
            <note from devin: should we rename this something more intuitive?>
            Default: None
        seed : optional, int
            Seed for the random number generator.  This should be used to
            ensure than the same noise is added each time the spectrum is
            regenerated, if desired.
            Default: None

        """
        np.random.seed(seed)
        if n_bins is None:
            n_bins = self.lambda_field.size
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

        **Parameters**

        function : string, optional
            Desired functional form for the applied LSF kernel.
            Valid options are currently "boxcar" or "gaussian"
            Default: None
        width : int, optional
            Width of the desired LSF kernel in bin elements
            Default: None
        filename : string, optional
            The filename of the user-supplied kernel for applying the LSF
            Default: None
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
        from astropy.convolution import convolve
        self.flux_field = convolve(self.flux_field, lsf.kernel)

    def load_spectrum(self, filename=None):
        """
        Load a previously generated spectrum.

        **Parameters**

        filename : string
            The HDF5 file from which the previously generated spectrum
            should be read.  Note: only HDF5 files can currently be reloaded.
            Default: None
        """
        if not filename.endswith(".h5"):
            raise RuntimeError("Only hdf5 format supported for loading spectra.")
        in_file = h5py.File(filename, "r")
        self.lambda_field = in_file['wavelength'].value
        self.flux_field = in_file['flux'].value
        in_file.close()

    def clear_spectrum(self):
        """
        Clears the existing spectrum's flux and tau fields and replaces them
        with ones and zeros respectively.  Clear the line list kept in
        the AbsorptionSpectrum object as well.  Also clear the line_subset
        stored by the LineDatabase.
        """
        # Set flux and tau to ones and zeros
        self.flux_field = np.ones(self.lambda_field.size)
        self.tau_field = np.zeros(self.lambda_field.size)

        # Clear out the line list that is stored in AbsorptionSpectrum
        self.line_list = []

        # Make sure we reset the line database as well
        self.line_database.lines_subset = []

    def set_instrument(self, instrument):
        """
        Sets the appropriate range of wavelengths and binsize for the
        output spectrum as well as the line spread function.

        set_instrument accepts either the name of a valid instrument or
        a fully specified Instrument object.

        Valid instruments are: %s

        **Parameters**

        instrument : instrument object
            The instrument object that should be used to create the spectrum.
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
        """
        Adds desired line to the current LineDatabase object.

        **Parameters**

        element : string
            The element of the transition using element's symbol on periodic table

        ion_state : string
            The roman numeral representing the ionic state of the transition

        wavelength : float
            The wavelength of the transition in angstroms

        gamma : float
            The gamma of the transition in Hertz

        f_value: float
            The oscillator strength of the transition

        field : string, optional
            The default yt field name associated with the ion responsible for
            this line
            Default: None

        identifier : string, optional
            An optional identifier for the transition
            Default: None
        """
        self.line_database.add_line(element, ion_state, wavelength,
                                    gamma, f_value, field=field,
                                    identifier=identifier)

    def save_spectrum(self, filename='spectrum.h5', format=None):
        """
        Save the current spectral data to an output file.  Unless specified,
        the output data format will be determined by the suffix of the filename
        provided ("h5":HDF5, "fits":FITS, all other:ASCII).

        """
        if format is None:
            if filename.endswith('.h5'):
                self._write_spectrum_hdf5(filename)
            elif filename.endswith('.fits'):
                self._write_spectrum_fits(filename)
            else:
                self._write_spectrum_ascii(filename)
        elif format == 'HDF5':
            self._write_spectrum_hdf5(filename)
        elif format == 'FITS':
            self._write_spectrum_fits(filename)
        elif format == 'ASCII':
            self._write_spectrum_ascii(filename)
        else:
            mylog.warn("Invalid format.  Must be 'HDF5', 'FITS', 'ASCII'. Defaulting to ASCII.")
            self._write_spectrum_ascii(filename)

    def plot_spectrum(self, filename="spectrum.png",
                      lambda_limits=None, flux_limits=None, step=False,
                      title=None, label=None, figsize=None, features=None):
        """
        Plot the spectrum from the SpectrumGenerator class.

        This is a convenience method that wraps the plot_spectrum standalone
        function for use with the data from the SpectrumGenerator itself.

        **Parameters**

        filename : string, optional
            Output filename of the plotted spectrum.  Will be a png file.

        lambda_limits : tuple or list of floats, optional
            The minimum and maximum of the lambda range (x-axis) for the plot
            in angstroms.  If specified as None, will use whole lambda range 
            of spectrum.
            Default: None

        flux_limits : tuple or list of floats, optional
            The minimum and maximum of the flux range (y-axis) for the plot.
            If specified as None, limits are automatically from 
            [0, 1.1*max(flux)].
            Default: None

        title : string, optional
            Title for plot

        label : string or list of strings, optional
            Label for each spectrum to be plotted
        """
        plot_spectrum(self.lambda_field, self.flux_field, filename=filename,
                      lambda_limits=lambda_limits, flux_limits=flux_limits,
                      step=step, title=title, figsize=figsize, features=features)

    def __repr__(self):
        disp = "<SpectrumGenerator>:\n"
        disp += "    lambda_field: %s\n" % self.lambda_field
        disp += "    tau_field: %s\n" % self.tau_field
        disp += "    flux_field: %s\n" % self.flux_field
        disp += "%s" % self.instrument
        return disp
