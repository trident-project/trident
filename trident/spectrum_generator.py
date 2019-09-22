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

import numpy as np
import os

from trident.absorption_spectrum.absorption_spectrum import \
    AbsorptionSpectrum
from yt.convenience import \
    load
from yt.data_objects.data_containers import \
    YTDataContainer
from yt.data_objects.static_output import \
    Dataset
from yt.funcs import \
    mylog, \
    YTArray

from trident.config import \
    ion_table_dir, \
    ion_table_file, \
    ion_table_filepath
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
from trident.config import \
    trident_path
from yt.utilities.on_demand_imports import \
    _h5py, \
    _astropy

# Valid instruments
valid_instruments = \
    {'COS' :
       Instrument(1150, 1450, dlambda=0.01, lsf_kernel='avg_COS_G130M.txt', name='COS'),
     'COS-G130M' :
       Instrument(1150, 1450, dlambda=0.01, lsf_kernel='avg_COS_G130M.txt', name='COS-G130M'),
     'COS-G140L' :
       Instrument(1130, 2330, dlambda=0.08, lsf_kernel='avg_COS_G140L.txt', name='COS-G140L'),
     'COS-G160M' :
       Instrument(1405, 1777, dlambda=0.012, lsf_kernel='avg_COS_G160M.txt', name='COS-G160M')
     }

class SpectrumGenerator(AbsorptionSpectrum):
    """
    Preferred class for generating, storing, and plotting absorption-line spectra.
    SpectrumGenerator is a subclass of yt's AbsorptionSpectrum class
    with additional functionality like line lists, post-processing to
    include Milky Way foreground, quasar backgrounds, applying line-spread
    functions, and plotting.

    User first specifies the telescope/instrument used for generating spectra
    (e.g. 'COS' for the Cosmic Origins Spectrograph aboard the
    Hubble Space Telescope).  This can be done by naming the
    :class:`~trident.Instrument`, or manually setting all of the spectral
    characteristics including ``lambda_min``, ``lambda_max``, ``lsf_kernel``,
    and ``n_lambda`` or ``dlambda``.  If none of these arguments are set,
    defaults to 'COS' as the default instrument covering 1150-1450 Angstroms
    with a binsize (``dlambda``) of 0.1 Angstroms.

    Once a :class:`~trident.SpectrumGenerator` has been initialized, pass it
    ``LightRay`` objects using :class:`~trident.SpectrumGenerator.make_spectrum`
    to actually generate the spectra themselves.  Then one can post-process,
    plot, or save them using
    :class:`~trident.SpectrumGenerator.add_milky_way_foreground`,
    :class:`~trident.SpectrumGenerator.add_qso_spectrum`,
    :class:`~trident.SpectrumGenerator.apply_lsf`,
    :class:`~trident.SpectrumGenerator.save_spectrum`, and
    :class:`~trident.SpectrumGenerator.plot_spectrum`.

    **Parameters**

    :instrument: string, optional

        The spectrograph to use.  Currently, the only named options are
        different observing modes of the Cosmic Origins Spectrograph 'COS':
        'COS-G130M', 'COS-G160M', and 'COS-G140L' as well as 'COS' which
        aliases to 'COS-G130M'. These automatically set the ``lambda_min``,
        ``lambda_max``, ``dlambda`` and ``lsf_kernel``s appropriately.
        If you're going to set ``lambda_min``, ``lambda_max``, et al manually,
        leave this set to None.
        Default: None

    :lambda_min: float, YTQuantity, or 'auto'

       lower wavelength bound in angstroms or velocity bound in km/s
       (if bin_space set to 'velocity'). If set to 'auto', the lower
       bound will be automatically adjusted to encompass all absorption
       lines. The window will not be expanded for continuum features,
       only absorption lines.

    :lambda_max: float, YTQuantity, or 'auto'

       upper wavelength bound in angstroms or velocity bound in km/s
       (if bin_space set to 'velocity'). If set to 'auto', the upper
       bound will be automatically adjusted to encompass all absorption
       lines. The window will not be expanded for continuum features,
       only absorption lines.

    :n_lambda: int

        The number of bins in the spectrum (inclusive), so if
        extrema = 10 and 20, and dlambda (binsize) = 1, then n_lambda = 11.
        Default: None

    :dlambda: float

        size of the wavelength bins in angstroms or velocity bins in km/s.
        Default: None

    :bin_space: 'wavelength' or 'velocity'

        Sets the dimension in which spectra are created. If set to
        wavelength, the resulting spectra are flux (or tau) vs. observed
        wavelength. If set to velocity, the spectra are flux vs.
        velocity offset from the rest wavelength of the absorption line.
        Default: wavelength

    :lsf_kernel: string, optional

        The filename for the LSF kernel. Files are found in
        trident.__path__/data/lsf_kernels or current working directory.
        Only necessary if you are applying an LSF to the spectrum in
        postprocessing.
        Default: None

    :line_database: string or :class:`~trident.LineDatabase`, optional

        A text file listing the various lines to insert into the line database,
        or a :class:`~trident.LineDatabase` object in memory. The line database
        provides a list of all possible lines that could be added to the
        spectrum. For a text file, it should have 4 tab-delimited columns of
        name (e.g. MgII), wavelength in angstroms, gamma of transition, and
        f-value of transition. See example datasets in
        trident.path/data/line_lists for examples.
        Default: lines.txt

    :ionization_table: hdf5 file, optional

        An HDF5 file used for computing the ionization fraction of the gas
        based on its density, temperature, metallicity, and redshift.  If
        set to None, will use the ion table defined in your Trident config
        file.
        Default: None

    **Example**

    Create a one-zone ray, and generate a COS spectrum from that ray.

    >>> import trident
    >>> ray = trident.make_onezone_ray()
    >>> sg = trident.SpectrumGenerator('COS')
    >>> sg.make_spectrum(ray)
    >>> sg.plot_spectrum('spec_raw.png')

    Create a one-zone ray at redshift 0.5, and generate a spectrum with 1
    angstrom spectral bins from 2000-4000 angstroms, then post-process by
    adding a MW foreground a QSO background at z=0.5 and add a boxcar line
    spread function of 100 angstroms width.  Plot it and save the figure to
    'spec_final.png'.

    >>> import trident
    >>> ray = trident.make_onezone_ray(redshift=0.5)
    >>> sg = trident.SpectrumGenerator(lambda_min=2000, lambda_max=4000,
    ... dlambda=1)
    >>> sg.make_spectrum(ray)
    >>> sg.add_qso_spectrum(emitting_redshift=.5)
    >>> sg.add_milky_way_foreground()
    >>> sg.apply_lsf(function='boxcar', width=100)
    >>> sg.plot_spectrum('spec_final.png')
    """
    def __init__(self, instrument=None, lambda_min=None, lambda_max=None,
                 n_lambda=None, dlambda=None, lsf_kernel=None,
                 line_database='lines.txt', ionization_table=None,
                 bin_space='wavelength'):
        if instrument is None and \
          ((lambda_min is None or lambda_max is None) or \
           (dlambda is None and n_lambda is None)):
            instrument = 'COS'
            mylog.info("No parameters specified, defaulting to COS instrument.")
        elif instrument is None:
            instrument = Instrument(lambda_min=lambda_min,
                                    lambda_max=lambda_max,
                                    n_lambda=n_lambda,
                                    dlambda=dlambda,
                                    lsf_kernel=lsf_kernel, name="Custom",
                                    bin_space=bin_space)
        self.observing_redshift = 0.
        self._set_instrument(instrument)
        mylog.info("Setting instrument to %s" % self.instrument.name)
        self.dlambda = self.instrument.dlambda

        AbsorptionSpectrum.__init__(self,
                                    self.instrument.lambda_min,
                                    self.instrument.lambda_max,
                                    n_lambda=self.instrument.n_lambda,
                                    dlambda=self.instrument.dlambda,
                                    bin_space=bin_space)

        if isinstance(line_database, LineDatabase):
            self.line_database = line_database
        else:
            # instantiate the LineDatabase
            self.line_database = LineDatabase(line_database)

        # Instantiate the spectrum to be zeros and ones for tau_field and
        # flux_field respectively.
        self.clear_spectrum()

        if ionization_table is not None:
            # figure out where the user-specified files lives
            if os.path.isfile(ion_table_file):
                self.ionization_table = ion_table_file
            elif os.path.isfile(ion_table_filepath):
                self.ionization_table = ion_table_filepath
            else:
                raise RuntimeError("ionization_table %s is not found in local "
                                   "directory or in %s" %
                                   (ion_table_file.split(os.sep)[-1],
                                    ion_table_dir))
        else:
            self.ionization_table = None

    def make_spectrum(self, ray, lines='all',
                      output_file=None,
                      output_absorbers_file=None,
                      use_peculiar_velocity=True,
                      observing_redshift=0.0,
                      ly_continuum=True,
                      store_observables=False,
                      min_tau=1e-3,
                      njobs="auto"):
        """
        Make a spectrum from ray data depositing the desired lines.  Make sure
        to pass this function a LightRay object and potentially also a list of
        strings representing what lines you'd like to actually have be
        deposited in your final spectrum.

        **Parameters**

        :ray: string, dataset, or data container

            If a string, the path to the ray dataset. As a dataset,
            this is the ray dataset loaded by yt. As a data container,
            this is a data object created from a ray dataset, such as
            a cut region.

        :lines: list of strings

            List of strings that determine which lines will be added
            to the spectrum.  List can include things like "C", "O VI",
            or "Mg II ####", where #### would be the integer wavelength
            value of the desired line.  If set to 'all', includes all lines
            in LineDatabase set in SpectrumGenerator.
            Default: 'all'

        :output_file: optional, string

            Filename of output if you wish to save the spectrum immediately
            without any further processing. File formats are chosen based on the
            filename extension.  ".h5" for HDF5, ".fits" for FITS,
            and everything else is ASCII.  Equivalent of calling
            :class:`~trident.SpectrumGenerator.save_spectrum`.
            Default: None

        :output_absorbers_file: optional, string

           Option to save a text file containing all of the absorbers and
           corresponding wavelength and redshift information.
           For parallel jobs, combining the lines lists can be slow so it
           is recommended to set to None in such circumstances.
           Default: None

        :use_peculiar_velocity: optional, bool

            If True, include the effects of doppler redshift of the gas
            in shifting lines in the final spectrum.
            Default: True

        :observing_redshift: optional, float

            This is the value of the redshift at which the observer of this
            spectrum exists.  In most cases, this will be a redshift of 0.
            Default: 0.

        :ly_continuum: optional, boolean

            If any H I lines are used in the line list, this assures a
            Lyman continuum will be included in the spectral generation.
            Lyman continuum begins at final Lyman line deposited (Ly 39 =
            912.32 A) not at formal Lyman Limit (911.76 A) so as to not have
            a gap between final Lyman lines and continuum.  Uses power law
            of index 3 and normalization to match opacity of final Lyman lines.
            Default: True

        :store_observables: optional, boolean

            If set to true, observable properties for each cell in the light
            ray will be saved for each line in the line list. Properties
            include the column density, tau, thermal b, and the wavelength
            where tau was deposited. Best applied for a reasonable number
            of lines.  These quantities will be saved to the SpectrumGenerator
            attribute: 'line_observables_dict'.
            Default: False

        :min_tau: optional, float
           This value determines size of the wavelength window used to
           deposit lines or continua.  The wavelength window is expanded
           until the optical depth at the edge is below this value.  If too
           high, this will result in features appearing cut off at the edges.
           Decreasing this will make features smoother but will also increase
           run time.  An increase by a factor of ten will result in roughly a
           2x slow down.
           Default: 1e-3.

        :njobs: optional, int or "auto"

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

        **Example**

        Make a one zone ray and generate a COS spectrum for it including
        only Oxygen VI, Mg II, and all Carbon lines, and plot to disk.

        >>> import trident
        >>> ray = trident.make_onezone_ray()
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray, lines=['O VI', 'Mg II', 'C'])
        >>> sg.plot_spectrum('spec_raw.png')
        """
        self.observing_redshift = observing_redshift

        if isinstance(ray, str):
            ray = load(ray)
        if isinstance(ray, Dataset):
            ad = ray.all_data()
        elif isinstance(ray, YTDataContainer):
            ad = ray
            ray = ad.ds
        else:
            raise RuntimeError("Unrecognized ray type.")

        # Clear out any previous spectrum that existed first
        self.clear_spectrum()

        active_lines = self.line_database.parse_subset(lines)

        # Make sure we've produced all the necessary
        # derived fields if they aren't native to the data
        for line in active_lines:
            # if successful, means line.field is in ds.derived_field_list
            try:
                ad._determine_fields(line.field)[0]
            # otherwise we probably need to add the field to the dataset
            except BaseException:
                my_ion = \
                  line.field[:line.field.find("number_density")]
                on_ion = my_ion.split("_")
                # Add the field if greater than level 1 ionization
                # because there is only one naming convention for these fields:
                # X_pY_number_density
                if on_ion[1]:
                    my_lev = int(on_ion[1][1:]) + 1
                    mylog.info("Creating %s from ray's density, "
                               "temperature, metallicity." % (line.field))
                    add_ion_number_density_field(on_ion[0], my_lev, ray,
                                     ionization_table=self.ionization_table)
                # If level 1 ionization, check to see if other name for
                # field is present in dataset
                else:
                    my_lev = 1
                    alias_field = ('gas', "".join([my_ion, 'p0_number_density']))
                    # Don't add the X_number_density if X_p0_number_density is
                    # in dataset already
                    if alias_field in ray.derived_field_list:
                        line.field = alias_field
                    # But add the field if neither X_number_density nor
                    # X_p0_number_density is in the dataset
                    else:
                        mylog.info("Creating %s from ray's density, "
                                   "temperature, metallicity." % (line.field))
                        add_ion_number_density_field(on_ion[0], my_lev, ray,
                                     ionization_table=self.ionization_table)

            self.add_line(line.identifier, line.field,
                          float(line.wavelength),
                          float(line.f_value),
                          float(line.gamma),
                          atomic_mass[line.element],
                          label_threshold=1e3)

        # If there are H I lines present, add a Lyman continuum source
        # Lyman continuum source starts at wavelength where last Lyman line
        # is deposited (Ly 40), as opposed to true Lyman Limit at 911.763 A
        # so there won't be a gap between lines and continuum.  Using
        # power law of index 3.0 and normalization to match the opacity of
        # the final Lyman line into the FUV.
        H_lines = self.line_database.select_lines(source_list=active_lines,
                                                  element='H', ion_state='I')
        if len(H_lines) > 0 and ly_continuum:
            self.add_continuum('Ly C', H_lines[0].field, 912.32336, 1.6e17, 3.0)

        AbsorptionSpectrum.make_spectrum(self, ad,
                                         output_file=None,
                                         line_list_file=None,
                                         output_absorbers_file=output_absorbers_file,
                                         use_peculiar_velocity=use_peculiar_velocity,
                                         observing_redshift=observing_redshift,
                                         store_observables=store_observables,
                                         min_tau=min_tau, njobs=njobs)

    def _get_qso_spectrum(self, emitting_redshift, observing_redshift,
                          filename=None):
        """
        Read in QSO spectrum and return an interpolated version. Interpolated
        version for fitting the desired wavelength interval and binning.
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
        Read in MW spectrum and return an interpolated version.
        Interpolated version for fitting the desired wavelength interval and
        binning.

        Source data come from Charles Danforth and are a median-combination of
        92 normalized COS/G130M+G160M AGN spectra valid from 1130-1800A.
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
        Postprocess a spectrum to add a Milky Way foreground.  Data
        from Charles Danforth. Median-filter of 92 normalized
        COS/G130M+G160M AGN spectra spanning the wavelength range of
        1130 to 1800 Angstroms in 0.07 Angstrom bin size.

        **Parameters**

        :flux_field: optional, array

            Array of flux values to which the Milky Way foreground is applied.
            Default: None

        :filename: string

            Filename where the Milky Way foreground values used to modify
            the flux are stored.
            Default: None

        **Example**

        Make a one zone ray and generate a COS spectrum for it.  Add
        MW foreground to it, and save it.

        >>> import trident
        >>> ray = trident.make_onezone_ray()
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray)
        >>> sg.add_milky_way_foreground()
        >>> sg.plot_spectrum('spec_mw_corrected.png')

        Plot a naked MW spectrum.

        >>> import trident
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.add_milky_way_foreground()
        >>> sg.plot_spectrum('spec_mw.png')
        """
        if flux_field is None:
            flux_field = self.flux_field
        MW_spectrum = self._get_milky_way_foreground(filename=filename)
        flux_field *= MW_spectrum

        # Negative fluxes don't make sense, so clip
        np.clip(flux_field, 0, np.inf, out=flux_field)

    def add_qso_spectrum(self, flux_field=None,
                         emitting_redshift=None,
                         observing_redshift=None,
                         filename=None):
        """
        Postprocess a spectrum to add a QSO spectrum background. Uses data from
        Telfer et al., ApJ, 565, 773 "The Rest-Frame Extreme Ultraviolet
        Spectral Properties of QSO". HST Radio Quiet composite for < 1275 Ang,
        SDSS composite > 2000 Ang, mean in between 8251 0

        **Parameters**

        :flux_field: array, optional

            Array of flux values to which the quasar background is applied.
            Default: None

        :emitting_redshift: float, optional

            Redshift value at which the QSO emitted its light.  If specified
            as None, use 0.
            Default: None

        :observing_redshift: float, optional

            Redshift value at which the quasar is observed.  If specified as
            None, use the observing_redshift value specified in make_spectrum()
            which defaults to 0.
            Default: None

        :filename: string

            Filename where the Milky Way foreground values used to modify
            the flux are stored.
            Default: None

        **Example**

        Make a one zone ray at redshift of .5 and generate a COS spectrum for
        it.  Add z=0.5 quasar background to it, and save it.

        >>> import trident
        >>> ray = trident.make_onezone_ray(redshift=0.5)
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray)
        >>> sg.add_qso_spectrum(emitting_redshift=0.5)
        >>> sg.plot_spectrum('spec_qso_corrected.png')

        Plot a naked QSO spectrum at z=.1

        >>> import trident
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.add_qso_spectrum(emitting_redshift=.1)
        >>> sg.plot_spectrum('spec_qso.png')
         """
        if flux_field is None:
            flux_field = self.flux_field
        qso_spectrum = self._get_qso_spectrum(emitting_redshift=emitting_redshift,
                                              observing_redshift=observing_redshift,
                                              filename=filename)
        flux_field *= qso_spectrum

        # Negative fluxes don't make sense, so clip
        np.clip(flux_field, 0, np.inf, out=flux_field)

    def add_gaussian_noise(self, snr, seed=None):
        """
        Postprocess a spectrum to add gaussian random noise of a given SNR.

        **Parameters**

        :snr: float

            The desired signal-to-noise ratio for determining the amount of
            gaussian noise

        :seed: optional, int

            Seed for the random number generator.  This should be used to
            ensure than the same noise is added each time the spectrum is
            regenerated, if desired.
            Default: None

        **Example**

        Make a one zone ray and generate a COS spectrum for it.  Add noise
        to the spectrum as though it were observed with a signal to noise
        ratio of 30.

        >>> import trident
        >>> ray = trident.make_onezone_ray(redshift=0.5)
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray)
        >>> sg.add_gaussian_noise(30)
        >>> sg.plot_spectrum('spec_noise_corrected.png')

        Plot a DLA with SNR of 10.

        >>> import trident
        >>> ray = trident.make_onezone_ray(column_densities={'H_p0_number_density':1e21})
        >>> sg = trident.SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.1)
        >>> sg.make_spectrum(ray, lines=['Ly a'])
        >>> sg.add_gaussian_noise(10)
        >>> sg.plot_spectrum('spec_noise.png')
        """
        self.snr = snr
        np.random.seed(seed)
        noise = np.random.normal(loc=0.0, scale=1/float(snr),
                                 size=self.flux_field.size)
        self.add_noise_vector(noise)

        # Negative fluxes don't make sense, so clip
        np.clip(self.flux_field, 0, np.inf, out=self.flux_field)

    def add_noise_vector(self, noise):
        """
        Add an array of noise to the spectrum.

        **Parameters**

        :noise: array of floats

            The array of noise values to be added to the spectrum.  This
            array must be of the same size as the flux array.

        **Example**

        >>> import numpy as np
        >>> import trident
        >>> ray = trident.make_onezone_ray(column_densities={'H_p0_number_density':1e21})
        >>> sg = trident.SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.1)
        >>> sg.make_spectrum(ray, lines=['Ly a'])
        >>> my_noise = np.random.normal(loc=0.0, scale=0.1, size=sg.flux_field.size)
        >>> sg.add_noise_vector(my_noise)
        >>> sg.plot_spectrum('spec_noise.png')
        """

        if not isinstance(noise, np.ndarray):
            raise SyntaxError(
                "Noise field must be a numpy array.")
        if self.flux_field.shape != noise.shape:
            raise SyntaxError(
                "Flux (%s) and noise (%s) vectors must have same shape." %
                (self.flux_field.shape, noise.shape))
        self.flux_field += noise
        self.snr = 1 / np.std(noise)

    def apply_lsf(self, function=None, width=None, filename=None):
        """
        Postprocess a spectrum to apply a line spread function.
        If the SpectrumGenerator already has an LSF_kernel set, it will
        be used when no keywords are supplied.  Otherwise, the user can
        specify a filename of a user-defined kernel or a function+width
        for a kernel.  Valid functions are: "boxcar" and "gaussian".

        For more information, see :class:`~trident.LSF` and
        :class:`~trident.Instrument`.

        **Parameters**

        :function: string, optional

            Desired functional form for the applied LSF kernel.
            Valid options are currently "boxcar" or "gaussian"
            Default: None

        :width: int, optional

            Width of the desired LSF kernel in bin elements
            Default: None

        :filename: string, optional

            The filename of the user-supplied kernel for applying the LSF
            Default: None

        **Example**

        Make a one zone ray and generate a COS spectrum for it.  Apply the
        COS line spread function to it.

        >>> import trident
        >>> ray = trident.make_onezone_ray()
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray)
        >>> sg.apply_lsf()
        >>> sg.plot_spectrum('spec_lsf_corrected.png')

        Make a one zone ray and generate a spectrum with bin width = 1 angstrom.
        Apply a boxcar LSF to it with width 50 angstroms.

        >>> import trident
        >>> ray = trident.make_onezone_ray()
        >>> sg = trident.SpectrumGenerator(lambda_min=1100, lambda_max=1200, dlambda=1)
        >>> sg.make_spectrum(ray)
        >>> sg.apply_lsf(function='boxcar', width=50)
        >>> sg.plot_spectrum('spec_lsf_corrected.png')
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

        # Negative fluxes don't make sense, so clip
        np.clip(self.flux_field, 0, np.inf, out=self.flux_field)

    def load_spectrum(self, lambda_field=None, tau_field=None, flux_field=None):
        """
        Load data arrays into an existing spectrum object.

        **Parameters**

        :lambda_field: array

            The array of valid wavelengths
            Default: None

        :tau_field: array

            The array of optical depths for the corresponding wavelengths
            Default: None

        :flux_field: array

            The array of flux values for the corresponding wavelengths
            Default: None

        **Example**

        Loading a custom set of data into an existing SpectrumGenerator object:

        >>> import trident
        >>> import numpy as np
        >>> lambda_field = np.arange(1200,1400)
        >>> flux_field = np.ones(len(lambda_field))
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.load_spectrum(lambda_field=lambda_field, flux_field=flux_field)
        >>> sg.plot_spectrum('temp.png')
        """
        if lambda_field is not None:
            self.lambda_field = lambda_field
            self.lambda_min = self.lambda_field[0]
            self.lambda_max = self.lambda_field[-1]
            self.n_lambda = len(self.lambda_field)
            self.dlambda = self.lambda_field[1] - self.lambda_field[0]
        if tau_field is not None:
            self.tau_field = tau_field
        if flux_field is not None:
            self.flux_field = flux_field
        if not len(self.flux_field) == len(self.lambda_field):
            raise RuntimeError("Loaded spectra must have the same dimensions"
                               "for lambda_field and flux_field.  Currently:"
                               "len(lambda_field) = %d, len(flux_field) = %d"
                               % (self.lambda_field, self.flux_field))

    def clear_spectrum(self):
        """
        Clear the current spectrum in the SpectrumGenerator.
        Clears the existing spectrum's flux and tau fields and replaces them
        with ones and zeros respectively.  Clear the line list kept in
        the AbsorptionSpectrum object as well.  Also clear the line_subset
        stored by the LineDatabase.
        """
        if self.lambda_field is not None:
            # Set flux and tau to ones and zeros
            self.flux_field = np.ones(self.lambda_field.size)
            self.tau_field = np.zeros(self.lambda_field.size)
        else:
            self.flux_field = None
            self.tau_field = None

        # Clear out the line list that is stored in AbsorptionSpectrum
        self.line_list = []

        # Make sure we reset the line database as well
        self.line_database.lines_subset = []

    def _set_instrument(self, instrument):
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
                raise RuntimeError("_set_instrument accepts only Instrument "
                                   "objects or the names of valid "
                                   "instruments: ", valid_instruments.keys())
            self.instrument = valid_instruments[instrument]
        elif isinstance(instrument, Instrument):
            self.instrument = instrument
        else:
            raise RuntimeError("_set_instrument accepts only Instrument "
                               "objects or the names of valid instruments: ",
                               valid_instruments.keys())

    def add_line_to_database(self, element, ion_state, wavelength, gamma,
                             f_value, field=None, identifier=None):
        """
        Adds desired line to the current :class:`~trident.LineDatabase` object.

        **Parameters**

        :element: string

            The element of the transition using element's symbol on periodic table

        :ion_state: string

            The roman numeral representing the ionic state of the transition

        :wavelength: float

            The wavelength of the transition in angstroms

        :gamma: float

            The gamma of the transition in Hertz

        :f_value: float

            The oscillator strength of the transition

        :field: string, optional

            The default yt field name associated with the ion responsible for
            this line
            Default: None

        :identifier: string, optional

            An optional identifier for the transition
            Default: None
        """
        self.line_database.add_line(element, ion_state, wavelength,
                                    gamma, f_value, field=field,
                                    identifier=identifier)

    def save_spectrum(self, filename='spectrum.h5', format=None):
        """
        Save the current spectrum data to an output file.  Unless specified,
        the output data format will be determined by the suffix of the filename
        provided ("h5":HDF5, "fits":FITS, all other:ASCII).

        ASCII data is stored as a tab-delimited text file.

        **Parameters**

        :filename: string, optional

            Output filename for storing the data.
            Default: 'spectrum.h5'

        :format: string, optional

            Data format of the output file.  Valid examples are "HDF5",
            "FITS", and "ASCII".  If None is set, selects based on suffix
            of filename.
            Default: None

        **Example**

        Save a spectrum to disk, load it from disk, and plot it.

        >>> import trident
        >>> ray = trident.make_onezone_ray()
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray)
        >>> sg.save_spectrum('temp.h5')
        >>> sg.clear_spectrum()
        >>> sg.load_spectrum('temp.h5')
        >>> sg.plot_spectrum('temp.png')
        """
        if format is None:
            if filename.endswith('.h5') or filename.endswith('hdf5'):
                self._write_spectrum_hdf5(filename)
            elif filename.endswith('.fits') or filename.endswith('FITS'):
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
            mylog.warning("Invalid format.  Must be 'HDF5', 'FITS', 'ASCII'. Defaulting to ASCII.")
            self._write_spectrum_ascii(filename)

    def plot_spectrum(self, filename="spectrum.png",
                      lambda_limits=None, flux_limits=None, step=False,
                      title=None, label=None, figsize=None, features=None,
                      axis_labels=None):
        """
        Plot the current spectrum and save to disk.

        This is a convenience method that wraps the
        :class:`~trident.plot_spectrum` standalone function for use with the
        data from the :class:`~trident.SpectrumGenerator` itself.

        **Parameters**

        :filename: string, optional

            Output filename of the plotted spectrum.  Will be a png file.
            Default: 'spectrum.png'

        :lambda_limits: tuple or list of floats, optional

            The minimum and maximum of the lambda range (x-axis) for the plot
            in angstroms.  If specified as None, will use whole lambda range
            of spectrum.  Example: (1200, 1400) for 1200-1400 Angstroms.
            Default: None

        :flux_limits: tuple or list of floats, optional

            The minimum and maximum of the flux range (y-axis) for the plot.
            If specified as None, limits are automatically from
            [0, 1.1*max(flux)].  Example: (0, 1) for normal flux range before
            postprocessing.
            Default: None

        :step: boolean, optional

            Plot the spectrum as a series of step functions.  Appropriate for
            plotting processed and noisy data.

        :title: string, optional

            Optional title for plot
            Default: None

        :label: string, optional

            Label for spectrum to be plotted.  Will automatically trigger a
            legend to be generated.
            Default: None

        :features: dict, optional

            Include vertical lines with labels to represent certain spectral
            features.  Each entry in the dictionary consists of a key string to
            be overplot and the value float as to where in wavelength space it
            will be plot as a vertical line with the corresponding label.

            Example: features={'Ly a' : 1216, 'Ly b' : 1026}

            Default: None

        :axis_labels: tuple of strings, optional

            Optionally set the axis labels directly.  If set to None, defaults to
            ('Wavelength [$\\rm\\AA$]', 'Relative Flux').
            Default: None

        **Example**

        Create a one-zone ray, and generate a COS spectrum from that ray. Plot
        the resulting spectrum highlighting the Lyman alpha feature.

        >>> import trident
        >>> ray = trident.make_onezone_ray()
        >>> sg = trident.SpectrumGenerator('COS')
        >>> sg.make_spectrum(ray)
        >>> sg.plot_spectrum('spec_raw.png', features={'Ly a' : 1216})
        """

        if self.tau_field is None:
            mylog.warning('Spectrum is totally empty, no plotting to be done.')
            return

        plot_spectrum(self.lambda_field, self.flux_field, filename=filename,
                      lambda_limits=lambda_limits, flux_limits=flux_limits,
                      title=title, label=label, figsize=figsize, step=step,
                      features=features, axis_labels=axis_labels)

    def __repr__(self):
        disp = "<SpectrumGenerator>:\n"
        disp += "    lambda_field: %s\n" % self.lambda_field
        disp += "    tau_field: %s\n" % self.tau_field
        disp += "    flux_field: %s\n" % self.flux_field
        disp += "%s" % self.instrument
        return disp

def load_spectrum(filename, format='auto', instrument=None, lsf_kernel=None,
                  line_database='lines.txt', ionization_table=None):
    """
    Load a previously saved spectrum from disk.

    **Parameters**

    :filename: string

        Filename of the saved spectrum.

    :format: string

        File format of the saved spectrum file.  Valid values are: "auto",
        "hdf5", "fits", and "ascii".  If you select "auto", the code will
        attempt to auto-detect the file format from the extension of the data
        file: ".h5" or ".hdf5" -> hdf5, ".fits" or ".FITS" -> fits, all other
        -> ascii.
        Default: "auto"

    :instrument: string, optional

        The telescope+instrument combination to use for the loaded spectrum.
        Default: None

    :lsf_kernel: string, optional

        The filename for the LSF kernel to use for the loaded spectrum.
        Default: None

    :line_database: string, optional

        A text file listing the various lines to insert into the line database
        to use for the loaded spectrum.
        Default: None

    :ionization_table: hdf5 file, optional

        An HDF5 file used for computing the ionization fraction of the gas
        based on its density, temperature, metallicity, and redshift.
        Default: None

    **Example**

    Create a simple spectrum, save it to disk, and load it back as a new
    SpectrumGenerator object.

    >>> import trident
    >>> ray = trident.make_onezone_ray()
    >>> sg = trident.SpectrumGenerator('COS')
    >>> sg.make_spectrum(ray)
    >>> sg.save_spectrum('spec.h5')
    >>> sg_copy = trident.load_spectrum('spec.h5')
    """
    if format == 'auto':
        if filename.endswith('.h5') or filename.endswith('.hdf5'):
            format = 'hdf5'
        elif filename.endswith('.fits') or filename.endswith('.FITS'):
            format = 'fits'
        else:
            format = 'ascii'
    if format == 'hdf5':
        f = _h5py.File(filename, 'r')
        lambda_field = f['wavelength'][()]
        flux_field = f['flux'][()]
        tau_field = f['tau'][()]
    elif format == 'fits':
        pyfits = _astropy.pyfits
        hdulist = pyfits.open(filename)
        data = hdulist[1].data
        lambda_field = data['wavelength']
        tau_field = None
        # Switch above line to tau_field = data['tau'] when yt PR #2314 is merged.
        flux_field = data['flux']
    elif format == 'ascii':
        data = np.genfromtxt(filename)
        lambda_field = data[:,0]
        tau_field = data[:,1]
        flux_field = data[:,2]
    else:
        raise RuntimeError("load_spectrum 'format' keyword must be 'hdf5', 'ascii', 'fits', or 'auto'")

    lambda_field = YTArray(lambda_field, "angstrom")
    lambda_min = lambda_field[0]
    lambda_max = lambda_field[-1]
    n_lambda = lambda_field.size
    sg = SpectrumGenerator(instrument=instrument, lambda_min=lambda_min,
                           lambda_max=lambda_max, n_lambda=n_lambda,
                           lsf_kernel=lsf_kernel, line_database=line_database,
                           ionization_table=ionization_table)
    if tau_field is not None:
        sg.load_spectrum(lambda_field=lambda_field, tau_field=tau_field,
                         flux_field=flux_field)
    else:
        sg.load_spectrum(lambda_field=lambda_field, flux_field=flux_field)
    return sg
