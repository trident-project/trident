"""
AbsorptionSpectrum class and member functions.



"""

from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2017, yt Development Team.
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

from yt.data_objects.data_containers import \
    YTDataContainer
from yt.data_objects.static_output import \
    Dataset
from yt.convenience import load
from yt.extern.six import string_types
from yt.funcs import get_pbar, mylog
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.on_demand_imports import _astropy
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    _get_comm, \
    parallel_objects, \
    parallel_root_only
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    speed_of_light_cgs

from trident.absorption_spectrum.absorption_line import \
    tau_profile

pyfits = _astropy.pyfits

_bin_space_units = {'wavelength': 'angstrom',
                    'velocity': 'km/s'}
c_kms = speed_of_light_cgs.to('km/s')

class AbsorptionSpectrum(object):
    r"""Base class for generating absorption spectra.  This code was originally
    based in yt and more restrictive in terms of what development was allowed,
    so the :class:`~trident.SpectrumGenerator` subclass has more advanced
    functionality built on top of this.  The base algorithm and functionality
    for spectral generation occurs here though.
    
    .. note::

        The preferred method for generating spectra is using
        :class:`~trident.SpectrumGenerator`.

    **Parameters**

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

    :n_lambda: optional, int

       number of bins. This cannot be set when setting either lambda_min
       or lambda_max to auto.

    :dlambda: optional, float or YTQuantity

      size of the wavelength bins in angstroms or velocity bins in km/s.

    :bin_space: 'wavelength' or 'velocity'

        Sets the dimension in which spectra are created. If set to
        wavelength, the resulting spectra are flux (or tau) vs. observed
        wavelength. If set to velocity, the spectra are flux vs.
        velocity offset from the rest wavelength of the absorption line.
        Default: wavelength
    """

    def __init__(self, lambda_min, lambda_max, n_lambda=None, dlambda=None,
                 bin_space='wavelength'):

        if bin_space not in _bin_space_units:
            raise RuntimeError(
                'Invalid bin_space value: "%s". Valid values are: "%s".' %
                (bin_space, '", "'.join(list(_bin_space_units))))
        self.bin_space = bin_space
        lunits = _bin_space_units[self.bin_space]

        if dlambda is not None:
            if not isinstance(dlambda, YTQuantity):
                dlambda = YTQuantity(dlambda, lunits)
            self.bin_width = dlambda

        if n_lambda is None and dlambda is None:
            raise RuntimeError(
                'Either n_lambda or dlambda must be given.')
        if n_lambda is not None and dlambda is not None:
            n_lambda = None

        if str(lambda_min) != 'auto':
            if not isinstance(lambda_min, YTQuantity):
                lambda_min = YTQuantity(
                    lambda_min, lunits)
            # round limits to bin size
            if dlambda is not None:
                lambda_min = np.round(lambda_min / dlambda) * dlambda
        self.lambda_min = lambda_min

        if str(lambda_max) != 'auto':
            if not isinstance(lambda_max, YTQuantity):
                lambda_max = YTQuantity(
                    lambda_max, lunits)
            # round limits to bin size
            if dlambda is not None:
                lambda_max = np.round(lambda_max / dlambda) * dlambda
        self.lambda_max = lambda_max

        self._auto_lambda = 'auto' in [str(self.lambda_min),
                                       str(self.lambda_max)]
        if self._auto_lambda and \
          (n_lambda is not None and n_lambda != 'auto'):
            raise RuntimeError(
                'Cannot set n_lambda when setting lambda_min or lambda_max to auto.')

        if dlambda is not None:
            self.bin_width = YTQuantity(dlambda, lunits)
            if not self._auto_lambda:
                n_lambda = \
                  self._get_field_size(self.lambda_min, self.lambda_max,
                                       self.bin_width)

        if self._auto_lambda:
            self.lambda_field = None
        else:
            if lambda_min >= lambda_max:
                raise RuntimeError(
                    'lambda_min (%f) cannot be greater than or equal to lambda_max (%f).' %
                    (lambda_min, lambda_max))

            if dlambda is None:
                n_lambda = int(n_lambda)
                self.bin_width = YTQuantity(
                    float(lambda_max - lambda_min) / (n_lambda - 1),
                    lunits)
            else:
                n_lambda = \
                  self._get_field_size(self.lambda_min, self.lambda_max,
                                       self.bin_width)

            self.lambda_field = \
              self._create_lambda_field(lambda_min, lambda_max, n_lambda)

        self.flux_field = None
        self.absorbers_list = None
        # a dictionary that will store spectral quantities for each index in the light ray
        self.line_observables_dict = None
        self.line_list = []
        self.continuum_list = []
        self.snr = 100  # default signal to noise ratio for error estimation

    def _get_field_size(self, lambda_min, lambda_max, dlambda):
        """
        Calculate number of bins.
        """

        if isinstance(lambda_min, YTQuantity):
            my_min = lambda_min.d
        else:
            my_min = lambda_min

        if isinstance(lambda_max, YTQuantity):
            my_max = lambda_max.d
        else:
            my_max = lambda_max

        if isinstance(dlambda, YTQuantity):
            my_dlambda = dlambda.d
        else:
            my_dlambda = dlambda

        n_lambda = (my_max - my_min) / my_dlambda + 1
        n_lambda = int(np.round(n_lambda))
        return n_lambda

    def _create_lambda_field(self, lambda_min, lambda_max, n_lambda,
                             units=None):
        """
        Create a lambda array with units.
        """

        if units is None:
            units = _bin_space_units[self.bin_space]

        if isinstance(lambda_min, YTQuantity):
            my_min = lambda_min.d
        else:
            my_min = lambda_min

        if isinstance(lambda_max, YTQuantity):
            my_max = lambda_max.d
        else:
            my_max = lambda_max

        return YTArray(np.linspace(my_min, my_max, n_lambda), units)

    _lambda_field = None
    @property
    def lambda_field(self):
        """
        The lambda field.
        """
        return self._lambda_field

    @lambda_field.setter
    def lambda_field(self, val):
        # This is useful for testing. If something goes wrong with
        # the lambda field, put something here to see every time
        # it gets set.
        self._lambda_field = val

    _tau_field = None
    @property
    def tau_field(self):
        """
        This is the total optical depth of all lines and continua.
        """
        if self.lambda_field is None:
            return None
        if self._tau_field is None:
            self._tau_field = np.zeros(self.lambda_field.size)
        return self._tau_field

    @tau_field.setter
    def tau_field(self, val):
        self._tau_field = val

    _current_tau_field = None
    @property
    def current_tau_field(self):
        """
        This is the optical depth array for the current absorption line
        being deposited. We will do the deposition of lines into this
        array, and then add it to self.tau_field at the end.
        """
        if self.lambda_field is None:
            return None
        if self._current_tau_field is None:
            self._current_tau_field = np.zeros(self.lambda_field.size)
        return self._current_tau_field

    @current_tau_field.setter
    def current_tau_field(self, val):
        self._current_tau_field = val

    def add_line(self, label, field_name, wavelength,
                 f_value, gamma, atomic_mass,
                 label_threshold=None):
        r"""Add an absorption line to the list of lines included in the spectrum.

        **Parameters**

        :label: string

           label for the line.

        :field_name: string

           field name from ray data for column densities.

        :wavelength: float

           line rest wavelength in angstroms.

        :f_value: float

           line f-value.

        :gamma: float

           line gamme value.

        :atomic_mass: float

           mass of atom in amu.
        """

        self.line_list.append({'label': label, 'field_name': field_name,
                               'wavelength': YTQuantity(wavelength, "angstrom"),
                               'f_value': f_value,
                               'gamma': gamma,
                               'atomic_mass': YTQuantity(atomic_mass, "amu"),
                               'label_threshold': label_threshold})

    def add_continuum(self, label, field_name, wavelength,
                      normalization, index):
        """
        Add a continuum feature that follows a power-law.

        **Parameters**

        :label: string

           label for the feature.

        :field_name: string

           field name from ray data for column densities.

        :wavelength: float

           line rest wavelength in angstroms.

        :normalization: float
        
           the column density normalization.

        :index: float

           the power-law index for the wavelength dependence.
        """

        self.continuum_list.append({'label': label, 'field_name': field_name,
                                    'wavelength': wavelength,
                                    'normalization': normalization,
                                    'index': index})

    def make_spectrum(self, input_object, output_file=None,
                      line_list_file=None, output_absorbers_file=None,
                      use_peculiar_velocity=True,
                      store_observables=False,
                      subgrid_resolution=10, observing_redshift=0.,
                      min_tau=1e-3, njobs="auto"):
        """
        Make spectrum from ray data using the line list.

        **Parameters**

        :input_object: string, dataset, or data container

           If a string, the path to the ray dataset. As a dataset,
           this is the ray dataset loaded by yt. As a data container,
           this is a data object created from a ray dataset, such as
           a cut region.

        :output_file: optional, string

           Option to save a file containing the wavelength, flux, and optical
           depth fields.  File formats are chosen based on the filename
           extension. ``.h5`` for hdf5, ``.fits`` for fits, and everything
           else is ASCII.
           Default: None

        :output_absorbers_file: optional, string

           Option to save a text file containing all of the absorbers and
           corresponding wavelength and redshift information.
           For parallel jobs, combining the lines lists can be slow so it
           is recommended to set to None in such circumstances.
           Default: None

        :use_peculiar_velocity: optional, bool
        
           if True, include peculiar velocity for calculating doppler redshift
           to shift lines.  Requires similar flag to be set in LightRay
           generation.
           Default: True

        :store_observables: optional, bool

           if True, stores observable properties of each cell along the line of
           sight for each line, such as tau, column density, and thermal b.
           These quantities will be saved to the AbsorptionSpectrum
           attribute: 'line_observables_dict'.
           Default: False

        :subgrid_resolution: optional, int

           When a line is being added that is unresolved (ie its thermal
           width is less than the spectral bin width), the voigt profile of
           the line is deposited into an array of virtual wavelength bins at
           higher resolution.  The optical depth from these virtual bins is
           integrated and then added to the coarser spectral wavelength bin.
           The subgrid_resolution value determines the ratio between the
           thermal width and the bin width of the virtual bins.  Increasing
           this value yields smaller virtual bins, which increases accuracy,
           but is more expensive.  A value of 10 yields accuracy to the 4th
           significant digit in tau.
           Default: 10

        :observing_redshift: optional, float

           This is the redshift at which the observer is observing
           the absorption spectrum.
           Default: 0

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
        self.snr = 100
        if line_list_file is not None:
            mylog.info("'line_list_file' keyword is deprecated. Please use " \
                       "'output_absorbers_file'.")
            output_absorbers_file = line_list_file

        input_fields = ['dl', 'redshift', 'temperature']
        field_units = {"dl": "cm", "redshift": "", "temperature": "K"}
        if use_peculiar_velocity:
            input_fields.append('velocity_los')
            input_fields.append('redshift_eff')
            field_units["velocity_los"] = "cm/s"
            field_units["redshift_eff"] = ""
        if observing_redshift != 0.:
            input_fields.append('redshift_dopp')
            field_units["redshift_dopp"] = ""
        for feature in self.line_list + self.continuum_list:
            if not feature['field_name'] in input_fields:
                input_fields.append(feature['field_name'])
                field_units[feature["field_name"]] = "cm**-3"

        if isinstance(input_object, string_types):
            input_ds = load(input_object)
            field_data = input_ds.all_data()
        elif isinstance(input_object, Dataset):
            input_ds = input_object
            field_data = input_ds.all_data()
        elif isinstance(input_object, YTDataContainer):
            input_ds = input_object.ds
            field_data = input_object

        # temperature field required to calculate voigt profile widths
        if ('temperature' not in input_ds.derived_field_list) and \
           (('gas', 'temperature') not in input_ds.derived_field_list):
            raise RuntimeError(
                "('gas', 'temperature') field required to be present in %s "
                "for AbsorptionSpectrum to function." % str(input_object))

        self.absorbers_list = []
        self.line_observables_dict = {}

        if njobs == "auto":
            comm = _get_comm(())
            njobs = min(comm.size, len(self.line_list))

        mylog.info("Creating spectrum")
        self._add_lines_to_spectrum(field_data, use_peculiar_velocity,
                                    output_absorbers_file,store_observables,
                                    subgrid_resolution=subgrid_resolution,
                                    observing_redshift=observing_redshift,
                                    min_tau=min_tau, njobs=njobs)
        self._add_continua_to_spectrum(field_data, use_peculiar_velocity,
                                       observing_redshift=observing_redshift,
                                       min_tau=min_tau)

        if self.tau_field is None:
            mylog.warning('Spectrum is totally empty!')
        else:
            self.flux_field = np.exp(-self.tau_field)

        if output_file is None:
            pass
        elif output_file.endswith('.h5') or output_file.endswith('.hdf5'):
            self._write_spectrum_hdf5(output_file)
        elif output_file.endswith('.fits'):
            self._write_spectrum_fits(output_file)
        else:
            self._write_spectrum_ascii(output_file)
        if output_absorbers_file is not None:
            self._write_absorbers_file(output_absorbers_file)

        del field_data
        return (self.lambda_field, self.flux_field)

    def error_func(self, flux):
        """
        Approximate the flux error for a spectrum.
        Many observational analysis programs require a flux error channel
        in addition to a flux channel.  So we create a zeroth order
        approximation of the flux error, simply by taking the square root
        of the flux.  Unfortunately, with flux normalized to be < 1, this
        would result in errors larger than the flux values themselves,
        so we normalize by an arbitrary signal-to-noise ratio, which by default
        is set to 100.  This yields a typical error for a normalized spectrum of
        sqrt(1.0*100)/100 = 0.1.  This assures our flux errors are smaller
        than our fluxes for most flux reasonable flux values.  Note that
        when a signal to noise ratio is specified for adding gaussian noise,
        it uses this updated value for estimating the errors.  SNR is set
        as an attribute of AbsorptionSpectrum directly (e.g., as.snr = N).

        **Parameters**

        :flux: array of floats

            The array of flux values
        """
        return np.sqrt(flux*self.snr)/self.snr

    def _apply_observing_redshift(self, field_data, use_peculiar_velocity,
                                 observing_redshift):
        """
        Change the redshifts of individual absorbers to account for the
        redshift at which the observer sits.

        The intermediate redshift that is seen by an observer
        at a redshift other than z=0 is z12, where z1 is the
        observing redshift and z2 is the emitted photon's redshift
        Hogg (2000) eq. 13:

        1 + z12 = (1 + z2) / (1 + z1)
        """
        if observing_redshift == 0.:
            # This is already assumed in the generation of the LightRay
            redshift = field_data['redshift']
            if use_peculiar_velocity:
                redshift_eff = field_data['redshift_eff']
        else:
            # The intermediate redshift that is seen by an observer
            # at a redshift other than z=0 is z12, where z1 is the
            # observing redshift and z2 is the emitted photon's redshift
            # Hogg (2000) eq. 13:
            # 1 + z12 = (1 + z2) / (1 + z1)
            redshift = ((1 + field_data['redshift']) / \
                        (1 + observing_redshift)) - 1.
            # Combining cosmological redshift and doppler redshift
            # into an effective redshift is found in Peacock's
            # Cosmological Physics eqn 3.75:
            # 1 + z_eff = (1 + z_cosmo) * (1 + z_doppler)
            if use_peculiar_velocity:
                redshift_eff = ((1 + redshift) * \
                                (1 + field_data['redshift_dopp'])) - 1.

        if not use_peculiar_velocity:
            redshift_eff = redshift

        return redshift, redshift_eff

    def _add_continua_to_spectrum(self, field_data, use_peculiar_velocity,
                                  observing_redshift=0., min_tau=1e-3):
        """
        Add continuum features to the spectrum.  Continuua are recorded as
        a name, associated field, wavelength, normalization value, and index.
        Continuua are applied at and below the denoted wavelength, where the
        optical depth decreases as a power law of desired index.  For positive
        index values, this means optical depth is highest at the denoted
        wavelength, and it drops with shorter and shorter wavelengths.
        Consequently, transmitted flux undergoes a discontinuous cutoff at the
        denoted wavelength, and then slowly increases with decreasing wavelength
        according to the power law.
        """

        if self._auto_lambda and self.lambda_field is None:
            mylog.warning(
                'Cannot add continuum with empty spectrum and lambda_min/max ' +
                'set to auto.')
            return

        # Change the redshifts of continuum sources to account for the
        # redshift at which the observer sits
        redshift, redshift_eff = self._apply_observing_redshift(field_data,
                                 use_peculiar_velocity, observing_redshift)

        # min_tau is the minimum optical depth value that warrants
        # accounting for an absorber.  for a single absorber, noticeable
        # continuum effects begin for tau = 1e-3 (leading to transmitted
        # flux of e^-tau ~ 0.999).  but we apply a cutoff to remove
        # absorbers with insufficient column_density to contribute
        # significantly to a continuum (see below).  because lots of
        # low column density absorbers can add up to a significant
        # continuum effect, we normalize min_tau by the n_absorbers.
        n_absorbers = field_data['dl'].size

        if n_absorbers == 0:
            mylog.info("No absorbers in path of LightRay.")
            return

        min_tau /= n_absorbers

        for continuum in self.continuum_list:

            # Normalization is in cm**-2, so column density must be as well
            column_density = (field_data[continuum['field_name']] *
                              field_data['dl']).in_units('cm**-2')
            if (column_density == 0).all():
                mylog.info("Not adding continuum %s: insufficient column density" % continuum['label'])
                continue

            # redshift_eff field combines cosmological and velocity redshifts
            if use_peculiar_velocity:
                delta_lambda = continuum['wavelength'] * redshift_eff
            else:
                delta_lambda = continuum['wavelength'] * redshift

            # right index of continuum affected area is wavelength itself
            this_wavelength = delta_lambda + continuum['wavelength']
            right_index = np.digitize(this_wavelength,
                                      self.lambda_field).clip(0, self.lambda_field.size)
            # left index of continuum affected area wavelength at which
            # optical depth reaches tau_min
            left_index = np.digitize((this_wavelength *
                              np.power((min_tau * continuum['normalization'] /
                                        column_density),
                                       (1. / continuum['index']))),
                              self.lambda_field).clip(0, self.lambda_field.size)

            # Only calculate the effects of continuua where normalized
            # column_density is greater than min_tau
            # because lower column will not have significant contribution
            valid_continuua = np.where(((column_density /
                                         continuum['normalization']) > min_tau) &
                                       (right_index - left_index > 1))[0]
            if valid_continuua.size == 0:
                mylog.info("Not adding continuum %s: insufficient column density or out of range" %
                    continuum['label'])
                continue

            pbar = get_pbar("Adding continuum - %s [%f A]: " % \
                                (continuum['label'], continuum['wavelength']),
                            valid_continuua.size)

            # Tau value is (wavelength / continuum_wavelength)**index /
            #              (column_dens / norm)
            # i.e. a power law decreasing as wavelength decreases

            # Step through the absorber list and add continuum tau for each to
            # the total optical depth for all wavelengths
            for i, lixel in enumerate(valid_continuua):
                cont_tau = \
                    np.power((self.lambda_field[left_index[lixel] :
                                                right_index[lixel]] /
                                   this_wavelength[lixel]), \
                              continuum['index']) * \
                    (column_density[lixel] / continuum['normalization'])
                self.tau_field[left_index[lixel]:right_index[lixel]] += cont_tau.d
                pbar.update(i)
            pbar.finish()

    def _add_lines_to_spectrum(self, field_data, use_peculiar_velocity,
                               output_absorbers_file, store_observables,
                               subgrid_resolution=10, observing_redshift=0.,
                               njobs=-1, min_tau=1e-3):
        """
        Add the absorption lines to the spectrum.
        """

        if len(self.line_list) == 0:
            return

        if self.bin_space == 'velocity':
            wavelength_zero_point = self.line_list[0]['wavelength']

        # Change the redshifts of individual absorbers to account for the
        # redshift at which the observer sits
        redshift, redshift_eff = self._apply_observing_redshift(field_data,
                                 use_peculiar_velocity, observing_redshift)

        # step through each ionic transition (e.g. HI, HII, MgII) specified
        # and deposit the lines into the spectrum
        for store, line in parallel_objects(self.line_list, njobs=njobs,
                                            storage=self.line_observables_dict):
            column_density = field_data[line['field_name']] * field_data['dl']
            if (column_density < 0).any():
                mylog.warning(
                    "Setting negative densities for field %s to 0! Bad!" % line['field_name'])
                np.clip(column_density, 0, np.inf, out=column_density)
            if (column_density == 0).all():
                mylog.info("Not adding line %s: insufficient column density" % line['label'])
                continue

            # redshift_eff field combines cosmological and velocity redshifts
            # so delta_lambda gives the offset in angstroms from the rest frame
            # wavelength to the observed wavelength of the transition
            if use_peculiar_velocity:
                delta_lambda = line['wavelength'] * redshift_eff
            else:
                delta_lambda = line['wavelength'] * redshift
            # lambda_obs is central wavelength of line after redshift
            lambda_obs = (line['wavelength'] + delta_lambda).to('angstrom')

            # either the observed wavelength or velocity offset
            if self.bin_space == 'wavelength':
                my_obs = lambda_obs[:]
            elif self.bin_space == 'velocity':
                my_obs = c_kms * \
                  (lambda_obs - wavelength_zero_point) / \
                  wavelength_zero_point
                my_obs.convert_to_units(_bin_space_units[self.bin_space])
            else:
                raise RuntimeError('What bin_space is this?')

            # the total number of absorbers per transition
            n_absorbers = len(lambda_obs)

            # thermal broadening b parameter
            thermal_b =  np.sqrt((2 * boltzmann_constant_cgs *
                                  field_data['temperature']) /
                                  line['atomic_mass'])

            # the actual thermal width of the lines
            thermal_width = (lambda_obs * thermal_b /
                             c_kms).to('angstrom')

            # Sanitize units for faster runtime of the tau_profile machinery.
            lambda_0 = line['wavelength'].d  # line's rest frame; angstroms
            cdens = column_density.in_units("cm**-2").d # cm**-2
            thermb = thermal_b.to('cm/s').d  # thermal b coefficient; cm / s
            dlambda = delta_lambda.d  # lambda offset; angstroms
            # Array to store sum of the tau values for each index in the
            # light ray that is deposited to the final spectrum
            if store_observables:
                tau_ray = np.zeros(cdens.size)
            if use_peculiar_velocity:
                vlos = field_data['velocity_los'].in_units("km/s").d # km/s
            else:
                vlos = np.zeros(field_data['temperature'].size)

            # When we actually deposit the voigt profile, sometimes we will
            # have underresolved lines (ie lines with smaller widths than
            # the spectral bin size).  Here, we create virtual wavelength bins
            # small enough in width to well resolve each line, deposit the
            # voigt profile into them, then numerically integrate their tau
            # values and sum them to redeposit them into the actual spectral
            # bins.

            # virtual bins (vbins) will be:
            # 1) <= the bin_width; assures at least as good as spectral bins
            # 2) <= 1/10th the thermal width; assures resolving voigt profiles
            #   (actually 1/subgrid_resolution value, default is 1/10)
            # 3) a bin width will be divisible by vbin_width times a power of
            #    10; this will assure we don't get spikes in the deposited
            #    spectra from uneven numbers of vbins per bin

            if self.bin_space == 'wavelength':
                my_width = thermal_width
            elif self.bin_space == 'velocity':
                my_width = thermal_b
            else:
                raise RuntimeError('What bin space is this?')

            resolution = my_width / self.bin_width
            n_vbins_per_bin = (10 ** (np.ceil( np.log10( subgrid_resolution /
                               resolution) ).clip(0, np.inf) ) ).astype('int')
            vbin_width = self.bin_width.d / n_vbins_per_bin

            # a note to the user about which lines components are unresolved
            if (my_width < self.bin_width).any():
                mylog.info("%d out of %d line components will be " +
                            "deposited as unresolved lines.",
                            (thermal_width < self.bin_width).sum(),
                            n_absorbers)

            # Keep track of the lambda field before depositing a new line
            # so we can add the current_tau_field and the tau_field together.
            last_lambda_field = self.lambda_field

            # provide a progress bar with information about lines processsed
            pbar = get_pbar("Adding line - %s [%f A]: " % \
                            (line['label'], line['wavelength']), n_absorbers)

            # for a given transition, step through each location in the
            # observed spectrum where it occurs and deposit a voigt profile
            for i in parallel_objects(np.arange(n_absorbers), njobs=-1):

                # if there is a ray element with temperature = 0 or column
                # density = 0, skip it
                if (thermal_b[i] == 0.) or (cdens[i] == 0.):
                    pbar.update(i)
                    continue

                # the virtual window into which the line is deposited initially
                # spans a region of 2 coarse spectral bins
                # (one on each side of the center_index) but the window
                # can expand as necessary.
                # it will continue to expand until the tau value in the far
                # edge of the wings is less than the min_tau value or it
                # reaches the edge of the spectrum
                window_width_in_bins = 2

                # Widen wavelength window until optical depth falls below min_tau
                # value at the ends to assure that the wings of a line have been
                # fully resolved.
                while True:

                    # calculate wavelength window
                    if self._auto_lambda and self.lambda_field is None:
                        my_lambda_min = my_obs[i] - \
                          window_width_in_bins * self.bin_width / 2
                        # round off to multiple of bin_width
                        my_lambda_min = self.bin_width * \
                          np.ceil(my_lambda_min / self.bin_width)
                        my_lambda = my_lambda_min + \
                          self.bin_width * np.arange(window_width_in_bins)

                    else:
                        my_lambda = self.lambda_field

                    # we want to know the bin index in the lambda_field array
                    # where each line has its central wavelength after being
                    # redshifted.  however, because we don't know a priori how wide
                    # a line will be (ie DLAs), we have to include bin indices
                    # *outside* the spectral range of the AbsorptionSpectrum
                    # object.  Thus, we find the "equivalent" bin index, which
                    # may be <0 or >the size of the array.  In the end, we deposit
                    # the bins that actually overlap with the AbsorptionSpectrum's
                    # range in lambda.

                    left_index, center_index, right_index = \
                      self._get_bin_indices(
                          my_lambda, self.bin_width,
                          my_obs[i], window_width_in_bins)
                    n_vbins = window_width_in_bins * n_vbins_per_bin[i]

                    # the array of virtual bins in lambda space
                    vbins = \
                        np.linspace(my_lambda.d[0] + self.bin_width.d * left_index,
                                    my_lambda.d[0] + self.bin_width.d * right_index,
                                    n_vbins, endpoint=False)

                    if self.bin_space == 'wavelength':
                        my_vbins = vbins
                    elif self.bin_space == 'velocity':
                        my_vbins = vbins * \
                          wavelength_zero_point.d / c_kms.d + \
                          wavelength_zero_point.d
                    else:
                        raise RuntimeError('What bin_space is this?')

                    # the virtual bins and their corresponding opacities
                    my_vbins, vtau = \
                        tau_profile(
                            lambda_0, line['f_value'], line['gamma'],
                            thermb[i], cdens[i],
                            delta_lambda=dlambda[i], lambda_bins=my_vbins)

                    # If tau has not dropped below min tau threshold by the
                    # edges (ie the wings), then widen the wavelength
                    # window and repeat process.
                    if (vtau[0] < min_tau and vtau[-1] < min_tau):
                        if self._auto_lambda:
                            self._create_auto_field_arrays(
                                left_index, right_index, my_lambda)
                            left_index, center_index, right_index = \
                              self._get_bin_indices(
                                  self.lambda_field, self.bin_width,
                                  my_obs[i], window_width_in_bins)

                        break
                    window_width_in_bins *= 2

                if center_index is None:
                    pbar.update(i)
                    continue

                # Numerically integrate the virtual bins to calculate a
                # virtual "equivalent width" of optical depth; then sum these
                # virtual equivalent widths in tau and deposit back into each
                # original spectral tau bin
                # Please note: this is not a true equivalent width in the
                # normal use of the word by observers.  It is an equivalent
                # with in tau, not in flux, and is only used internally in
                # this subgrid deposition as EW_tau.
                vEW_tau = vtau * vbin_width[i]
                EW_tau = np.zeros(right_index - left_index)
                EW_tau_indices = np.arange(left_index, right_index)
                for k, val in enumerate(EW_tau_indices):
                    EW_tau[k] = vEW_tau[n_vbins_per_bin[i] * k:
                                        n_vbins_per_bin[i] * (k + 1)].sum()
                EW_tau = EW_tau/self.bin_width.d

                # only deposit EW_tau bins that actually intersect the original
                # spectral wavelength range (i.e. lambda_field)

                # if EW_tau bins don't intersect the original spectral range at
                # all then skip the deposition
                if ((left_index >= self.lambda_field.size) or \
                    (right_index < 0)):
                    pbar.update(i)
                    continue

                # otherwise, determine how much of the original spectrum
                # is intersected by the expanded line window to be deposited,
                # and deposit the Equivalent Width in tau into that intersecting
                # window in the original spectrum's tau array
                else:
                    intersect_left_index = max(left_index, 0)
                    intersect_right_index = min(right_index, self.lambda_field.size)
                    EW_tau_deposit = EW_tau[(intersect_left_index - left_index): \
                                            (intersect_right_index - left_index)]
                    self.current_tau_field[intersect_left_index:intersect_right_index] \
                        += EW_tau_deposit
                    if store_observables:
                        tau_ray[i] = np.sum(EW_tau_deposit)
                # write out absorbers to file if the column density of
                # an absorber is greater than the specified "label_threshold"
                # of that absorption line
                if output_absorbers_file and \
                   line['label_threshold'] is not None and \
                   cdens[i] >= line['label_threshold']:

                    if use_peculiar_velocity:
                        peculiar_velocity = vlos[i]
                    else:
                        peculiar_velocity = 0.0
                    self.absorbers_list.append({'label': line['label'],
                                                'wavelength': (lambda_0 + dlambda[i]),
                                                'column_density': column_density[i],
                                                'b_thermal': thermal_b[i],
                                                'redshift': redshift[i],
                                                'redshift_eff': redshift_eff[i],
                                                'v_pec': peculiar_velocity})
                pbar.update(i)
            pbar.finish()

            # Expand the tau_field array to match the updated wavelength
            # array from the last line deposition.
            self._adjust_field_array(last_lambda_field, self.lambda_field,
                                     "tau_field")

            if self.current_tau_field is not None:
                # Now add the current_tau_field.
                self.tau_field += self.current_tau_field

            ## Check keyword before storing any observables
            if store_observables:
                # If running in parallel, make sure that the observable
                # quantities for the dictionary are combined correctly.
                comm = _get_comm(())

                if self._auto_lambda:
                    global_lambda_field = self._get_global_lambda_field(comm=comm)
                    self._adjust_field_array(self.lambda_field, global_lambda_field,
                                             "current_tau_field")

                if comm.size > 1:
                    obs_dict_fields = \
                      [column_density, tau_ray, self.current_tau_field,
                       delta_lambda, lambda_obs, thermal_b, thermal_width]
                    obs_dict_fields = [comm.mpi_allreduce(field,op="sum")
                                       for field in obs_dict_fields]

                # Calculate the flux decrement equivalent width (the true
                # equivalent width!) for use in post-processing
                if self.current_tau_field is None:
                    EW = 0.
                else:
                    EW = np.sum(1-np.exp(-self.current_tau_field))*self.bin_width
                # Update the line_observables_dict with values for this line
                obs_dict = {"column_density":column_density,
                            "tau_ray":tau_ray,
                            "EW":EW,
                            "delta_lambda":delta_lambda,
                            "lambda_obs":lambda_obs,
                            "thermal_b":thermal_b,
                            "thermal_width":thermal_width}
                if self.bin_space == 'velocity':
                    obs_dict['velocity_offset'] = my_obs
                store.result_id = line['label']
                store.result = obs_dict
                ## Can only delete these if in this statement:
                del obs_dict, tau_ray

            self.current_tau_field = None

            # These always need to be deleted
            del column_density, delta_lambda, lambda_obs, my_obs, \
                thermal_b, thermal_width, cdens, thermb, dlambda, \
                vlos, resolution, vbin_width, n_vbins, n_vbins_per_bin


        comm = _get_comm(())
        if self._auto_lambda:
            new_lambda = self._get_global_lambda_field(comm=comm)
            self._adjust_field_array(self.lambda_field, new_lambda,
                                     "tau_field")
            self.lambda_field = new_lambda
        self.tau_field = comm.mpi_allreduce(self.tau_field, op="sum")
        if output_absorbers_file:
            self.absorbers_list = comm.par_combine_object(
                self.absorbers_list, "cat", datatype="list")


    def _get_bin_indices(self, lambda_field, dlambda, lambda_obs,
                         window_width_in_bins):
        """
        Return the indices of the lambda field corresponding to
        the lower limit, line center, and upper limit.
        """

        if lambda_field is None or lambda_field.size == 0:
            return None, None, None

        # this equation gives us the "equivalent" bin index for each line
        # if it were placed into the self.lambda_field array
        center_index = ((lambda_obs - lambda_field[0]) /
                        dlambda).d
        center_index = int(np.ceil(center_index))
        left_index = (center_index - window_width_in_bins//2)
        right_index = (center_index + window_width_in_bins//2)

        return left_index, center_index, right_index

    def _get_global_lambda_field(self, comm=None):
        """
        Get the maximum lambda field bounds and create a new array.
        """
        if comm is None:
            comm = _get_comm(())

        if self.lambda_field is None:
            my_min = np.inf
            my_max = -np.inf
        else:
            my_min = self.lambda_field[0]
            my_max = self.lambda_field[-1]

        lf_min = comm.mpi_allreduce(my_min, op="min")
        lf_max = comm.mpi_allreduce(my_max, op="max")

        if lf_min != np.inf:
            lf_min = np.round(lf_min / self.bin_width) * self.bin_width
            lf_max = np.round(lf_max / self.bin_width) * self.bin_width
            n_lambda = self._get_field_size(lf_min, lf_max, self.bin_width)
            new_lambda = self._create_lambda_field(lf_min, lf_max, n_lambda)
        else:
            new_lambda = None

        return new_lambda

    def _create_auto_field_arrays(self, left_index, right_index,
                                  my_lambda):
        """
        Adjust the current wavelength window to encompass the full
        window of the next feature we are depositing.

        my_lambda is the current wavelength window and the left and
        right index are the bounds of the desired wavelength window.
        If the bounds are outside the current window, we enlarge it
        to fit and adjust the tau field accordingly.
        """

        if left_index >= 0 and right_index < my_lambda.size:
            return

        dlambda = my_lambda[1] - my_lambda[0]
        new_lambda_min = my_lambda[0] + \
          dlambda * min(0, left_index)
        new_lambda_max = my_lambda[0] + \
          dlambda * max(my_lambda.size-1, right_index)

        if str(self.lambda_min) != 'auto':
            new_lambda_min = max(new_lambda_min, self.lambda_min)
        if str(self.lambda_max) != 'auto':
            new_lambda_max = min(new_lambda_max, self.lambda_max)
        if new_lambda_min >= new_lambda_max:
            return

        n_lambda = self._get_field_size(
            new_lambda_min, new_lambda_max, dlambda)
        new_lambda = \
          self._create_lambda_field(new_lambda_min, new_lambda_max,
                                    n_lambda)

        if self._current_tau_field is not None:
            self._adjust_field_array(self.lambda_field, new_lambda,
                                     "current_tau_field")

        self.lambda_field = new_lambda

    def _adjust_field_array(self, old_lambda, new_lambda, array_name):
        """
        Adjust the field array associated with the old wavelength array
        so that it lines up correctly with the new wavelength array.
        """
        if old_lambda is None or new_lambda is None:
            return

        start_index = np.digitize(old_lambda[0], new_lambda) - 1
        new_array = np.zeros(new_lambda.size)

        old_array = getattr(self, array_name)
        new_array[start_index:start_index+old_array.size] = old_array
        setattr(self, array_name, new_array)

    @parallel_root_only
    def _write_absorbers_file(self, filename):
        """
        Write out ASCII list of all substantial absorbers found in spectrum
        """
        if filename is None:
            return
        if self.tau_field is None:
            return
        mylog.info("Writing absorber list: %s.", filename)
        self.absorbers_list.sort(key=lambda obj: obj['wavelength'])
        f = open(filename, 'w')
        f.write('#%-14s %-14s %-12s %-14s %-15s %-9s %-10s\n' %
                ('Wavelength', 'Line', 'N [cm^-2]', 'b [km/s]', 'z_cosmo', \
                 'z_eff', 'v_pec [km/s]'))
        for line in self.absorbers_list:
            f.write('%-14.6f %-14ls %e %e % e % e % e\n' % (line['wavelength'], \
                line['label'], line['column_density'], line['b_thermal'], \
                line['redshift'], line['redshift_eff'], line['v_pec']))
        f.close()

    @parallel_root_only
    def _write_spectrum_ascii(self, filename):
        """
        Write spectrum to an ascii file.
        """
        if self.tau_field is None:
            return
        mylog.info("Writing spectrum to ascii file: %s.", filename)
        f = open(filename, 'w')
        f.write("# wavelength[A] tau flux flux_error\n")
        for i in range(self.lambda_field.size):
            f.write("%e %e %e %e\n" % (self.lambda_field[i],
                                    self.tau_field[i],
                                    self.flux_field[i],
                                    self.error_func(self.flux_field[i])))
        f.close()

    @parallel_root_only
    def _write_spectrum_fits(self, filename):
        """
        Write spectrum to a fits file.
        """
        if self.tau_field is None:
            return
        mylog.info("Writing spectrum to fits file: %s.", filename)
        col1 = pyfits.Column(name='wavelength', format='E', array=self.lambda_field)
        col2 = pyfits.Column(name='tau', format='E', array=self.tau_field)
        col3 = pyfits.Column(name='flux', format='E', array=self.flux_field)
        col4 = pyfits.Column(name='flux_error', format='E', array=self.error_func(self.flux_field))
        cols = pyfits.ColDefs([col1, col2, col3, col4])
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(filename, overwrite=True)

    @parallel_root_only
    def _write_spectrum_hdf5(self, filename):
        """
        Write spectrum to an hdf5 file.

        """
        if self.tau_field is None:
            return
        mylog.info("Writing spectrum to hdf5 file: %s.", filename)
        output = h5py.File(filename, 'w')
        output.create_dataset('wavelength', data=self.lambda_field)
        output.create_dataset('tau', data=self.tau_field)
        output.create_dataset('flux', data=self.flux_field)
        output.create_dataset('flux_error', data=self.error_func(self.flux_field))
        output.close()
