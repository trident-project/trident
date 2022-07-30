"""
FireflyGenerator class and member functions.

"""

import numpy as np
import os

from yt.loaders import \
    load
from yt.data_objects.data_containers import \
    YTDataContainer
from yt.data_objects.static_output import \
    Dataset
from yt.funcs import \
    mylog, \
    YTArray
from yt.utilities.exceptions import \
    YTFieldNotFound
from yt.utilities.on_demand_imports import \
    _firefly

from trident.config import \
    ion_table_dir, \
    ion_table_file, \
    ion_table_filepath
from trident.ion_balance import \
    add_ion_number_density_field, \
    atomic_mass
from trident.line_database import \
    LineDatabase

class FireflyGenerator( object ):
    """
    ## OLD TEXT COMMENTED OUT ##
    ## Keeping it so that I can easily use consistent formatting. ##
    # Preferred class for generating, storing, and plotting absorption-line spectra.
    # SpectrumGenerator is a subclass of yt's AbsorptionSpectrum class
    # with additional functionality like line lists, post-processing to
    # include Milky Way foreground, quasar backgrounds, applying line-spread
    # functions, and plotting.

    # User first specifies the telescope/instrument used for generating spectra
    # (e.g. 'COS' for the Cosmic Origins Spectrograph aboard the
    # Hubble Space Telescope).  This can be done by naming the
    # :class:`~trident.Instrument`, or manually setting all of the spectral
    # characteristics including ``lambda_min``, ``lambda_max``, ``lsf_kernel``,
    # and ``n_lambda`` or ``dlambda``.  If none of these arguments are set,
    # defaults to 'COS' as the default instrument covering 1150-1450 Angstroms
    # with a binsize (``dlambda``) of 0.1 Angstroms.

    # Once a :class:`~trident.SpectrumGenerator` has been initialized, pass it
    # ``LightRay`` objects using :class:`~trident.SpectrumGenerator.make_spectrum`
    # to actually generate the spectra themselves.  Then one can post-process,
    # plot, or save them using
    # :class:`~trident.SpectrumGenerator.add_milky_way_foreground`,
    # :class:`~trident.SpectrumGenerator.add_qso_spectrum`,
    # :class:`~trident.SpectrumGenerator.apply_lsf`,
    # :class:`~trident.SpectrumGenerator.save_spectrum`, and
    # :class:`~trident.SpectrumGenerator.plot_spectrum`.

    # **Parameters**

    # :instrument: string, optional

    #     The spectrograph to use.  Currently, the only named options are
    #     different observing modes of the Cosmic Origins Spectrograph 'COS':
    #     'COS-G130M', 'COS-G160M', and 'COS-G140L' as well as 'COS' which
    #     aliases to 'COS-G130M'. These automatically set the ``lambda_min``,
    #     ``lambda_max``, ``dlambda`` and ``lsf_kernel``s appropriately.
    #     If you're going to set ``lambda_min``, ``lambda_max``, et al manually,
    #     leave this set to None.
    #     Default: None

    # :lambda_min: float, YTQuantity, or 'auto'

    #    lower wavelength bound in angstroms or velocity bound in km/s
    #    (if bin_space set to 'velocity'). If set to 'auto', the lower
    #    bound will be automatically adjusted to encompass all absorption
    #    lines. The window will not be expanded for continuum features,
    #    only absorption lines.

    # :lambda_max: float, YTQuantity, or 'auto'

    #    upper wavelength bound in angstroms or velocity bound in km/s
    #    (if bin_space set to 'velocity'). If set to 'auto', the upper
    #    bound will be automatically adjusted to encompass all absorption
    #    lines. The window will not be expanded for continuum features,
    #    only absorption lines.

    # :n_lambda: int

    #     The number of bins in the spectrum (inclusive), so if
    #     extrema = 10 and 20, and dlambda (binsize) = 1, then n_lambda = 11.
    #     Default: None

    # :dlambda: float

    #     size of the wavelength bins in angstroms or velocity bins in km/s.
    #     Default: None

    # :bin_space: 'wavelength' or 'velocity'

    #     Sets the dimension in which spectra are created. If set to
    #     wavelength, the resulting spectra are flux (or tau) vs. observed
    #     wavelength. If set to velocity, the spectra are flux vs.
    #     velocity offset from the rest wavelength of the absorption line.
    #     Default: wavelength

    # :lsf_kernel: string, optional

    #     The filename for the LSF kernel. Files are found in
    #     trident.__path__/data/lsf_kernels or current working directory.
    #     Only necessary if you are applying an LSF to the spectrum in
    #     postprocessing.
    #     Default: None

    # :line_database: string or :class:`~trident.LineDatabase`, optional

    #     A text file listing the various lines to insert into the line database,
    #     or a :class:`~trident.LineDatabase` object in memory. The line database
    #     provides a list of all possible lines that could be added to the
    #     spectrum. For a text file, it should have 4 tab-delimited columns of
    #     name (e.g. MgII), wavelength in angstroms, gamma of transition, and
    #     f-value of transition. See example datasets in
    #     trident.path/data/line_lists for examples.
    #     Default: lines.txt

    # :ionization_table: hdf5 file, optional

    #     An HDF5 file used for computing the ionization fraction of the gas
    #     based on its density, temperature, metallicity, and redshift.  If
    #     set to None, will use the ion table defined in your Trident config
    #     file.
    #     Default: None

    # **Example**

    # Create a one-zone ray, and generate a COS spectrum from that ray.

    # >>> import trident
    # >>> ray = trident.make_onezone_ray()
    # >>> sg = trident.SpectrumGenerator('COS')
    # >>> sg.make_spectrum(ray)
    # >>> sg.plot_spectrum('spec_raw.png')

    # Create a one-zone ray at redshift 0.5, and generate a spectrum with 1
    # angstrom spectral bins from 2000-4000 angstroms, then post-process by
    # adding a MW foreground a QSO background at z=0.5 and add a boxcar line
    # spread function of 100 angstroms width.  Plot it and save the figure to
    # 'spec_final.png'.

    # >>> import trident
    # >>> ray = trident.make_onezone_ray(redshift=0.5)
    # >>> sg = trident.SpectrumGenerator(lambda_min=2000, lambda_max=4000,
    # ... dlambda=1)
    # >>> sg.make_spectrum(ray)
    # >>> sg.add_qso_spectrum(emitting_redshift=.5)
    # >>> sg.add_milky_way_foreground()
    # >>> sg.apply_lsf(function='boxcar', width=100)
    # >>> sg.plot_spectrum('spec_final.png')
    """
    def __init__(self, line_database='lines.txt',
                ionization_table=None ):

        if isinstance(line_database, LineDatabase):
            self.line_database = line_database
        else:
            # instantiate the LineDatabase
            self.line_database = LineDatabase(line_database)

        self.observing_redshift = 0.

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

    def create_particle_group(self, ray,
                    fields_to_include=None,
                    fields_units=None,
                    coordinate_units="kpc",
                    velocity_units="km/s",
                    center = None,
                    observing_redshift=0.0,
                    lines='all',
                    ):

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

        # temporary fix for yt-4.0 ytdata selection issue
        # ray.domain_left_edge = ray.domain_left_edge.to('code_length')
        # ray.domain_right_edge = ray.domain_right_edge.to('code_length')

        active_lines = self.line_database.parse_subset(lines)

        # Make sure we've produced all the necessary
        # derived fields if they aren't native to the data
        for line in active_lines:
            # try to determine if field for line is present in dataset
            # if successful, means line.field is in ds.derived_field_list
            # otherwise we probably need to add the field to the dataset
            try:
                ad._determine_fields(line.field)[0]
            except YTFieldNotFound:
                fname = line.field[1] # grab field string from field name tuple
                my_ion = \
                  fname[:fname.find("number_density")]
                on_ion = my_ion.split("_")
                # Add the field using ray's density, temperature, & metallicity
                my_lev = int(on_ion[1][1:]) + 1
                mylog.info("Creating %s from ray's fields." % (line.field[1]))
                add_ion_number_density_field(on_ion[0], my_lev, ray,
                                 ionization_table=self.ionization_table)


        # Get coordinates and velocities
        coordinates = np.array([
            ad['grid', x_i].in_units( coordinate_units ) for x_i in [ 'x', 'y', 'z' ]
        ]).transpose() * ray.quan( 1, coordinate_units )
        mylog.info( 'Performing a temporary fix until we figure out why the coordinates are off-center!' )
        coordinates -= coordinates[-1]
        if center is not None:
            coordinates += center
        velocities = np.array([
            ad['grid', 'relative_velocity_' + x_i ].in_units( velocity_units ) for x_i in [ 'x', 'y', 'z' ]
        ]).transpose() * ray.quan( 1, velocity_units )

        ## handle default arguments
        if fields_to_include is None:
            fields_to_include = []
        if fields_units is None:
            fields_units = []

        ## handle input validation, if any
        if len(fields_units) != len(fields_to_include):
            raise RuntimeError("Each requested field must have units.")

        ## explicitly go after the fields we want
        field_arrays = []
        field_names = []
        unavailable_fields = []
        for field, units in zip(fields_to_include, fields_units):
            ## determine if you want to take the log of the field for Firefly
            log_flag = "log(" in units

            ## read the field array from the dataset
            try:
                this_field_array = ad['grid', field]
            except YTFieldNotFound:
                unavailable_fields.append( field )
                continue

            ## fix the units string and prepend 'log' to the field for
            ##  the UI name
            if log_flag:
                units = units[len("log(") : -1]
                field = f"log{field}"

            ## perform the unit conversion and take the log if
            ##  necessary.
            this_field_array.in_units(units)
            if log_flag:
                this_field_array = np.log10(this_field_array)

            ## add this array to the tracked arrays
            field_arrays.append( this_field_array )
            field_names.append( field )

        # Print fields skipped because they were unavailable
        if len( unavailable_fields ) > 0:
            print(
                'For ray, unable to retrieve these fields: {}'.format(
                    unavailable_fields
                )
            )

        ## create a firefly ParticleGroup for this particle type
        # Include fields if available
        if len( field_arrays ) != 0:
            ParticleGroup_kwargs = {
                'field_arrays': field_arrays,
                'field_names': field_names,
            }
        else:
            ParticleGroup_kwargs = {}
        pg = _firefly.data_reader.ParticleGroup(
            UIname='ray',
            coordinates=coordinates,
            velocities=velocities,
            **ParticleGroup_kwargs
        )

        return pg