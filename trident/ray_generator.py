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

import six

from trident.light_ray import \
    LightRay
from yt.convenience import \
    load, \
    simulation
from trident.config import \
    ion_table_filepath
from trident.line_database import \
    LineDatabase, \
    uniquify
from trident.roman import \
    from_roman
from yt.data_objects.static_output import \
    Dataset
from trident.ion_balance import \
    add_ion_number_density_field, \
    atomic_number
from trident.utilities import \
    _check_sampling_types_match

def make_simple_ray(dataset_file, start_position, end_position,
                    lines=None, ftype="gas", fields=None,
                    solution_filename=None, data_filename=None,
                    trajectory=None, redshift=None, field_parameters=None,
                    setup_function=None, load_kwargs=None,
                    line_database=None, ionization_table=None):
    """
    Create a yt LightRay object for a single dataset (eg CGM).  This is a
    wrapper function around yt's LightRay interface to reduce some of the
    complexity there.

    A simple ray is a straight line passing through a single dataset
    where each gas cell intersected by the line is sampled for the desired
    fields and stored.  Several additional fields are created and stored
    including ``dl`` to represent the path length in space
    for each element in the ray, ``v_los`` to represent the line of
    sight velocity along the ray, and ``redshift``, ``redshift_dopp``, and
    ``redshift_eff`` to represent the cosmological redshift, doppler redshift
    and effective redshift (combined doppler and cosmological) for each
    element of the ray.

    A simple ray is typically specified by its start and end positions in the
    dataset volume.  Because a simple ray only probes a single output, it
    lacks foreground absorbers between the observer at z=0 and the redshift
    of the dataset that one would naturally encounter.  Thus it is usually
    only appropriate for studying the circumgalactic medium rather than
    the intergalactic medium.

    This function can accept a yt dataset already loaded in memory,
    or it can load a dataset if you pass it the dataset's filename and
    optionally any load_kwargs or setup_function necessary to load/process it
    properly before generating the LightRay object.

    The :lines: keyword can be set to automatically add all fields to the
    resulting ray necessary for later use with the SpectrumGenerator class.
    If the necessary fields do not exist for your line of choice, they will
    be added to your dataset before adding them to the ray.

    If using the :lines: keyword with an SPH dataset, it is very important
    to set the :ftype: keyword appropriately, or you may end up calculating
    ion fields by interpolating on data already smoothed to the grid.  This is
    generally not desired.

    **Parameters**

    :dataset_file: string or yt Dataset object

        Either a yt dataset or the filename of a dataset on disk.  If you are
        passing it a filename, consider usage of the ``load_kwargs`` and
        ``setup_function`` kwargs.

    :start_position, end_position: list of floats or YTArray object

        The coordinates of the starting and ending position of the desired
        ray.  If providing a raw list, coordinates are assumed to be in
        code length units, but if providing a YTArray, any units can be
        specified.

                    lines=None, fields=None, solution_filename=None,
                    data_filename=None, trajectory=None, redshift=None,
                    line_database=None, ftype="gas",
                    setup_function=None, load_kwargs=None,
                    ionization_table=None):

    :lines: list of strings, optional

        List of strings that determine which fields will be added to the ray
        to support line deposition to an absorption line spectrum.  List can
        include things like "C", "O VI", or "Mg II ####", where #### would be
        the integer wavelength value of the desired line.  If set to 'all',
        includes all possible ions from H to Zn. :lines: can be used
        in conjunction with :fields: as they will not override each other.
        If using the :lines: keyword with an SPH dataset, it is very important
        to set the :ftype: keyword appropriately, or you may end up calculating
        ion fields by interpolating on data already smoothed to the grid.
        This is generally not desired.
        Default: None

    :ftype: string, optional

        For use with the :lines: keyword.  It is the field type of the fields to
        be added.  It is the first string in the field tuple e.g. "gas" in
        ("gas", "O_p5_number_density"). For SPH datasets, it is important to
        set this to the field type of the gas particles in your dataset
        (e.g. 'PartType0'), as it determines the source data for the ion
        fields to be added. If you leave it set to "gas", it will calculate
        the ion fields based on the hydro fields already smoothed on the grid,
        which is usually not desired.
        Default: "gas"

    :fields: list of strings, optional

        The list of which fields to store in the output LightRay.
        See :lines: keyword for additional functionality that will add fields
        necessary for creating absorption line spectra for certain line
        features.
        Default: None

    :solution_filename: string, optional

        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    :data_filename: string, optional

        Output filename for ray data stored as an HDF5 file.  Note that
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    :trajectory: list of floats, optional

        The (r, theta, phi) direction of the LightRay.  Use either end_position
        or trajectory, but not both.
        Default: None

    :redshift: float, optional

        Sets the highest cosmological redshift of the ray.  By default, it will
        use the cosmological redshift of the dataset, if set, and if not set,
        it will use a redshift of 0.
        Default: None

    :field_parameters: optional, dict
        Used to set field parameters in light rays. For example,
        if the 'bulk_velocity' field parameter is set, the relative
        velocities used to calculate peculiar velocity will be adjusted
        accordingly.
        Default: None.

    :setup_function: function, optional

        A function that will be called on the dataset as it is loaded but
        before the LightRay is generated.  Very useful for adding derived
        fields and other manipulations of the dataset prior to LightRay
        creation.
        Default: None

    :load_kwargs: dict, optional

        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    :line_database: string, optional

        For use with the :lines: keyword. If you want to limit the available
        ion fields to be added to those available in a particular subset,
        you can use a :class:`~trident.LineDatabase`.  This means when you
        set :lines:='all', it will only use those ions present in the
        corresponding LineDatabase.  If :LineDatabase: is set to None,
        and :lines:='all', it will add every ion of every element up to Zinc.
        Default: None

    :ionization_table: string, optional

        For use with the :lines: keyword.  Path to an appropriately formatted
        HDF5 table that can be used to compute the ion fraction as a function
        of density, temperature, metallicity, and redshift.  When set to None,
        it uses the table specified in ~/.trident/config
        Default: None

    **Example**

    Generate a simple ray passing from the lower left corner to the upper
    right corner through some Gizmo dataset where gas particles are
    ftype='PartType0':

    >>> import trident
    >>> import yt
    >>> ds = yt.load('path/to/dataset')
    >>> ray = trident.make_simple_ray(ds,
    ... start_position=ds.domain_left_edge, end_position=ds.domain_right_edge,
    ... lines=['H', 'O', 'Mg II'], ftype='PartType0')
    """
    if load_kwargs is None:
        load_kwargs = {}
    if fields is None:
        fields = []
    if data_filename is None:
        data_filename = 'ray.h5'

    if isinstance(dataset_file, str):
        ds = load(dataset_file, **load_kwargs)
    elif isinstance(dataset_file, Dataset):
        ds = dataset_file

    lr = LightRay(ds, load_kwargs=load_kwargs)

    if ionization_table is None:
        ionization_table = ion_table_filepath

    # Include some default fields in the ray to assure it's processed correctly.

    fields = _add_default_fields(ds, fields)

    # If 'lines' kwarg is set, we need to get all the fields required to
    # create the desired absorption lines in the grid format, since grid-based
    # fields are what are directly probed by the LightRay object.

    # We first determine what fields are necessary for the desired lines, and
    # inspect the dataset to see if they already exist.  If so, we add them
    # to the field list for the ray.  If not, we have to create them.

    if lines is not None:

        ion_list = _determine_ions_from_lines(line_database, lines)

        sampling_type = _check_sampling_types_match(ds, ftype)

        fields, fields_to_add_to_ds = _determine_fields_from_ions(ds, ion_list,
                                        fields, ftype, sampling_type)

        # actually add the fields we need to add to the dataset
        for atom, ion_state in fields_to_add_to_ds:
            add_ion_number_density_field(atom, ion_state, ds,
                                         ftype=ftype,
                                         ionization_table=ionization_table,
                                         sampling_type=sampling_type)

    # To assure there are no fields that are double specified or that collide
    # based on being specified as "density" as well as ("gas", "density"),
    # we will just assume that all non-tuple fields requested are ftype "gas".
    for i in range(len(fields)):
        if isinstance(fields[i], str):
            fields[i] = ('gas', fields[i])
    fields = uniquify(fields)

    return lr.make_light_ray(start_position=start_position,
                             end_position=end_position,
                             trajectory=trajectory,
                             fields=fields,
                             setup_function=setup_function,
                             solution_filename=solution_filename,
                             data_filename=data_filename,
                             field_parameters=field_parameters,
                             redshift=redshift)

def make_compound_ray(parameter_filename, simulation_type,
                      near_redshift, far_redshift,
                      lines=None, ftype='gas', fields=None,
                      solution_filename=None, data_filename=None,
                      use_minimum_datasets=True, max_box_fraction=1.0,
                      deltaz_min=0.0, minimum_coherent_box_fraction=0.0,
                      seed=None, setup_function=None, load_kwargs=None,
                      line_database=None, ionization_table=None,
                      field_parameters = None):
    """
    Create a yt LightRay object for multiple consecutive datasets (eg IGM).
    This is a wrapper function around yt's LightRay interface to reduce some
    of the complexity there.

    .. note::

        The compound ray functionality has only been implemented for the
        Enzo and Gadget code.  If you would like to help us implement
        this functionality for your simulation code, please contact us
        about this on the mailing list.

    A compound ray is a series of straight lines passing through multiple
    consecutive outputs from a single cosmological simulation to approximate
    a continuous line of sight to high redshift.

    Because a single continuous ray traversing a simulated volume can only
    cover a small range in redshift space (e.g. 100 Mpc only covers the
    redshift range from z=0 to z=0.023), the compound ray passes rays through
    multiple consecutive outputs from the same simulation to approximate the
    path of a single line of sight to high redshift.  By probing all of the
    foreground material out to any given redshift, the compound ray is
    appropriate for studies of the intergalactic medium and circumgalactic
    medium.

    By default, it selects a random starting location and trajectory in
    each dataset it traverses, to assure that the same cosmological structures
    are not being probed multiple times from the same direction.  In doing
    this, the ray becomes discontinuous across each dataset.

    The compound ray requires the parameter_filename of the simulation run.
    This is *not* the dataset filename from a single output, but the parameter
    file that was used to run the simulation itself.  It is in this parameter
    file that the output frequency, simulation volume, and cosmological
    parameters are described to assure full redshift coverage can be achieved
    for a compound ray.  It also requires the simulation_type of the simulation.

    Unlike the simple ray, which is specified by its start and end positions
    in the dataset volume, the compound ray requires the near_redshift and
    far_redshift to determine which datasets to use to get full coverage
    in redshift space as the ray propagates from near_redshift to far_redshift.

    Like the simple ray produced by :class:`~trident.make_simple_ray`,
    each gas cell intersected by the LightRay is sampled for the desired
    fields and stored.  Several additional fields are created and stored
    including ``dl`` to represent the path length in space
    for each element in the ray, ``v_los`` to represent the line of
    sight velocity along the ray, and ``redshift``, ``redshift_dopp``, and
    ``redshift_eff`` to represent the cosmological redshift, doppler redshift
    and effective redshift (combined doppler and cosmological) for each
    element of the ray.

    The :lines: keyword can be set to automatically add all fields to the
    resulting ray necessary for later use with the SpectrumGenerator class.
    If the necessary fields do not exist for your line of choice, they will
    be added to your datasets before adding them to the ray.

    If using the :lines: keyword with SPH datasets, it is very important
    to set the :ftype: keyword appropriately, or you may end up calculating
    ion fields by interpolating on data already smoothed to the grid.  This is
    generally not desired.

    **Parameters**

    :parameter_filename: string

        The simulation parameter file *not* the dataset filename

    :simulation_type: string

        The simulation type of the parameter file.  At present, this
        functionality only works with "Enzo" and "Gadget" yt frontends.

    :near_redshift, far_redshift: floats

        The near and far redshift bounds of the LightRay through the
        simulation datasets.

    :lines: list of strings, optional

        List of strings that determine which fields will be added to the ray
        to support line deposition to an absorption line spectrum.  List can
        include things like "C", "O VI", or "Mg II ####", where #### would be
        the integer wavelength value of the desired line.  If set to 'all',
        includes all possible ions from H to Zn. :lines: can be used
        in conjunction with :fields: as they will not override each other.
        If using the :lines: keyword with an SPH dataset, it is very important
        to set the :ftype: keyword appropriately, or you may end up calculating
        ion fields by interpolating on data already smoothed to the grid.
        This is generally not desired.
        Default: None

    :ftype: string, optional

        For use with the :lines: keyword.  It is the field type of the fields to
        be added.  It is the first string in the field tuple e.g. "gas" in
        ("gas", "O_p5_number_density"). For SPH datasets, it is important to
        set this to the field type of the gas particles in your dataset
        (e.g. 'PartType0'), as it determines the source data for the ion
        fields to be added. If you leave it set to "gas", it will calculate
        the ion fields based on the hydro fields already smoothed on the grid,
        which is usually not desired.
        Default: "gas"

    :fields: list of strings, optional

        The list of which fields to store in the output LightRay.
        See :lines: keyword for additional functionality that will add fields
        necessary for creating absorption line spectra for certain line
        features.
        Default: None

    :solution_filename: string, optional

        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    :data_filename: string, optional

        Output filename for ray data stored as an HDF5 file.  Note that
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    :use_minimum_datasets: bool, optional

        Use the minimum number of datasets to make the ray continuous
        through the supplied datasets from the near_redshift to the
        far_redshift.  If false, the LightRay solution will contain as many
        datasets as possible to enable the light ray to traverse the
        desired redshift interval.
        Default: True

    :max_box_fraction: float, optional

        The maximum length a light ray segment can be in order to span the
        redshift interval from one dataset to another in units of the domain
        size.  Values larger than 1.0 will result in LightRays crossing the
        domain of a given dataset more than once, which is generally undesired.
        Zoom-in simulations can use a value equal to the length of the
        high-resolution region so as to limit ray segments to that size.  If
        the high-resolution region is not cubical, the smallest size should b
        used.
        Default: 1.0 (the size of the box)

    :deltaz_min: float, optional

        The minimum delta-redshift value between consecutive datasets used
        in the LightRay solution.
        Default: 0.0

    :minimum_coherent_box_fraction: float, optional

        When use_minimum_datasets is set to False, this parameter specifies
        the fraction of the total box width to be traversed before
        rerandomizing the ray location and trajectory.
        Default: 0.0

    :seed: int, optional

        Sets the seed for the random number generator used to determine the
        location and trajectory of the LightRay as it traverses the
        simulation datasets.  For consistent results between LightRays,
        use the same seed value.
        Default: None

    :setup_function: function, optional

        A function that will be called on the dataset as it is loaded but
        before the LightRay is generated.  Very useful for adding derived
        fields and other manipulations of the dataset prior to LightRay
        creation.
        Default: None

    :load_kwargs: dict, optional

        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    :line_database: string, optional

        For use with the :lines: keyword. If you want to limit the available
        ion fields to be added to those available in a particular subset,
        you can use a :class:`~trident.LineDatabase`.  This means when you
        set :lines:='all', it will only use those ions present in the
        corresponding LineDatabase.  If :LineDatabase: is set to None,
        and :lines:='all', it will add every ion of every element up to Zinc.
        Default: None

    :ionization_table: string, optional

        For use with the :lines: keyword.  Path to an appropriately formatted
        HDF5 table that can be used to compute the ion fraction as a function
        of density, temperature, metallicity, and redshift.  When set to None,
        it uses the table specified in ~/.trident/config
        Default: None

    :field_parameters: optional, dict
        Used to set field parameters in light rays. For example,
        if the 'bulk_velocity' field parameter is set, the relative
        velocities used to calculate peculiar velocity will be adjusted
        accordingly.
        Default: None.

    **Example**

    Generate a compound ray passing from the redshift 0 to redshift 0.05
    through a multi-output enzo simulation.

    >>> import trident
    >>> fn = 'path/to/simulation/parameter/file'
    >>> ray = trident.make_compound_ray(fn, simulation_type='Enzo',
    ... near_redshift=0.0, far_redshift=0.05, ftype='gas',
    ... lines=['H', 'O', 'Mg II'])

    Generate a compound ray passing from the redshift 0 to redshift 0.05
    through a multi-output gadget simulation.

    >>> import trident
    >>> fn = 'path/to/simulation/parameter/file'
    >>> ray = trident.make_compound_ray(fn, simulation_type='Gadget',
    ... near_redshift=0.0, far_redshift=0.05,
    ... lines=['H', 'O', 'Mg II'], ftype='PartType0')
    """
    if load_kwargs is None:
        load_kwargs = {}
    if fields is None:
        fields = []
    if data_filename is None:
        data_filename = 'ray.h5'

    lr = LightRay(parameter_filename,
                  simulation_type=simulation_type,
                  near_redshift=near_redshift,
                  far_redshift=far_redshift,
                  use_minimum_datasets=use_minimum_datasets,
                  max_box_fraction=max_box_fraction,
                  deltaz_min=deltaz_min,
                  minimum_coherent_box_fraction=minimum_coherent_box_fraction,
                  load_kwargs=load_kwargs)

    if ionization_table is None:
        ionization_table = ion_table_filepath

    # We use the final dataset from the simulation in order to test it for
    # what fields are present, etc.  This all assumes that the fields present
    # in this output will be present in ALL outputs.  Hopefully this is true,
    # because testing each dataset is going to be slow and a pain.

    sim = simulation(parameter_filename, simulation_type)
    ds = load(sim.all_outputs[-1]['filename'])

    # Include some default fields in the ray to assure it's processed correctly.

    fields = _add_default_fields(ds, fields)

    # If 'lines' kwarg is set, we need to get all the fields required to
    # create the desired absorption lines in the grid format, since grid-based
    # fields are what are directly probed by the LightRay object.

    # We first determine what fields are necessary for the desired lines, and
    # inspect the dataset to see if they already exist.  If so, we add them
    # to the field list for the ray.  If not, we have to create them.

    # We use the final dataset from the simulation in order to test it for
    # the sampling_type of the field, what fields are present, etc.  This
    # all assumes that the fields present in this output will be present in
    # ALL outputs.  Hopefully this is true, because testing each dataset
    # is going to be slow and a pain.

    if lines is not None:

        ion_list = _determine_ions_from_lines(line_database, lines)

        sampling_type = _check_sampling_types_match(ds, ftype)

        fields, fields_to_add_to_ds = _determine_fields_from_ions(ds, ion_list,
                                        fields, ftype, sampling_type)

        # actually add the fields we need to add to the dataset
        # by adding the fields to the setup_function passed to each dataset
        # as it is loaded by make_light_ray()

        def setup_ds(ds):

            for atom, ion_state in fields_to_add_to_ds:
                add_ion_number_density_field(atom, ion_state, ds,
                                            ftype=ftype,
                                            ionization_table=ionization_table,
                                            sampling_type=sampling_type)
            if setup_function is not None:
                setup_function(ds)

    else:
        # Define setup_ds in cases when no lines are specified.

        def setup_ds(ds):
            pass

    # To assure there are no fields that are double specified or that collide
    # based on being specified as "density" as well as ("gas", "density"),
    # we will just assume that all non-tuple fields requested are ftype "gas".
    for i in range(len(fields)):
        if isinstance(fields[i], str):
            fields[i] = ('gas', fields[i])
    fields = uniquify(fields)

    return lr.make_light_ray(seed=seed,
                             fields=fields,
                             setup_function=setup_ds,
                             solution_filename=solution_filename,
                             data_filename=data_filename,
                             redshift=None, njobs=-1,
                             field_parameters = field_parameters)

def _determine_ions_from_lines(line_database, lines):
    """
    Figure out what ions are necessary to produce the desired lines
    """
    if line_database is not None:
        line_database = LineDatabase(line_database)
        ion_list = line_database.parse_subset_to_ions(lines)
    else:
        ion_list = []
        if lines == 'all' or lines == ['all']:
            for k,v in six.iteritems(atomic_number):
                for j in range(v+1):
                    ion_list.append((k, j+1))
        else:
            for line in lines:
                linen = line.split()
                if len(linen) >= 2:
                    ion_list.append((linen[0], from_roman(linen[1])))
                elif len(linen) == 1:
                    num_states = atomic_number[linen[0]]
                    for j in range(num_states+1):
                        ion_list.append((linen[0], j+1))
                else:
                    raise RuntimeError("Cannot add a blank ion.")

    return uniquify(ion_list)

def _determine_fields_from_ions(ds, ion_list, fields, ftype, sampling_type):
    """
    Figure out what fields need to be added based on the ions present
    and the particle type

    For sph (determined by sampling_type):
    Identify if the number_density fields for the desired ions already
    exist on the dataset, and if so, just add them to the list
    of fields to be added to the ray.
    If not, add these ion fields to the dataset as particle fields,
    which prompts them being smoothed to the grid, and add the resulting
    smoothed fields (e.g. ("gas", "x_number_density")) to the list of fields
    to be added to the ray.  Include the 'temperature' field for
    calculating the width of voigt profiles in the absorption spectrum.

    For grid:
    Check if the number_density fields for these ions exist, and if so, add
    them to field list. If not, leave them off, as they'll be generated
    on the fly by SpectrumGenerator as long as we include the 'density',
    'temperature', and appropriate 'metallicity' fields.
    """
    fields_to_add_to_ds = []

    for ion in ion_list:
        atom = ion[0].capitalize()
        ion_state = ion[1]
        nuclei_field = "%s_nuclei_mass_density" % atom
        metallicity_field = "%s_metallicity" % atom
        field = "%s_p%d_number_density" % (atom, ion_state-1)

        if ion_state == 1:
            alias_field = "%s_number_density" % atom
        else:
            alias_field = field

        # This is ugly, but I couldn't find a way around it to hit
        # all 6 cases of when fields were present or not and particle
        # type or not.
        if (ftype, field) not in ds.derived_field_list:
            if (ftype, alias_field) not in ds.derived_field_list:
                if sampling_type == 'particle':
                    fields_to_add_to_ds.append((atom, ion_state))
                    fields.append(("gas", field))
                else:
                    # If this is a  grid-based field where the ion field
                    # doesn't yet exist, just append the density and
                    # appropriate metal field for ion field calculation
                    # on the ray itself instead of adding it to the full ds
                    fields.append(('gas', 'density'))
                    if ('gas', metallicity_field) in ds.derived_field_list:
                        fields.append(('gas', metallicity_field))
                    elif ('gas', nuclei_field) in ds.derived_field_list:
                        fields.append(('gas', nuclei_field))
                    elif atom != 'H':
                        fields.append(('gas', 'metallicity'))
                    else:
                        # Don't need metallicity field if we're just looking
                        # at hydrogen
                        pass

            else:
                fields.append(("gas", alias_field))
        else:
            fields.append(("gas", field))


    return fields, fields_to_add_to_ds

def _add_default_fields(ds, fields):
    """
    Add some default fields to rays to assure they can be processed correctly.
    """
    if ("gas", "temperature") in ds.derived_field_list:
        fields.append(("gas", 'temperature'))

    # H_nuclei_density should be added if possible to assure that the _log_nH
    # field, which is used as "density" in the ion_balance interpolation to
    # produce ion fields, is calculated as accurately as possible.
    if ('gas', 'H_nuclei_density') in ds.derived_field_list:
        fields.append(('gas', 'H_nuclei_density'))

    return fields
