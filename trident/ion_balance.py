"""
Ion fraction fields using Cloudy data.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_detector import \
    FieldDetector
from yt.utilities.linear_interpolators import \
    TrilinearFieldInterpolator, \
    UnilinearFieldInterpolator
from yt.utilities.physical_constants import mh
from yt.funcs import mylog
import numpy as np
import h5py
import copy
import os
from trident.config import \
    ion_table_filepath
from trident.line_database import \
    LineDatabase, \
    uniquify
from trident.roman import \
    from_roman

H_mass_fraction = 0.76
to_nH = H_mass_fraction / mh

# set fractions to 0 for values lower than 1e-9,
# which is what is used in Sutherland & Dopita (1993).
fraction_zero_point = 1.e-9
zero_out_value = -30.

table_store = {}

class IonBalanceTable(object):
    def __init__(self, filename=None, atom=None):
        """
        Base class for building additional ion fields

        Used to load in an HDF5 file that contains the values for the
        elemental ionization state of the gas as a function of density,
        temperature, metallcity, and redshift (for metagalactic photoionizing
        radiation field).

        **Parameters**

        :filename: string, optional

            Name of the HDF5 file that contains the ionization table data.

            Default: it uses the table specified in ~/.trident/config

        :atom: string, optional

            The atomic species for which you want to create an IonBalanceTable

            Default: None
        """
        if filename is None:
            filename = ion_table_filepath
        self.filename = filename
        self.parameters = []
        self.ion_fraction = []
        self._load_hdf5_table(atom)

    def _load_hdf5_table(self, atom):
        """
        Read in the HDF5 ion balance table
        """

        input = h5py.File(self.filename, 'r')
        self.ion_fraction = input[atom][()]
        self.ion_fraction[self.ion_fraction < np.log10(fraction_zero_point)] = zero_out_value
        for par in range(1, len(self.ion_fraction.shape) - 1):
            name = "Parameter%d" % par
            self.parameters.append(input[atom].attrs[name])
        self.parameters.append(input[atom].attrs['Temperature'])
        input.close()

def _log_nH(field, data):
    """
    One index of ion balance table is in log of density, so this translates
    dataset's density values into the same format for indexing the table

    N.B. All datasets *should* have an H_nuclei_density field defined if
    created in the standard way in yt.  Ray objects will also include
    H_nuclei_density from the parent dataset to assure identical behavior
    when ions are added to a ray as when they are added to an original dataset
    before then being included in a ray.
    """
    if ("gas", "H_nuclei_density") in data.ds.derived_field_list:
        log_nH_field = np.log10(data["gas", "H_nuclei_density"])
    else:
        log_nH_field = np.log10(data["gas", "density"] * to_nH)
    return log_nH_field

def _redshift(field, data):
    """
    One index of ion balance table is in redshift, so this translates
    dataset's redshift values into the same format for indexing the table

    Note that if there already exists a "redshift" field on the dataset (e.g.,
    on a LightRay dataset), that redshift field will be used instead.  This can
    lead to slight differences (1 part in 1e8) in the calculation of ion fields
    when added to a LightRay than when added to a dataset because redshift
    is continually varying (slightly) along the ray, whereas it is fixed for
    a standard dataset.
    """
    # Assure that redshift is defined for dataset--if not, assume z=0
    try:
        current_redshift = data.ds.current_redshift
    except AttributeError:
        current_redshift = 0.
    return current_redshift * \
        np.ones(data["gas", "density"].shape, dtype=data["gas", "density"].dtype)

def _log_T(field, data):
    """
    One index of ion balance table is in log of temperature, so this translates
    dataset's temperature values into the same format for indexing the table
    """
    return np.log10(data["gas", "temperature"])

def add_ion_fields(ds, ions, ftype='gas',
                   ionization_table=None,
                   field_suffix=False,
                   line_database=None,
                   sampling_type='local',
                   particle_type=None,
                   metal_source = 'best',
                   custom_metal_function = None):
    """
    Preferred method for adding ion fields to a yt dataset.

    Select ions based on the selection indexing set up in
    :class:`~trident.LineDatabase.parse_subset_to_ions` function, that is,
    by specifying a list of strings where each string represents an ion or
    line.  Strings are of one of three forms:

        * <element>
        * <element> <ion state>
        * <element> <ion state> <line_wavelength>

    If a line_database is selected, then the ions chosen will be a subset
    of the ions present in the equivalent :class:`~trident.LineDatabase`,
    nominally located in ``trident.__path__/data/line_lists``.

    For each ion species selected, four fields will be added (example for
    Mg II):

        * Ion fraction field. e.g. ("gas", 'Mg_p1_ion_fraction')
        * Number density field. e.g. ("gas", 'Mg_p1_number_density')
        * Density field. e.g. ("gas", 'Mg_p1_density')
        * Mass field. e.g. ("gas", 'Mg_p1_mass')

    This function is the preferred method for adding ion fields to one's
    dataset, but for more fine-grained control, one can also employ the
    :class:`~trident.add_ion_fraction_field`,
    :class:`~trident.add_ion_number_density_field`,
    :class:`~trident.add_ion_density_field`,
    :class:`~trident.add_ion_mass_field` functions individually.

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **Parameters**

    :ds: yt dataset object

        This is the dataset to which the ion fraction field will be added.

    :ions: list of strings

            List of strings matching possible lines.  Strings can be of the
            form:
            * Atom - Examples: "H", "C", "Mg"
            * Ion - Examples: "H I", "H II", "C IV", "Mg II"
            * Line - Examples: "H I 1216", "C II 1336", "Mg II 1240"

            If set to 'all', creates **all** ions for the first 30 elements:
            (ie hydrogen to zinc).  If set to 'all' with ``line_database``
            keyword set, then creates **all** ions associated with the lines
            specified in the equivalent :class:`~trident.LineDatabase`.

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        metallicity, and redshift.  When set to None, it uses the table
        specified in ~/.trident/config
        Default: None

    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field name that
        indicates what ionization table was used.  Useful when using generating
        ion_fields that already exist in a dataset.

    :line_database: string, optional

        Ions are selected out of the set of ions present in the line_database
        constructed from the line list filename specified here.  See
        :class:`~trident.LineDatabase` for more information.

    :ftype: string, optional

        This is deprecated and no longer necessary since all relevant
        fields are aliased to the 'gas' ftype.
        Default: 'gas'

    :sampling_type: string, optional

        This is deprecated and no longer necessary.
        Default: 'local'

    :particle_type: boolean, optional

        This is deprecated and no longer necessary.
        Default: 'auto'
        
    :metal_source: string, optional
    
        Specifies what source to use for metal nucleus
        abundances. Options are 'best', 'nuclei_mass_density',
        'nuclei_density', 'species_metallicity', 'metallicity',
        or 'custom'. 'best' finds most specific source available,
        all others besides 'custom' require that source, and throw
        FieldNotFound error if missing. 'custom' requires a
        custom_metal_function to be specified instead.
        Default: 'best'
        
    :custom_metal_function: func, optional
    
        if 'metal_source' is 'custom', this function will be used to
        determine metal abundances. Args are 'field' and 'data', like
        other function definitions (what ion is being created should be
        accessed from 'field').
        Default: None

    **Example**

    To add ionized hydrogen, doubly-ionized Carbon, and all of the Magnesium
    species fields to a dataset, you would run:

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_fields(ds, ions=['H II', 'C III', 'Mg'])
    """
    ion_list = []

    if ionization_table is None:
        ionization_table = ion_table_filepath

    # Parse the ions given following the LineDatabase syntax

    # If line_database is set, then use the underlying file as the line list
    # to select ions from.
    if line_database is not None:
        line_database = LineDatabase(line_database)
        ion_list = line_database.parse_subset_to_ions(ions)

    # Otherwise, any ion can be selected (not just ones in the line list).
    else:
        if ions == 'all' or ions == ['all']:
            for k, v in atomic_number.items():
                for j in range(v+1):
                    ion_list.append((k, j+1))
        else:
            for ion in ions:
                ionn = ion.split()
                if len(ionn) >= 2:
                    ion_list.append((ionn[0], from_roman(ionn[1])))
                elif len(ionn) == 1:
                    num_states = atomic_number[ionn[0]]
                    for j in range(num_states+1):
                        ion_list.append((ionn[0], j+1))
                else:
                    raise RuntimeError("Cannot add a blank ion.")

    # make sure ion list is unique
    ion_list = uniquify(ion_list)

    for (atom, ion) in ion_list:
        add_ion_fraction_field(atom, ion, ds, ftype, ionization_table,
            field_suffix=field_suffix, sampling_type=sampling_type,
            metal_source=metal_source, custom_metal_function=custom_metal_function)
        add_ion_number_density_field(atom, ion, ds, ftype, ionization_table,
            field_suffix=field_suffix, sampling_type=sampling_type,
            metal_source=metal_source, custom_metal_function=custom_metal_function)
        add_ion_density_field(atom, ion, ds, ftype, ionization_table,
            field_suffix=field_suffix, sampling_type=sampling_type,
            metal_source=metal_source, custom_metal_function=custom_metal_function)
        add_ion_mass_field(atom, ion, ds, ftype, ionization_table,
            field_suffix=field_suffix, sampling_type=sampling_type,
            metal_source=metal_source, custom_metal_function=custom_metal_function)

def add_ion_fraction_field(atom, ion, ds, ftype="gas",
                           ionization_table=None,
                           field_suffix=False,
                           sampling_type='local',
                           particle_type=None,
                           metal_source = 'best',
                           custom_metal_function = None):
    """
    Add ion fraction field to a yt dataset for the desired ion.

    .. note::

        The preferred method for adding ion fields to a dataset is using
        :class:`~trident.add_ion_fields`,

    For example, add_ion_fraction_field('O', 6, ds) creates a field
    called O_p5_ion_fraction for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI = 'O', 6).

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **Parameters**

    :atom: string
        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    :ion: integer
        Ion number for desired species (e.g. 1 = neutral, 2 = singly ionized,
        3 = doubly ionized, etc.)

    :ds: yt dataset object
        This is the dataset to which the ion fraction field will be added.

    :ftype: string, optional

        This is deprecated and no longer necessary since all relevant
        fields are aliased to the 'gas' ftype.
        Default: 'gas'

    :ionization_table: string, optional
        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config

    :field_suffix: boolean, optional
        Determines whether or not to append a suffix to the field name that
        indicates what ionization table was used

    :sampling_type: string, optional

        This is deprecated and no longer necessary.
        Default: 'local'

    :particle_type: boolean, optional

        This is deprecated and no longer necessary.
        Default: 'auto'
        
    :metal_source: string, optional
    
        In this function, this is deprecated and just used 
        for symmetry with the other calls.
        
    :custom_metal_function: func, optional
    
        In this function, this is deprecated and just used 
        for symmetry with the other calls.

    **Example**

    Add C IV (triply-ionized carbon) ion fraction field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_fraction_field('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_ion_fraction').save()
    """

    if ionization_table is None:
        ionization_table = ion_table_filepath

    if ("gas", "log_nH") not in ds.derived_field_list:
        _add_field(ds, ("gas", "log_nH"), function=_log_nH, units="",
                   sampling_type=sampling_type)

    if ("gas", "redshift") not in ds.derived_field_list:
        _add_field(ds, ("gas", "redshift"), function=_redshift, units="",
                   sampling_type=sampling_type)

    if ("gas", "log_T") not in ds.derived_field_list:
        _add_field(ds, ("gas", "log_T"), function=_log_T, units="",
                   sampling_type=sampling_type)

    atom = atom.capitalize()

    field = "%s_p%d_ion_fraction" % (atom, ion-1)

    if field_suffix:
        field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]

    if field not in table_store:
        ionTable = IonBalanceTable(ionization_table, atom)
        table_store[field] = {'fraction': copy.deepcopy(ionTable.ion_fraction[ion-1]),
                              'parameters': copy.deepcopy(ionTable.parameters)}
        del ionTable

    # if on-disk fields exist for calculation ion_fraction, use them
    if (("gas", "%s_p%d_number_density" % (atom, ion-1)) in ds.derived_field_list) and \
       (("gas", "%s_nuclei_density" % atom) in ds.derived_field_list):
        _add_field(ds, ("gas", field), function=_internal_ion_fraction_field, units="",
                   sampling_type=sampling_type)
    # otherwise, calculate ion_fraction from ion_balance table
    else:
        _add_field(ds, ("gas", field), function=_ion_fraction_field, units="",
                   sampling_type=sampling_type)

def add_ion_number_density_field(atom, ion, ds, ftype="gas",
                                 ionization_table=None,
                                 field_suffix=False,
                                 sampling_type='local',
                                 particle_type=None,
                                 metal_source = 'best',
                                 custom_metal_function = None):
    """
    Add ion number density field to a yt dataset for the desired ion.

    .. note::

        The preferred method for adding ion fields to a dataset is using
        :class:`~trident.add_ion_fields`,

    For example, add_ion_number_density_field('O', 6, ds) creates a field
    called O_p5_number_density for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **Parameters**

    :atom: string

        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    :ion: integer

        Ion number for desired species (e.g. 1 = neutral, 2 = singly ionized,
        3 = doubly ionized, etc.)

    :ds: yt dataset object

        This is the dataset to which the ion fraction field will be added.

    :ftype: string, optional

        This is deprecated and no longer necessary since all relevant
        fields are aliased to the 'gas' ftype.
        Default: 'gas'

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config

    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used

    :sampling_type: string, optional

        This is deprecated and no longer necessary.
        Default: 'local'

    :particle_type: boolean, optional

        This is deprecated and no longer necessary.
        Default: 'auto'

    :metal_source: string, optional
    
        Specifies what source to use for metal nucleus
        abundances. Options are 'best', 'nuclei_mass_density',
        'nuclei_density', 'species_metallicity', 'metallicity',
        or 'custom'. 'best' finds most specific source available,
        all others besides 'custom' require that source, and throw
        FieldNotFound error if missing. 'custom' requires a
        custom_metal_function to be specified instead.
        Default: 'best'
        
    :custom_metal_function: func, optional
    
        if 'metal_source' is 'custom', this function will be used to
        determine metal abundances. Args are 'field' and 'data', like
        other function definitions (what ion is being created should be
        accessed from 'field').
        Default: None

    **Example**

    Add C IV (triply-ionized carbon) number density field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_number_density('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_number_density').save()
    """

    if ionization_table is None:
        ionization_table = ion_table_filepath
    atom = atom.capitalize()

    field = "%s_p%d_number_density" % (atom, ion-1)
    fraction_field = "%s_p%d_ion_fraction" % (atom, ion-1)

    if field_suffix:
        field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]
        fraction_field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]

    if (ftype,fraction_field) not in ds.derived_field_list:
        mylog.info("fraction field %s not in dataset. Adding now."%fraction_field)
        add_ion_fraction_field(atom, ion, ds, ftype=ftype,
                                 ionization_table=ionization_table,
                                 field_suffix=field_suffix,
                                 sampling_type=sampling_type,
                                 particle_type=particle_type,
                                 metal_source = metal_source,
                                 custom_metal_function = custom_metal_function)
    
    _ion_number_density = _ion_number_density_wrapper(atom, ftype, ds, 
                                                      metal_source,
                                                      custom_metal_function)
    _add_field(ds, ("gas", field),function=_ion_number_density,
               units="cm**-3", sampling_type=sampling_type)

def add_ion_density_field(atom, ion, ds, ftype="gas",
                          ionization_table=None,
                          field_suffix=False,
                          sampling_type='local',
                          particle_type=None,
                          metal_source = 'best',
                          custom_metal_function = None):
    """
    Add ion mass density field to a yt dataset for the desired ion.

    .. note::

        The preferred method for adding ion fields to a dataset is using
        :class:`~trident.add_ion_fields`,

    For example, add_ion_density_field('O', 6, ds) creates a field
    called O_p5_density for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **Parameters**

    :atom: string

        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    :ion: integer

        Ion number for desired species (e.g. 1 = neutral, 2 = singly ionized,
        3 = doubly ionized, etc.)

    :ds: yt dataset object

        This is the dataset to which the ion fraction field will be added.

    :ftype: string, optional

        This is deprecated and no longer necessary since all relevant
        fields are aliased to the 'gas' ftype.
        Default: 'gas'

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config

    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used

    :sampling_type: string, optional

        This is deprecated and no longer necessary.
        Default: 'local'

    :particle_type: boolean, optional

        This is deprecated and no longer necessary.
        Default: 'auto'
        
    :metal_source: string, optional
    
        This is passed on to "add_ion_number_density" if that is called
        (which happens if the number density field doesn't already exist)
        
    :custom_metal_function: func, optional
    
        This is passed on to "add_ion_number_density" if that is called
        (which happens if the number density field doesn't already exist)

    **Example**

    Add C IV (triply-ionized carbon) mass density field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_density_field('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_density').save()
    """

    if ionization_table is None:
        ionization_table = ion_table_filepath
    atom = atom.capitalize()

    field = "%s_p%d_density" % (atom, ion-1)
    number_density_field = "%s_p%d_number_density" % (atom, ion-1)
    
    if field_suffix:
        field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]
        number_density_field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]
        
    if (ftype,number_density_field) not in ds.derived_field_list:
        mylog.info("number density field %s not in dataset. Adding now."%number_density_field)
        add_ion_number_density_field(atom, ion, ds, ftype=ftype,
                                 ionization_table=ionization_table,
                                 field_suffix=field_suffix,
                                 sampling_type=sampling_type,
                                 particle_type=particle_type,
                                 metal_source = metal_source,
                                 custom_metal_function = custom_metal_function)

    _add_field(ds, ("gas", field), function=_ion_density,
               units="g/cm**3", sampling_type=sampling_type)

def add_ion_mass_field(atom, ion, ds, ftype="gas",
                       ionization_table=None,
                       field_suffix=False,
                       sampling_type='local',
                       particle_type=None,
                       metal_source = 'best',
                       custom_metal_function = None):
    """
    Add ion mass field to a yt dataset for the desired ion.

    .. note::

        The preferred method for adding ion fields to a dataset is using
        :class:`~trident.add_ion_fields`,

    For example, add_ion_mass_field('O', 6, ds) creates a field
    called O_p5_mass for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **Parameters**

    :atom: string

        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    :ion: integer

        Ion number for desired species (e.g. 1 = neutral, 2 = singly ionized,
        3 = doubly ionized, etc.)

    :ds: yt dataset object

        This is the dataset to which the ion fraction field will be added.
        will be added.

    :ftype: string, optional

        This is deprecated and no longer necessary since all relevant
        fields are aliased to the 'gas' ftype.
        Default: 'gas'

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config

    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used

    :sampling_type: string, optional

        This is deprecated and no longer necessary.
        Default: 'local'

    :particle_type: boolean, optional

        This is deprecated and no longer necessary.
        Default: 'auto'
    
    :metal_source: string, optional
    
        This is passed on to "add_ion_density" if that is called
        (which happens if the density field doesn't already exist)
        
    :custom_metal_function: func, optional
    
        This is passed on to "add_ion_density" if that is called
        (which happens if the density field doesn't already exist)

    **Example**

    Add C IV (triply-ionized carbon) mass field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_mass_field('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_mass').save()
    """

    if ionization_table is None:
        ionization_table = ion_table_filepath
    atom = atom.capitalize()

    field = "%s_p%s_mass" % (atom, ion-1)
    density_field = "%s_p%s_density" % (atom, ion-1)
    
    if field_suffix:
        field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]
        density_field += "_%s" % ionization_table.split(os.sep)[-1].split(".h5")[0]
           
    if (ftype,density_field) not in ds.derived_field_list:
        mylog.info("density field %s not in dataset. Adding now."%density_field)
        add_ion_density_field(atom, ion, ds, ftype=ftype,
                                 ionization_table=ionization_table,
                                 field_suffix=field_suffix,
                                 sampling_type=sampling_type,
                                 particle_type=particle_type,
                                 metal_source = metal_source,
                                 custom_metal_function = custom_metal_function)

    _add_field(ds, ("gas", field), function=_ion_mass, units=r"g",
               sampling_type=sampling_type)

def _ion_mass(field, data):
    """
    Creates the function for a derived field for following the mass
    of an ion over a dataset given that the specified ion's ion_density field
    exists in the dataset.
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype = "gas"
        field_name = field.name
    atom = field_name.split("_")[0]
    prefix = field_name.split("_mass")[0]
    suffix = field_name.split("_mass")[-1]
    effective_volume = data[ftype,'mass']/data[ftype,'density']
    return data[ftype,'%s_density%s'%(prefix,suffix)]*effective_volume

def _ion_density(field, data):
    """
    Creates the function for a derived field for following the density of an
    ion over a dataset given that the specified ion's number_density field
    exists in the dataset.
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype = "gas"
        field_name = field.name
    atom = field_name.split("_")[0]
    prefix = field_name.split("_density")[0]
    suffix = field_name.split("_density")[-1]
    number_density_field_name = "%s_number_density%s" % (prefix, suffix)
    # the "mh" makes sure that the units work out
    return atomic_mass[atom] * data[ftype, number_density_field_name] * mh

def _determine_best_metal_source(atom, ftype, ds):
    """
    Decide on the most specific metal information available for 
    a particular atom. Returns a string which tells _ion_number_density
    what field to use.
    """
    if atom == 'H' or atom == 'He':
        return None
    if (ftype, "%s_nuclei_mass_density" % atom) in ds.derived_field_list:
        to_return = "nuclei_mass_density"
    elif (ftype, "%s_nuclei_density" % atom) in ds.derived_field_list:
        to_return = "nuclei_density"
    elif (ftype, "%s_metallicity" % atom) in ds.derived_field_list:
        to_return = "species_metallicity"
    else:
        assert (ftype, "metallicity") in ds.derived_field_list, \
                    "No usable metal source found in dataset."
        to_return = "metallicity"
    mylog.info("using %s as 'best' metal source"%to_return)
    return to_return

def _ion_number_density_wrapper(atom, ftype, ds, metal_source, custom_metal_function):
    """
    This wrapper exists so that ion_number_density can know what
    metal_source to use inside its local scope when it runs.
    """
    if metal_source == 'best':
        metal_source = _determine_best_metal_source(atom, ftype, ds)
    elif metal_source == 'custom':
        assert custom_metal_function is not None, \
            "custom metal function required if metal_source = 'custom'"
    def _ion_number_density(field, data):
        """
        Creates the function for a derived field for following the number_density
        of an ion over a dataset given that the specified ion's ion_fraction field
        exists in the dataset.
        """
        if isinstance(field.name, tuple):
            ftype = field.name[0]
            field_name = field.name[1]
        else:
            ftype = "gas"
            field_name = field.name
        atom = field_name.split("_")[0]
        prefix = field_name.split("_number_density")[0]
        suffix = field_name.split("_number_density")[-1]
        fraction_field_name = "%s_ion_fraction%s" % (prefix, suffix)
        
        if (ftype, "H_nuclei_density") in data.ds.derived_field_list:
            H_density = data[ftype, "H_nuclei_density"]
        else:
            H_density = data[ftype, "density"] * to_nH
        if atom == 'H' or atom == 'He':
            relative_number_density = solar_abundance[atom] * data[ftype, fraction_field_name]
            return relative_number_density * H_density
        
        mylog.info("generating metal fields using '%s' as metal_source"%metal_source)
        if metal_source == "nuclei_mass_density":
            nuclei_mass_density_field = "%s_nuclei_mass_density" % atom
            return data[ftype, fraction_field_name] * \
              data[(ftype, nuclei_mass_density_field)] / (atomic_mass[atom] * mh)
        if metal_source == "nuclei_density":
            nuclei_density_field = "%s_nuclei_density" % atom
            return data[ftype, fraction_field_name] * \
              data[(ftype, nuclei_density_field)]
        if metal_source == "species_metallicity":
            species_metallicity_field = "%s_metallicity" % atom
            return data[ftype, fraction_field_name] * \
              data[ftype, "density"] * \
              data[ftype, species_metallicity_field] / \
              (atomic_mass[atom] * mh)
        if metal_source == 'metallicity':
            relative_number_density = data.ds.quan(solar_abundance[atom], "1.0/Zsun") * \
              data[ftype, fraction_field_name] * \
              data[ftype, "metallicity"]
            return relative_number_density * H_density
        if metal_source == 'custom':
            return custom_metal_function(field,data)
        mylog.error('metal_source %s not understood.'%metal_source)
        
    return _ion_number_density

def _ion_fraction_field(field, data):
    """
    Creates the function for a derived field for following the ion_fraction
    of an ion over a dataset by plugging in the density, temperature,
    metallicity and redshift of the output into the ionization table.
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype = "gas"
        field_name = field.name
    n_parameters = len(table_store[field_name]['parameters'])

    if n_parameters == 1:
        ionFraction = table_store[field_name]['fraction']
        t_param = table_store[field_name]['parameters'][0]
        bds = t_param.astype("=f8")

        interp = UnilinearFieldInterpolator(ionFraction, bds, 'log_T', truncate=True)

    elif n_parameters == 3:
        ionFraction = table_store[field_name]['fraction']
        n_param = table_store[field_name]['parameters'][0]
        z_param = table_store[field_name]['parameters'][1]
        t_param = table_store[field_name]['parameters'][2]
        bds = [n_param.astype("=f8"), z_param.astype("=f8"), t_param.astype("=f8")]

        interp = TrilinearFieldInterpolator(ionFraction, bds,
                                            [(ftype, "log_nH"),
                                             (ftype, "redshift"),
                                             (ftype, "log_T")],
                                            truncate=True)

    else:
        raise RuntimeError("This data file format is not supported.")

    fraction = np.power(10, interp(data))
    fraction[fraction <= fraction_zero_point] = 0.0
    if not isinstance(data, FieldDetector) and (fraction > 1.0).any():
        greater_than = fraction > 1.0
        mylog.warning("%s > 1 was calculated. Capping values at 1." % field_name)
        mylog.warning("%d offenders: median = %f; maximum = %f" % (len(fraction[greater_than]), np.median(fraction[greater_than]), np.max(fraction[greater_than])))
        fraction = np.clip(fraction, 0.0, 1.0)
    return fraction

def _internal_ion_fraction_field(field, data):
    """
    Creates the function for a derived field for following the ion_fraction
    of an ion over a dataset when the number density and nuclei_density for
    that field are already defined on disk.  Example: H I ion fraction is
    H I number density / H number density from on-disk fields.
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype = "gas"
        field_name = field.name

    field_array = field_name.split('_')
    ion = '_'.join(field_array[:2])
    atom = field_array[0]
    return data[(ftype, "%s_number_density" % ion)] / \
            data[(ftype, "%s_nuclei_density" % atom)]


def _add_field(ds, name, function, units, sampling_type):
    """
    Private function for adding fields that wraps the yt add_field function.
    First it checks to see if field exists and if so, it does not attempt to
    add it and gives a warning message to user.  This avoids a misleading
    warning message from yt about using force_override=True, which no longer
    applies in Trident.
    """
    if name in ds.derived_field_list:
        mylog.warning("Field %s already exists. Not clobbering." % str(name))
        return
    else:
        return ds.add_field(name, function=function, units=units,
                            sampling_type=sampling_type)


def _alias_field(ds, alias_name, name):
    """
    Private function for aliasing fields that wraps the yt alias functionality.
    First it checks to see if alias field exists and if so, it does not attempt
    to alias it and gives a warning message to user.  This avoids extra
    adds of existing aliased fields.
    """
    if alias_name in ds.derived_field_list:
        mylog.warning("Field %s already exists. Not clobbering." % str(alias_name))
    else:
        ds.field_info.alias(alias_name, name)
        ds.derived_field_list.append(alias_name)
    return


# Taken from Cloudy documentation.
solar_abundance = {
    'H' : 1.00e+00, 'He': 1.00e-01, 'Li': 2.04e-09,
    'Be': 2.63e-11, 'B' : 6.17e-10, 'C' : 2.45e-04,
    'N' : 8.51e-05, 'O' : 4.90e-04, 'F' : 3.02e-08,
    'Ne': 1.00e-04, 'Na': 2.14e-06, 'Mg': 3.47e-05,
    'Al': 2.95e-06, 'Si': 3.47e-05, 'P' : 3.20e-07,
    'S' : 1.84e-05, 'Cl': 1.91e-07, 'Ar': 2.51e-06,
    'K' : 1.32e-07, 'Ca': 2.29e-06, 'Sc': 1.48e-09,
    'Ti': 1.05e-07, 'V' : 1.00e-08, 'Cr': 4.68e-07,
    'Mn': 2.88e-07, 'Fe': 2.82e-05, 'Co': 8.32e-08,
    'Ni': 1.78e-06, 'Cu': 1.62e-08, 'Zn': 3.98e-08}

atomic_mass = {
    'H' : 1.00794,   'He': 4.002602,  'Li': 6.941,
    'Be': 9.012182,  'B' : 10.811,    'C' : 12.0107,
    'N' : 14.0067,   'O' : 15.9994,   'F' : 18.9984032,
    'Ne': 20.1797,   'Na': 22.989770, 'Mg': 24.3050,
    'Al': 26.981538, 'Si': 28.0855,   'P' : 30.973761,
    'S' : 32.065,    'Cl': 35.453,    'Ar': 39.948,
    'K' : 39.0983,   'Ca': 40.078,    'Sc': 44.955910,
    'Ti': 47.867,    'V' : 50.9415,   'Cr': 51.9961,
    'Mn': 54.938049, 'Fe': 55.845,    'Co': 58.933200,
    'Ni': 58.6934,   'Cu': 63.546,    'Zn': 65.409}

atomic_number = {
    'H' : 1,  'He': 2,  'Li': 3,
    'Be': 4,  'B' : 5,  'C' : 6,
    'N' : 7,  'O' : 8,  'F' : 9,
    'Ne': 10, 'Na': 11, 'Mg': 12,
    'Al': 13, 'Si': 14, 'P' : 15,
    'S' : 16, 'Cl': 17, 'Ar': 18,
    'K' : 19, 'Ca': 20, 'Sc': 21,
    'Ti': 22, 'V' : 23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27,
    'Ni': 28, 'Cu': 29, 'Zn': 30}
