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
from yt.fields.local_fields import add_field
from yt.utilities.linear_interpolators import \
    TrilinearFieldInterpolator, \
    UnilinearFieldInterpolator
from yt.utilities.physical_constants import mh
from yt.funcs import mylog
import numpy as np
import string
import h5py
import copy
import os
from trident.line_database import \
    LineDatabase
from trident.utilities import \
    ion_table_filepath

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
            filename = ion_table_filepath()
        self.filename = filename
        self.parameters = []
        self.ion_fraction = []
        self._load_hdf5_table(atom)

    def _load_hdf5_table(self, atom):
        """
        Read in the HDF5 ion balance table
        """

        input = h5py.File(self.filename, 'r')
        self.ion_fraction = input[atom].value
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
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
    else:
        ftype = "gas"
    return np.log10(data[ftype, "density"] * to_nH)

def _redshift(field, data):
    """
    One index of ion balance table is in redshift, so this translates
    dataset's redshift values into the same format for indexing the table
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
    else:
        ftype = "gas"
    return data.ds.current_redshift * \
        np.ones(data[ftype, "density"].shape, dtype=data[ftype, "density"].dtype)

def _log_T(field, data):
    """
    One index of ion balance table is in log of temperature, so this translates
    dataset's temperature values into the same format for indexing the table
    """
    if isinstance(field.name, tuple):
        ftype = field.name[0]
    else:
        ftype = "gas"
    return np.log10(data[ftype, "temperature"])

def add_ion_fields(ds, ions=None, ftype='gas', 
                   ionization_table=None, 
                   field_suffix=False, 
                   line_database='lines.txt',
                   force_override=False):
    """
    Preferred method for adding ion fields to a yt dataset.

    Select ions based on the selection indexing set up in 
    :class:`~trident.LineDatabase.parse_subset_to_ions` function, that is, 
    by specifying a list of strings where each string represents an ion or
    line.  Selected ions will be a subset of the ions necessary to
    produce the lines present in the file specified in :line_database:, 
    nominally located in ``trident.__path__/data/line_lists``.

    For each ion species selected, 4 fields will be added (example for Mg II):

        * Ion fraction field. e.g. (ftype, 'Mg_p1_ion_fraction')
        * Number density field. e.g. (ftype, 'Mg_p1_number_density') 
        * Density field. e.g. (ftype, 'Mg_p1_density')
        * Mass field. e.g. (ftype, 'Mg_p1_mass')

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

    :ions: list of strings, optional

            List of strings matching possible lines.  Strings can be of the
            form:
            * Atom - Examples: "H", "C", "Mg"
            * Ion - Examples: "H I", "H II", "C IV", "Mg II"
            * Line - Examples: "H I 1216", "C II 1336", "Mg II 1240"
            * Identifier - Examples: "Ly a", "Ly b"

            If set to None, selects ions necessary to produce **all** lines
            in :class:`~trident.LineDatabase` specified by :line_database:.
            Default: None

    :ftype: string, optional

        The field type of the field to add.  it is the first string in the 
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field name that 
        indicates what ionization table was used.  Useful when using generating
        ion_fields that already exist in a dataset.

    :line_database: string, optional

        Ions are selected out of the set of ions present in the line_database
        constructed from the line list filename specified here.  See
        :class:`~trident.LineDatabase` for more information.

    :force_override: boolean, optional
        
        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    **Example**

    To add ionized hydrogen, doubly-ionized Carbon, and all of the Magnesium
    species fields to a dataset, you would run:

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_fields(ds, ions=['H II', 'C III', 'Mg'])
    """

    # Determine whether the user is trying to add a particle field 
    # based on the nature of other fields of that ftype in the dataset
    try:
        field_list_arr = np.asarray(ds.derived_field_list)
        mask = field_list_arr[:,0] == ftype
        valid_field = tuple(field_list_arr[mask][0])
        particle_type = ds.field_info[valid_field].particle_type
    except IndexError:
        raise RuntimeError('ftype %s not found in dataset %s' % (ftype, ds))

    if ionization_table is None:
        ionization_table = ion_table_filepath()

    # Parse the ions given following the LineDatabase syntax
    line_database = LineDatabase(line_database)
    ions = line_database.parse_subset_to_ions(ions)

    # adding X_p#_ion_mass field triggers the addition of:
    # - X_P#_ion_fraction
    # - X_P#_number_density
    # - X_P#_density
    for (atom, ion) in ions:
        add_ion_mass_field(atom, ion, ds, ftype, ionization_table,
            field_suffix=field_suffix, force_override=force_override)

def add_ion_fraction_field(atom, ion, ds, ftype="gas",
                           ionization_table=None,
                           field_suffix=False,
                           force_override=False):
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
        The field type of the field to add.  it is the first string in the 
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :ionization_table: string, optional
        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    :field_suffix: boolean, optional
        Determines whether or not to append a suffix to the field name that 
        indicates what ionization table was used

    :force_override: boolean, optional
        
        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    **Example**

    Add C IV (triply-ionized carbon) ion fraction field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_fraction_field('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_ion_fraction').save()
    """

    # Determine whether the user is trying to add a particle field 
    # based on the nature of other fields of that ftype in the dataset
    try:
        field_list_arr = np.asarray(ds.derived_field_list)
        mask = field_list_arr[:,0] == ftype
        valid_field = tuple(field_list_arr[mask][0])
        particle_type = ds.field_info[valid_field].particle_type
    except IndexError:
        raise RuntimeError('ftype %s not found in dataset %s' % (ftype, ds))

    if ionization_table is None:
        ionization_table = ion_table_filepath()

    if (ftype, "log_nH") not in ds.derived_field_list:
        ds.add_field((ftype, "log_nH"), function=_log_nH, units="",
                     particle_type=particle_type, 
                     force_override=force_override)

    if (ftype, "redshift") not in ds.derived_field_list:
        ds.add_field((ftype, "redshift"), function=_redshift, units="",
                     particle_type=particle_type,
                     force_override=force_override)

    if (ftype, "log_T") not in ds.derived_field_list:
        ds.add_field((ftype, "log_T"), function=_log_T, units="",
                     particle_type=particle_type,
                     force_override=force_override)

    atom = string.capitalize(atom)

    # if neutral ion field, alias X_p0_ion_fraction to X_ion_fraction field
    if ion == 1:
        field = "%s_ion_fraction" % atom
        alias_field = "%s_p0_ion_fraction" % atom
    else:
        field = "%s_p%d_ion_fraction" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
        if ion == 1:
            alias_field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]

    if not table_store.has_key(field):
        ionTable = IonBalanceTable(ionization_table, atom)
        table_store[field] = {'fraction': copy.deepcopy(ionTable.ion_fraction[ion-1]),
                              'parameters': copy.deepcopy(ionTable.parameters)}
        del ionTable

    ds.add_field((ftype, field), function=_ion_fraction_field, units="",
                 particle_type=particle_type, force_override=force_override)
    if ion == 1: # add aliased field too
        ds.field_info.alias((ftype, alias_field), (ftype, field))
        ds.derived_field_list.append((ftype, alias_field))

    # if ion particle field, add a smoothed deposited version to gas fields
    if particle_type:
        new_field = ds.add_smoothed_particle_field((ftype, field))
        if ftype != "gas":
            ds.field_info.alias(("gas", field), new_field)
            ds.derived_field_list.append(("gas", field))
            if ion == 1: # add aliased field too
                ds.field_info.alias(("gas", alias_field), new_field)
                ds.derived_field_list.append(("gas", alias_field))

def add_ion_number_density_field(atom, ion, ds, ftype="gas",
                                 ionization_table=None,
                                 field_suffix=False,
                                 force_override=False):
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

        The field type of the field to add.  it is the first string in the 
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used

    :force_override: boolean, optional
        
        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    **Example**

    Add C IV (triply-ionized carbon) number density field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_number_density('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_number_density').save()
    """
    # Determine whether the user is trying to add a particle field 
    # based on the nature of other fields of that ftype in the dataset
    try:
        field_list_arr = np.asarray(ds.derived_field_list)
        mask = field_list_arr[:,0] == ftype
        valid_field = tuple(field_list_arr[mask][0])
        particle_type = ds.field_info[valid_field].particle_type
    except IndexError:
        raise RuntimeError('ftype %s not found in dataset %s' % (ftype, ds))

    if ionization_table is None:
        ionization_table = ion_table_filepath()
    atom = string.capitalize(atom)
    # if neutral ion field, alias X_p0_number_density to X_number_density field
    if ion == 1:
        field = "%s_number_density" % atom
        alias_field = "%s_p0_number_density" % atom
    else:
        field = "%s_p%d_number_density" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
        if ion == 1:
            alias_field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]

    add_ion_fraction_field(atom, ion, ds, ftype, ionization_table,
                           field_suffix=field_suffix, 
                           force_override=force_override)
    ds.add_field((ftype, field),function=_ion_number_density,
                 units="cm**-3", particle_type=particle_type, 
                 force_override=force_override)
    if ion == 1: # add aliased field too
        ds.field_info.alias((ftype, alias_field), (ftype, field))
        ds.derived_field_list.append((ftype, alias_field))

    # if ion particle field, add a smoothed deposited version to gas fields
    if particle_type:
        new_field = ds.add_smoothed_particle_field((ftype, field))
        if ftype != "gas":
            ds.field_info.alias(("gas", field), new_field)
            ds.derived_field_list.append(("gas", field))
            if ion == 1: # add aliased field too
                ds.field_info.alias(("gas", alias_field), new_field)
                ds.derived_field_list.append(("gas", alias_field))

def add_ion_density_field(atom, ion, ds, ftype="gas",
                          ionization_table=None,
                          field_suffix=False,
                          force_override=False):
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

        The field type of the field to add.  it is the first string in the 
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used

    :force_override: boolean, optional
        
        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    **Example**

    Add C IV (triply-ionized carbon) mass density field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_density_field('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_density').save()
    """
    # Determine whether the user is trying to add a particle field 
    # based on the nature of other fields of that ftype in the dataset
    try:
        field_list_arr = np.asarray(ds.derived_field_list)
        mask = field_list_arr[:,0] == ftype
        valid_field = tuple(field_list_arr[mask][0])
        particle_type = ds.field_info[valid_field].particle_type
    except IndexError:
        raise RuntimeError('ftype %s not found in dataset %s' % (ftype, ds))

    if ionization_table is None:
        ionization_table = ion_table_filepath()
    atom = string.capitalize(atom)

    # if neutral ion field, alias X_p0_number_density to X_number_density field
    if ion == 1:
        field = "%s_density" % atom
        alias_field = "%s_p0_density" % atom
    else:
        field = "%s_p%d_density" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
        if ion == 1:
            alias_field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]

    add_ion_number_density_field(atom, ion, ds, ftype, ionization_table,
                                 field_suffix=field_suffix,
                                 force_override=force_override)
    ds.add_field((ftype, field), function=_ion_density,
                 units="g/cm**3", particle_type=particle_type,
                 force_override=force_override)
    if ion == 1: # add aliased field too
        ds.field_info.alias((ftype, alias_field), (ftype, field))
        ds.derived_field_list.append((ftype, alias_field))

    # if ion particle field, add a smoothed deposited version to gas fields
    if particle_type:
        new_field = ds.add_smoothed_particle_field((ftype, field))
        if ftype != "gas":
            ds.field_info.alias(("gas", field), new_field)
            ds.derived_field_list.append(("gas", field))
            if ion == 1: # add aliased field too
                ds.field_info.alias(("gas", alias_field), new_field)
                ds.derived_field_list.append(("gas", alias_field))

def add_ion_mass_field(atom, ion, ds, ftype="gas",
                       ionization_table=None,
                       field_suffix=False,
                       force_override=False):
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

        The field type of the field to add.  it is the first string in the 
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :ionization_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used

    :force_override: boolean, optional
        
        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    **Example**

    Add C IV (triply-ionized carbon) mass field to dataset

    >>> import yt
    >>> import trident
    >>> ds = yt.load('path/to/file')
    >>> trident.add_ion_mass_field('C', 4, ds)
    >>> yt.ProjectionPlot(ds, 'x', 'C_p3_mass').save()
    """
    # Determine whether the user is trying to add a particle field 
    # based on the nature of other fields of that ftype in the dataset
    try:
        field_list_arr = np.asarray(ds.derived_field_list)
        mask = field_list_arr[:,0] == ftype
        valid_field = tuple(field_list_arr[mask][0])
        particle_type = ds.field_info[valid_field].particle_type
    except IndexError:
        raise RuntimeError('ftype %s not found in dataset %s' % (ftype, ds))

    if ionization_table is None:
        ionization_table = ion_table_filepath()
    atom = string.capitalize(atom)

    # if neutral ion field, alias X_p0_number_density to X_number_density field
    if ion == 1:
        field = "%s_mass" % atom
        alias_field = "%s_p0_mass" % atom
    else:
        field = "%s_p%s_mass" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
        if ion == 1:
            alias_field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]

    add_ion_density_field(atom, ion, ds, ftype, ionization_table,
                          field_suffix=field_suffix,
                          force_override=force_override)
    ds.add_field((ftype, field), function=_ion_mass, units=r"g",
                 particle_type=particle_type,
                 force_override=force_override)
    if ion == 1: # add aliased field too
        ds.field_info.alias((ftype, alias_field), (ftype, field))
        ds.derived_field_list.append((ftype, alias_field))

    # if ion particle field, add a smoothed deposited version to gas fields
    if particle_type:
        new_field = ds.add_smoothed_particle_field((ftype, field))
        if ftype != "gas":
            ds.field_info.alias(("gas", field), new_field)
            ds.derived_field_list.append(("gas", field))
            if ion == 1: # add aliased field too
                ds.field_info.alias(("gas", alias_field), new_field)
                ds.derived_field_list.append(("gas", alias_field))

def _ion_mass(field, data):
    """
    Creates the function for a derived field to follow the total mass of an 
    ion over a dataset given that the specified ion's density field exists 
    in the dataset.
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
    density_field_name = "%s_density%s" % (prefix, suffix)
    if data.ds.field_info[(ftype, density_field_name)].particle_type:
        fraction_field_name = "%s_ion_fraction%s" % (prefix, suffix)
        return data[ftype, fraction_field_name] * \
          data[ftype, "particle_mass"]
    else:
        return data[ftype, density_field_name] * \
          data[ftype, "cell_volume"]

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

    # try the atom-specific density field first
    nuclei_field = "%s_nuclei_mass_density" % atom
    if (ftype, nuclei_field) in data.ds.field_info:
        return data[fraction_field_name] * \
          data[(ftype, nuclei_field)] / atomic_mass[atom] / mh

    # try the species metallicity
    metallicity_field = "%s_metallicity" % atom
    if (ftype, metallicity_field) in data.ds.field_info:
        return data[fraction_field_name] * \
                  data[ftype, "density"] * \
                  data[ftype, metallicity_field] / \
                  atomic_mass[atom] / mh

    if atom == 'H' or atom == 'He':
        field = solar_abundance[atom] * data[fraction_field_name] * \
                data[ftype, "density"]
    else:
        field = data.ds.quan(solar_abundance[atom], "1.0/Zsun") * \
          data[ftype, fraction_field_name] * \
          data[ftype, "metallicity"] * \
          data[ftype, "density"]
        # Ideally we'd like to use the following line
        # but it is very slow to compute.
        # If we get H_nuclei_density sped up
        # then we will want to remove the "to_nH" below
        # (this applies above as well)
        # data['H_nuclei_density']
    field[field <= 0.0] = 1.e-50
    # the "to_nH", does the final conversion to number density
    return field * to_nH

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
        mylog.warning("An ion fraction greater than 1 was calculated. ")
        mylog.warning("Bad interpolation: capping at 1. ")
        mylog.warning("Original values: %s" % fraction[greater_than])
        fraction = np.clip(fraction, 0.0, 1.0)
    return fraction

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
