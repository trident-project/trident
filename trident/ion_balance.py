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
from yt.utilities.linear_interpolators import TrilinearFieldInterpolator, \
                                              UnilinearFieldInterpolator
from yt.utilities.physical_constants import mh
from yt.funcs import mylog
import numpy as np
import string
import h5py
import copy
import os
from trident.utilities import \
    parse_config

H_mass_fraction = 0.76
to_nH = H_mass_fraction / mh

# set fractions to 0 for values lower than 1e-9,
# which is what is used in Sutherland & Dopita (1993).
fraction_zero_point = 1.e-9
zero_out_value = -30.

table_store = {}
ion_table_dir, ion_table_file = parse_config()
ion_table_filepath = os.path.join(ion_table_dir, ion_table_file)

class IonBalanceTable(object):
    def __init__(self, filename=None, atom=None):
        """
        IonBalanceTable class

        Used to load in an HDF5 file that contains the values for the 
        elemental ionization state of the gas as a function of density, 
        temperature, and metallcity.

        **Parameters**

        filename: string, optional
            Name of the HDF5 file that contains the ionization table data.  

            default: it uses the table specified in ~/.trident/config
 
        atom: string, optional
            The atomic species for which you want to create an IonBalanceTable 

            default: None
        """
        if filename is None:
            filename = ion_table_filepath
        self.filename = filename
        self.parameters = []
        self.ion_fraction = []
        self._load_hdf5_table(atom)

    def _load_hdf5_table(self, atom):
        "Read in ion balance table from hdf5."

        input = h5py.File(self.filename, 'r')
        self.ion_fraction = input[atom].value
        self.ion_fraction[self.ion_fraction < np.log10(fraction_zero_point)] = zero_out_value
        for par in range(1, len(self.ion_fraction.shape) - 1):
            name = "Parameter%d" % par
            self.parameters.append(input[atom].attrs[name])
        self.parameters.append(input[atom].attrs['Temperature'])
        input.close()

def _log_nH(field, data):
    return np.log10(data["gas", "density"] * to_nH)

def _redshift(field, data):
    return data.ds.current_redshift * \
        np.ones(data["gas", "density"].shape, dtype=data["gas", "density"].dtype)

def _log_T(field, data):
    return np.log10(data["gas", "temperature"])

def add_ion_fraction_field(atom, ion, ds, ionization_table=None,
                           field_suffix=False):
    """
    Add ion fraction field to a yt dataset for the desired ion.

    For example, add_ion_fraction_field('O', 6, ds) creates a field
    called O_p5_ion_fraction for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    **Parameters**

    atom: string
        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    ion : integer
        Ion number for desired species (e.g. 0 = neutral, 1 = singly ionized,
        2 = doubly ionized, etc.)

    ds : yt dataset object
        This is the dataset to which the ion fraction field will be added.

    ionization_table : string, optional
        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    field_suffix : boolean, optional
        Determines whether or not to append a suffix to the field name that 
        indicates what ionization table was used
    """

    if ionization_table is None:
        ionization_table = ion_table_filepath

    if ("gas", "log_nH") not in ds.derived_field_list:
        ds.add_field(("gas", "log_nH"), function=_log_nH, units="")

    if ("gas", "redshift") not in ds.derived_field_list:
        ds.add_field(("gas", "redshift"), function=_redshift, units="")

    if ("gas", "log_T") not in ds.derived_field_list:
        ds.add_field(("gas", "log_T"), function=_log_T, units="")

    atom = string.capitalize(atom)
    if ion == 1:
        field = "%s_ion_fraction" % atom
    else:
        field = "%s_p%d_ion_fraction" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]

    data_file = ionization_table

    if not table_store.has_key(field):
        ionTable = IonBalanceTable(data_file, atom)
        table_store[field] = {'fraction': copy.deepcopy(ionTable.ion_fraction[ion-1]),
                              'parameters': copy.deepcopy(ionTable.parameters)}
        del ionTable

    ds.add_field(("gas", field),function=_ion_fraction_field, units="")

def add_ion_number_density_field(atom, ion, ds, ionization_table=None,
                                 field_suffix=False):
    """
    Add ion number density field to a yt data object.

    For example, add_ion_fraction_field('O', 6, ds) creates a field
    called O_p5_ion_fraction for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    **Parameters**

    atom: string
        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    ion : integer
        Ion number for desired species (e.g. 0 = neutral, 1 = singly ionized,
        2 = doubly ionized, etc.)

    ds : yt dataset object
        This is the dataset to which the ion fraction field will be added.

    ionization_table : string, optional
        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    field_suffix : boolean, optional
        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used
    """
    if ionization_table is None:
        ionization_table = ion_table_filepath
    atom = string.capitalize(atom)
    if ion == 1:
        field = "%s_number_density" % atom
    else:
        field = "%s_p%d_number_density" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
    add_ion_fraction_field(atom, ion, ds, ionization_table, 
                           field_suffix=field_suffix)
    ds.add_field(("gas", field),function=_ion_number_density,
              units="1.0/cm**3")

def add_ion_density_field(atom, ion, ds, ionization_table=None,
                          field_suffix=False):
    """
    Add ion mass density field to a yt data object.

    For example, add_ion_fraction_field('O', 6, ds) creates a field
    called O_p5_ion_fraction for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    **Parameters**

    atom: string
        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    ion : integer
        Ion number for desired species (e.g. 0 = neutral, 1 = singly ionized,
        2 = doubly ionized, etc.)

    ds : yt dataset object
        This is the dataset to which the ion fraction field will be added.

    ionization_table : string, optional
        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    field_suffix : boolean, optional
        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used
    """
    if ionization_table is None:
        ionization_table = ion_table_filepath
    atom = string.capitalize(atom)
    if ion == 1:
        field = "%s_density" % atom
    else:
        field = "%s_p%d_density" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
    add_ion_number_density_field(atom, ion, ds, ionization_table,
                                 field_suffix=field_suffix)
    ds.add_field(("gas", field),function=_ion_density,
              units="g/cm**3")

def add_ion_mass_field(atom, ion, ds, ionization_table=None,
                       field_suffix=False):
    """
    Add ion mass fields (g and Msun) to a yt data object.

    For example, add_ion_fraction_field('O', 6, ds) creates a field
    called O_p5_ion_fraction for dataset ds, which represents 5-ionized
    oxygen (O plus 5 = O VI).

    **Parameters**

    atom: string
        Atomic species for desired ion fraction (e.g. 'H', 'C', 'Mg')

    ion : integer
        Ion number for desired species (e.g. 0 = neutral, 1 = singly ionized,
        2 = doubly ionized, etc.)

    ds : yt dataset object
        This is the dataset to which the ion fraction field will be added.
        will be added.

    ionization_table : string, optional
        Path to an appropriately formatted HDF5 table that can be used to 
        compute the ion fraction as a function of density, temperature, 
        metallicity, and redshift.  By default, it uses the table specified in
        ~/.trident/config
 
    field_suffix : boolean, optional
        Determines whether or not to append a suffix to the field
        name that indicates what ionization table was used
    """
    if ionization_table is None:
        ionization_table = ion_table_filepath
    atom = string.capitalize(atom)
    if ion == 1:
        field = "%s_mass" % atom
    else:
        field = "%s_p%s_mass" % (atom, ion-1)
    if field_suffix:
        field += "_%s" % ionization_table.split("/")[-1].split(".h5")[0]
    add_ion_density_field(atom, ion, ds, ionization_table,
                          field_suffix=field_suffix)
    ds.add_field(("gas", field),function=_ion_mass, units=r"g")

def _ion_mass(field, data):
    """
    Creates the function for a derived field to follow the total mass of an 
    ion over a dataset given that the specified ion's density field exists 
    in the dataset.
    """
    if isinstance(field.name, tuple):
        field_name = field.name[1]
    else:
        field_name = field.name
    atom = field_name.split("_")[0]
    prefix = field_name.split("_mass")[0]
    suffix = field_name.split("_mass")[-1]
    densityField = "%s_density%s" %(prefix, suffix)
    return data[densityField] * data['cell_volume']

def _ion_density(field,data):
    """
    Creates the function for a derived field for following the density of an 
    ion over a dataset given that the specified ion's number_density field 
    exists in the dataset.
    """
    if isinstance(field.name, tuple):
        field_name = field.name[1]
    else:
        field_name = field.name
    atom = field_name.split("_")[0]
    prefix = field_name.split("_density")[0]
    suffix = field_name.split("_density")[-1]
    numberDensityField = "%s_number_density%s" %(prefix, suffix)
    # the "mh" makes sure that the units work out
    return atomic_mass[atom] * data[numberDensityField] * mh

def _ion_number_density(field,data):
    """
    Creates the function for a derived field for following the number_density 
    of an ion over a dataset given that the specified ion's ion_fraction field 
    exists in the dataset.
    """
    if isinstance(field.name, tuple):
        field_name = field.name[1]
    else:
        field_name = field.name
    atom = field_name.split("_")[0]
    prefix = field_name.split("_number_density")[0]
    suffix = field_name.split("_number_density")[-1]
    fractionField = "%s_ion_fraction%s" %(prefix, suffix)
    if atom == 'H' or atom == 'He':
        field = solar_abundance[atom] * data[fractionField] * \
                data["gas", "density"]
    else:
        field = data.ds.quan(solar_abundance[atom], "1.0/Zsun") * \
                data[fractionField] * data["gas", "metallicity"] * \
                data["gas", "density"]
                # Ideally we'd like to use the following line
                # but it is very slow to compute.
                # If we get H_nuclei_density spread up
                # then we will want to remove the "to_nH" below
                # (this applies above as well)
                #data['H_nuclei_density']
    field[field <= 0.0] = 1.e-50
    # the "to_nH", does the final conversion to number density
    return field * to_nH

def _ion_fraction_field(field,data):
    """
    Creates the function for a derived field for following the ion_fraction
    of an ion over a dataset by plugging in the density, temperature, 
    metallicity and redshift of the output into the ionization table.
    """
    if isinstance(field.name, tuple):
        field_name = field.name[1]
    else:
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
                                            [("gas", "log_nH"),
                                             ("gas", "redshift"),
                                             ("gas", "log_T")],
                                            truncate=True)

    else:
        raise RuntimeError("This data file format is not supported.")

    fraction = np.power(10, interp(data))
    fraction[fraction <= fraction_zero_point] = 0.0
    if not isinstance(data, FieldDetector) and (fraction > 1.0).any():
        mylog.warning("An ion fraction greater than 1 was calculated.  " +
                      "This is wrong!")
    return fraction

# Taken from Cloudy documentation.
solar_abundance = {'H':1.00e+00,'He':1.00e-01,'Li':2.04e-09,
                   'Be':2.63e-11,'B':6.17e-10,'C':2.45e-04,
                   'N':8.51e-05,'O':4.90e-04,'F':3.02e-08,
                   'Ne':1.00e-04,'Na':2.14e-06,'Mg':3.47e-05,
                   'Al':2.95e-06,'Si':3.47e-05,'P':3.20e-07,
                   'S':1.84e-05,'Cl':1.91e-07,'Ar':2.51e-06,
                   'K':1.32e-07,'Ca':2.29e-06,'Sc':1.48e-09,
                   'Ti':1.05e-07,'V':1.00e-08,'Cr':4.68e-07,
                   'Mn':2.88e-07,'Fe':2.82e-05,'Co':8.32e-08,
                   'Ni':1.78e-06,'Cu':1.62e-08,'Zn':3.98e-08}

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
