"""
Ion fraction fields using Cloudy data.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Britton Smith.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.fields.local_fields import add_field
from yt.utilities.linear_interpolators import TrilinearFieldInterpolator, \
                                              UnilinearFieldInterpolator
from yt.utilities.physical_constants import mh
import numpy as na
import string
import h5py
import copy
import os
import pdb

H_mass_fraction = 0.76
to_nH = H_mass_fraction / mh

# set fractions to 0 for values lower than 1e-7, 
# which is what is used in Sutherland & Dopita (1993).
fraction_zero_point = 1.e-7
zero_out_value = -30.

table_store = {}

# Reads in comma separated ionization balance tables from 
# Sutherland & Dopita (1993).
class IonBalanceTable(object):
    def __init__(self, filename, atom=None):
        self.filename = filename
        self.parameters = []
        self.ion_fraction = []
        self._load_hdf5_table(atom)

    def _load_hdf5_table(self, atom):
        "Read in ion balance table from hdf5."
        input = h5py.File(self.filename, 'r')
        self.ion_fraction = input[atom].value
        self.ion_fraction[self.ion_fraction < na.log10(fraction_zero_point)] = zero_out_value
        for par in range(1, len(self.ion_fraction.shape) - 1):
            name = "Parameter%d" % par
            self.parameters.append(input[atom].attrs[name])
        self.parameters.append(input[atom].attrs['Temperature'])
        input.close()

def add_ion_fraction_field(atom, ion, model, data_file=None):
    """
    Add ion fraction field to yt.
    For example, add_ion_fraction_field('O',6) creates a field 
    called O_p5_"model"_ion_fraction.
    """
    atom = string.capitalize(atom)
    field = "%s_p%s_%s_ion_fraction" % (atom, ion-1, model)

    if data_file is None:
        data_file = "%s_ion_balance.h5" %model

    tableFile = "%s/tables/%s" % (os.path.dirname(__file__), data_file)

    if not table_store.has_key(field):
        ionTable = IonBalanceTable(tableFile, atom)
        table_store[field] = {'fraction': copy.deepcopy(ionTable.ion_fraction[ion-1]),
                              'parameters': copy.deepcopy(ionTable.parameters)}
        del ionTable

    add_field(field,function=_ion_fraction_field, units="")

def add_ion_number_density_field(atom, ion, model, **kwargs):
    """
    Add ion number density field to yt.
    For example, add_ion_number_density_field('O',6) creates a field 
    called O_p5_"model"_number_density.
    """
    atom = string.capitalize(atom)
    field = "%s_p%s_%s_number_density" % (atom, ion-1, model)
    add_ion_fraction_field(atom, ion, model, **kwargs)
    add_field(field,function=_ion_number_density,
              units="1.0/cm**3")

def add_ion_density_field(atom, ion, model, **kwargs):
    """
    Add ion mass density field to yt.
    For example, add_ion_density_field('O',6) creates a field 
    called O_p5_"model"_density.
    """
    atom = string.capitalize(atom)
    field = "%s_p%s_%s_density" % (atom, ion-1, model)
    add_ion_number_density_field(atom, ion, model, **kwargs)
    add_field(field,function=_ion_density,
              units="g/cm**3")

def add_ion_mass_field(atom, ion, model, **kwargs):
    """
    Add ion mass fields (g and Msun) to yt.
    For example, add_ion_density_field('O',6) creates a field 
    called O_p5_"model"_mass.
    """
    atom = string.capitalize(atom)
    field = "%s_p%s_%s_mass" % (atom, ion-1, model)
    add_ion_density_field(atom, ion, model, **kwargs)
    add_field(field,function=_ion_mass, units=r"g")

def _ion_mass(field,data):
    atom = field.name[1].split("_")[0]
    prefix = field.name[1].split("_mass")[0]
    densityField = "%s_density" % prefix
    return data[densityField] * data['cell_volume']

def _ion_density(field,data):
    atom = field.name[1].split("_")[0]
    prefix = field.name[1].split("_density")[0]
    numberDensityField = "%s_number_density" % prefix
    # the "mh" makes sure that the units work out
    return atomicMass[atom] * data[numberDensityField] * mh

def _ion_number_density(field,data):
    atom = field.name[1].split("_")[0]
    prefix = field.name[1].split("_number_density")[0]
    fractionField = "%s_ion_fraction" % prefix
    if atom == 'H' or atom == 'He':
        field = solarAbundance[atom] * data[fractionField] * \
                data['density']
    else:    
        field = data.ds.quan(solarAbundance[atom], "1.0/Zsun") * \
                data[fractionField] * data['metallicity'] * \
                data['density']
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

    def _log_nH(field, data):
        return na.log10(data['density'] * to_nH)
    if 'log_nH' not in data.ds.derived_field_list:
        data.ds.add_field('log_nH', function=_log_nH, units="")

    def _redshift(field, data):
        return data.ds.current_redshift * \
            na.ones(data['density'].shape, dtype=data['density'].dtype)
    if 'redshift' not in data.ds.derived_field_list:
        data.ds.add_field('redshift', function=_redshift, units="")

    def _log_T(field, data):
        return na.log10(data['temperature'])
    if 'log_T' not in data.ds.derived_field_list:
        data.ds.add_field('log_T', function=_log_T, units="")

    n_parameters = len(table_store[field.name[1]]['parameters'])

    if n_parameters == 1:
        ionFraction = table_store[field.name[1]]['fraction']
        t_param = table_store[field.name[1]]['parameters'][0]
        bds = (t_param[0], t_param[-1])
        
        interp = UnilinearFieldInterpolator(ionFraction, bds, 'log_T', truncate=True)

    elif n_parameters == 3:
        ionFraction = table_store[field.name[1]]['fraction']
        n_param = table_store[field.name[1]]['parameters'][0]
        z_param = table_store[field.name[1]]['parameters'][1]
        t_param = table_store[field.name[1]]['parameters'][2]
        bds = na.array([n_param[0], n_param[-1], z_param[0], z_param[-1], 
                    t_param[0], t_param[-1]])

        interp = TrilinearFieldInterpolator(ionFraction, bds, 
                                            ['log_nH', 'redshift', 'log_T'],
                                            truncate=True)

    else:
        raise RuntimeError("This data file format is not supported.")

    fraction = na.power(10, interp(data))
    fraction[fraction <= fraction_zero_point] = 0.0
    if (fraction > 1.0).any():
        print "WARNING! An ion fraction greater than 1 was calculated.  This is wrong!"
    return fraction

# Taken from Cloudy documentation.
solarAbundance = {'H':1.00e+00,'He':1.00e-01,'Li':2.04e-09,
                  'Be':2.63e-11,'B':6.17e-10,'C':2.45e-04,
                  'N':8.51e-05,'O':4.90e-04,'F':3.02e-08,
                  'Ne':1.00e-04,'Na':2.14e-06,'Mg':3.47e-05,
                  'Al':2.95e-06,'Si':3.47e-05,'P':3.20e-07,
                  'S':1.84e-05,'Cl':1.91e-07,'Ar':2.51e-06,
                  'K':1.32e-07,'Ca':2.29e-06,'Sc':1.48e-09,
                  'Ti':1.05e-07,'V':1.00e-08,'Cr':4.68e-07,
                  'Mn':2.88e-07,'Fe':2.82e-05,'Co':8.32e-08,
                  'Ni':1.78e-06,'Cu':1.62e-08,'Zn':3.98e-08}

atomicMass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941,
              'Be': 9.012182, 'B': 10.811, 'C': 12.0107,
              'N': 14.0067, 'O': 15.9994, 'F': 18.9984032,
              'Ne': 20.1797, 'Na': 22.989770, 'Mg': 24.3050,
              'Al': 26.981538, 'Si': 28.0855, 'P': 30.973761,
              'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
              'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955910,
              'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961,
              'Mn': 54.938049, 'Fe': 55.845, 'Co': 58.933200,
              'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409}
