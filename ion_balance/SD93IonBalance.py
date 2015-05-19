"""
Sutherland & Dopita (1993) ion fraction and number density field 
generator for yt.

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
from yt.utilities.linear_interpolators import UnilinearFieldInterpolator
from yt.utilities.physical_constants import mh
import numpy as na
import string
import roman
import h5py
import copy
import os

H_mass_fraction = 0.76
to_nH = H_mass_fraction / mh

# set fractions to 0 for values lower than 1e-7, 
# which is the lowest value in SD93.
fraction_zero_point = 1.e-7
zero_out_value = -30.

SD93_table_store = {}

# Reads in comma separated ionization balance tables from 
# Sutherland & Dopita (1993).
class SD93IonBalanceTable(object):
    def __init__(self, filename, atom=None):
        self.filename = filename
        self.temperature = []
        self.ion_fraction = []
        if atom is None:
            self._load_ascii_table()
        else:
            self._load_hdf5_table(atom)

    def _load_ascii_table(self):
        "Read in comma separated Sutherland & Dopita (1993) ion balance table."
        f = open(self.filename)
        for line in f.xreadlines():
            line.strip()
            if (line[0] != '#'):
                onLine = line.split(',')
                self.temperature.append(float(onLine.pop(0)))
                for q in range(len(onLine)):
                    onLine[q] = -1.0 * float(onLine[q])
                self.ion_fraction.append(onLine)
        self.temperature = na.array(self.temperature)
        self.ion_fraction = na.array(self.ion_fraction)
        self.ion_fraction = na.transpose(self.ion_fraction)
        f.close()

    def _load_hdf5_table(self, atom):
        "Read in ion balance table from hdf5."
        input = h5py.File(self.filename, 'r')
        self.ion_fraction = input[atom].value
        self.ion_fraction[self.ion_fraction < na.log10(fraction_zero_point)] = zero_out_value
        self.temperature = input[atom].attrs['Temperature']
        input.close()

def add_SD93_ion_fraction_field(atom,ion):
    """
    Add ion fraction field to yt.
    For example, add_SD93_ion_fraction_field('O',6) creates a field 
    called OVI_Ion_Fraction.
    """
    atom = string.capitalize(atom)
    field = "%s%s_SD93_eq_Ion_Fraction" % (atom,roman.toRoman(ion))

    tableFile = "%s/tables/SD93_ion_balance.h5" % os.path.dirname(__file__)

    if not SD93_table_store.has_key(field):
        ionTable = SD93IonBalanceTable(tableFile, atom=atom)
        SD93_table_store[field] = {'fraction': copy.deepcopy(ionTable.ion_fraction[ion-1]),
                                   'temperature': copy.deepcopy(ionTable.temperature),
                                   't_min': ionTable.temperature.min(),
                                   't_max': ionTable.temperature.max()}
        del ionTable

    add_field(field,function=_ion_fraction_field,units=r"")

def add_SD93_ion_number_density_field(atom,ion):
    """
    Add ion number density field to yt.
    For example, add_SD93_ion_number_density_field('O',6) creates a field 
    called OVI_NumberDensity.
    """
    atom = string.capitalize(atom)
    field = "%s%s_SD93_eq_NumberDensity" % (atom,roman.toRoman(ion))
    add_SD93_ion_fraction_field(atom,ion)
    add_field(field,function=_ion_number_density,
              units="1.0/cm**3")

def add_SD93_ion_density_field(atom,ion):
    """
    Add ion mass density field to yt.
    For example, add_SD93_ion_density_field('O',6) creates a field 
    called OVI_Density.
    """
    atom = string.capitalize(atom)
    field = "%s%s_SD93_eq_Density" % (atom,roman.toRoman(ion))
    add_SD93_ion_number_density_field(atom,ion)
    add_field(field,function=_ion_density,
              units="g/cm**3")

def add_SD93_ion_mass_field(atom,ion):
    """
    Add ion mass fields (g and Msun) to yt.
    For example, add_SD93_ion_density_field('O',6) creates a field 
    called OVI_Density.
    """
    atom = string.capitalize(atom)
    field = "%s%s_SD93_eq_Mass" % (atom,roman.toRoman(ion))
    field_msun = "%s%s_SD93_eq_MassMsun" % (atom,roman.toRoman(ion))
    add_SD93_ion_density_field(atom,ion)
    add_field(field,function=_ion_mass, units=r"g")
    add_field(field_msun,function=_ion_mass, 
              units="Msun")

def _ion_mass(field,data):
    species = field.name[1].split("_")[0]
    if species[1] == string.lower(species[1]):
        atom = species[0:2]
    else:
        atom = species[0]

    densityField = "%s_SD93_eq_Density" % species
    return data[densityField] * data['cell_volume']

def _ion_density(field,data):
    species = field.name[1].split("_")[0]
    if species[1] == string.lower(species[1]):
        atom = species[0:2]
    else:
        atom = species[0]

    numberDensityField = "%s_SD93_eq_NumberDensity" % species
    # the "mh" makes sure that the units work out
    return atomicMass[atom] * data[numberDensityField] * mh

def _ion_number_density(field,data):
    species = field.name[1].split("_")[0]
    if species[1] == string.lower(species[1]):
        atom = species[0:2]
    else:
        atom = species[0]

    fractionField = "%s_SD93_eq_Ion_Fraction" % species
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

def _ion_fraction_field(field, data):
    ionFraction = SD93_table_store[field.name[1]]['fraction']

    def _log_T(field, data):
        return na.log10(data['temperature'])
    data.ds.add_field('log_T', function=_log_T, units="")

    bds = (SD93_table_store[field.name[1]]['t_min'], 
           SD93_table_store[field.name[1]]['t_max'])

    interp = UnilinearFieldInterpolator(ionFraction, bds, 'log_T', truncate=True)

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
