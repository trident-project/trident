import h5py
import numpy as np
import glob as glob
from CloudyIonBalance import *
import sys
from yt.units.yt_array import YTArray


def add_ion_field(data):
    print ""


fns = [sys.argv[1]]

atom_ion_count = {'H': 2, 'He': 3, 'Li': 4,
                  'Be': 5, 'B': 6, 'C': 6,
                  'N': 8, 'O': 9, 'F': 10,
                  'Ne': 11, 'Na': 12, 'Mg': 13,
                  'Al': 14, 'Si': 15, 'P': 16,
                  'S': 17, 'Cl': 18, 'Ar': 19,
                  'K': 20, 'Ca': 21, 'Sc': 22,
                  'Ti': 23, 'V': 24, 'Cr': 25,
                  'Mn': 26, 'Fe': 27}

verbose = False

for file_name in fns:

    print "Modifying %s..." %file_name

    lr = h5py.File(file_name, "r")

    data = {}

    for key in lr.keys():

        # Read in the data fields as YTArrays preserving the units output
        # in the hdf5 "units" attribute for each array
        # By default these are all in cgs, but we need metallicity in Zsun
        data[key] = YTArray(lr[key].value, lr[key].attrs['units'])

    for key in atom_ion_count.keys():
        for i in range(atom_ion_count[key]):
            if verbose:
                print "%s: %i of %i" %(key, i+1, atom_ion_count[key])
            add_Cloudy_ion_number_density_field(data, key, i+1, data_file="cloudy_ion_balance.h5")

    lr.close()

    new_file = h5py.File("mod_%s" %file_name, "w")

    for key in data.keys():
        new_file.create_dataset(key, dtype="float64", data=data[key])

    new_file.close()
