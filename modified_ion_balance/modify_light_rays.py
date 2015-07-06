import h5py
import numpy as np
import glob as glob
from cloudy_ion_balance import *
import sys
from yt.units.yt_array import YTArray
import glob


def _add_ion_fields(ray, data_filename, ion_fields, verbose, save_filename):
    """
    Actually adds the ion fields from Cloudy to the LightRay object
    """
    atom_ion_count = {'H': 2, 'He': 3, 'Li': 4,
                      'Be': 5, 'B': 6, 'C': 6,
                      'N': 8, 'O': 9, 'F': 10,
                      'Ne': 11, 'Na': 12, 'Mg': 13,
                      'Al': 14, 'Si': 15, 'P': 16,
                      'S': 17, 'Cl': 18, 'Ar': 19,
                      'K': 20, 'Ca': 21, 'Sc': 22,
                      'Ti': 23, 'V': 24, 'Cr': 25,
                      'Mn': 26, 'Fe': 27}
    for key in atom_ion_count.keys():
        for i in range(atom_ion_count[key]):
            if verbose:
                print "Adding %s: %i of %i" %(key, i+1, atom_ion_count[key])
            add_Cloudy_ion_number_density_field(ray, key, i+1, data_file=data_filename)

    if save_filename is not None:
        if os.path.isfile(save_filename): os.remove(save_filename)
        out_file = h5py.File(save_filename, "w")
        for key in ray.keys():
            # If our field is a tuple, then just use the final element as the 
            # key name
            if isinstance(key, tuple):
                keyname = key[1]
            else:
                keyname = key
            out_file.create_dataset(keyname, dtype="float64", data=ray[key])
        out_file.close()

    return ray

def add_ion_fields_to_ray(ray, 
                          data_filename=None,
                          ion_fields=None, verbose=False, save_filename=None):
    """
    This function modifies a LightRay object to add in ion fields as tabulated
    in the specified data file.  If ion_fields=None, all ions are included.
    """
    return _add_ion_fields(ray, data_filename, ion_fields, verbose, \
                           save_filename)
    
def add_ion_fields_to_file(filename,
                           data_filename="../data/ion_balance/default.h5",
                           ion_fields=None, verbose=False, save_filename=None):
    """
    This function modifies a LightRay object contained in an hdf5 file to add 
    in ion fields as tabulated in the specified data file.  If ion_fields=None, 
    all ions are included.
    """
    ray_file = h5py.File(filename, "r")

    ray = {}
    for key in ray_file.keys():

        # Read in the data fields as YTArrays preserving the units output
        # in the hdf5 "units" attribute for each array
        # By default these are all in cgs, but we need metallicity in Zsun
        ray[key] = YTArray(ray_file[key].value, ray_file[key].attrs['units'])

    ray_file.close()

    return _add_ion_fields(ray, data_filename, ion_fields, verbose, \
                           save_filename)


# If this script is called directly on an hdf5 file containing a ray (or
# a collection of them), then add the cloudy ions to the ray and output
# as mod_filename.h5

if __name__ == '__main__':

    fns = glob.glob(sys.argv[1])
    
    for filename in fns:
        print "Modifying %s..." % filename
        mod_light_ray_file(filename, save_filename="mod_%s" % filename)
