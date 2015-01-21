# Example of a currently working script using trident to generate a COS 
# spectrum from a non-cosmological dataset

import yt
import trident as tri
import h5py as h5
from yt.analysis_modules.cosmological_observation.api import LightRay
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1 import host_subplot
import math
import numpy as np
import sys
from CloudyIonBalance import *


# Load the dataset and define the coordinates of the start/end of the ray
fn = 'enzo_cosmology_plus/RD0009/RD0009'
fn = '/Users/chummels/src/yt-data/enzo_cosmology_plus/RD0009/RD0009'
fn = 'enzo_cosmology_plus/AMRCosmology.enzo'
#ds = yt.load(fn)
ray_start = [0,0,0]
ray_end = [1,1,1]

# Make a yt LightRay object for the temperature, density, and metallicity 
# fields of our dataset at redshift_start = redshift_end = 0.0.  Include true 
# HI to trump over Cloudy estimation of HI.  Save LightRay to ray.h5 and 
# ray.txt
lr = LightRay(fn, 'Enzo', 0.0, 0.0)
lr.make_light_ray(seed=12345, start_position=ray_start, end_position=ray_end,
                  fields=['temperature', 'density', 'H_number_density', 'metallicity'],
                  #fields=['temperature', 'density', 'metallicity'],
                  solution_filename="ray.txt", data_filename="ray.h5",
                  get_los_velocity=True)

# Reload the light ray from the h5 file so we can now update it to incorporate
# Cloudy approximations of ions.
lr = h5.File("ray.h5")

# Update the yt LightRay object to incorporate ionic species
# Because of a bug in 3.0's LightRay, we have to div all the metallicity values
# by 0.0204 before continuing.  This should eventually go away.
metal = np.array(lr['metallicity'])
metal/=0.0204
del lr['metallicity']
lr['metallicity'] = metal

# The code below is applying the modified_ion_balance correction to the ray.h5
# file to create a mod_ray.h5 file which has all of the desired ion
# fields included for it.  i've added this code inline so as to create a single
# script without having to run things multiple times.  it requires a couple
# of files to be present (or in your pythonpath): CloudyIonBalance.py,
# and cloudy_ion_balance.h5
# Eventually this should get included in line as a "add_cloudy_fields()" 
# function. --CBH
atom_ion_count = {'H': 2, 'He': 3, 'Li': 4,
                  'Be': 5, 'B': 6, 'C': 6,
                  'N': 8, 'O': 9, 'F': 10,
                  'Ne': 11, 'Na': 12, 'Mg': 13,
                  'Al': 14, 'Si': 15, 'P': 16,
                  'S': 17, 'Cl': 18, 'Ar': 19,
                  'K': 20, 'Ca': 21, 'Sc': 22,
                  'Ti': 23, 'V': 24, 'Cr': 25,
                  'Mn': 26, 'Fe': 27}
data = {}
for key in lr.keys():
    data[key] = lr[key].value

for key in atom_ion_count.keys():
    for i in range(atom_ion_count[key]):
        #print "%s: %i of %i" %(key, i+1, atom_ion_count[key])
        add_Cloudy_ion_number_density_field(data, key, i+1, data_file="/Users/devinsilvia/Research/code/modified_ion_balance/cloudy_ion_balance.h5")
lr.close()
mod_fn = "mod_ray.h5"
if os.path.isfile(mod_fn): os.remove(mod_fn)
f = h5py.File(mod_fn)
for key in data.keys():
    f.create_dataset(key, dtype="float64", data=data[key])
f.close()

# Now use the mod_ray h5 file to actually generate an absorption spectrum
# Set dlambda consistent with COS (need to do this with .set_instrument() 
# function based on spectral resolution R
lambda_min = 1200
lambda_max = 1400
dlambda = 0.01
n_lambda = (lambda_max - lambda_min) / dlambda + 1
sg = tri.SpectrumGenerator(lambda_min, lambda_max, n_lambda)
sg.load_line_list(filename=None) 
## THINGS GET HERE BEFORE BREAKING.
sg.make_spectrum("mod_ray.h5", output_file="spec.h5",
                line_list_file="line_list.txt")

# "Final" spectrum with added quasar, MW background, and gaussian noise (SNR=30)
sg.add_qso_spectrum(redshift=1.0)
sg.add_milky_way_foreground()
sg.apply_lsf('COS')
sg.add_gaussian_noise(10)
#sg.flux_field.clip(0, sg.flux_field.max(), out=sg.flux_field)
tri.plot_spectrum(sg.lambda_bins, sg.flux_field, 'spec.png')
