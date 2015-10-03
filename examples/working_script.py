# Example of a currently working script using trident to generate a COS
# spectrum from a non-cosmological dataset

import os
import yt
import trident as tri

# Define the dataset and the coordinates of the start/end of the ray
fn = 'enzo_cosmology_plus/RD0009/RD0009'
ray_start = [0,0,0]
ray_end = [1,1,1]

# Make a yt LightRay object for the temperature, density, and metallicity
# fields of our dataset at redshift_start = redshift_end = 0.0.  Include true
# HI to override the Cloudy estimation of HI.  Save LightRay to ray.h5 and
# ray.txt and use it internally as "ray".  If it already exists, then load it
# as a dataset.

output_file = "ray.h5"
if os.path.exists(output_file):
    ray = yt.load(output_file)
else:
    lr = tri.LightRay(fn)
    ray = lr.make_light_ray(start_position=ray_start, end_position=ray_end, 
                solution_filename="ray.txt", data_filename=output_file,
                fields=['temperature', 'density', 'H_number_density', 
                'metallicity'])

# Now use the ray object to actually generate an absorption spectrum
# Use the settings (spectral range, LSF, and spectral resolution for COS)
sg = tri.SpectrumGenerator('COS')
sg.make_spectrum(ray, output_file="spec.h5")

# "Final" spectrum with added quasar, MW background, and gaussian noise 
# (SNR=30)
sg.add_qso_spectrum(redshift=0.0)
sg.add_milky_way_foreground()
sg.apply_lsf()
sg.add_gaussian_noise(30)
tri.plot_spectrum(sg.lambda_bins, sg.flux_field, 'spec.png')
