# Example of a currently working script using trident to generate a COS
# spectrum from a non-cosmological dataset

import os
import trident as tri

# Define the dataset and the coordinates of the start/end of the ray
fn = 'enzo_cosmology_plus/RD0009/RD0009'
ray_start = [0,0,0]
ray_end = [1,1,1]

# Make a LightRay object for the temperature, density, and metallicity
# fields of our dataset at redshift_start = redshift_end = 0.0.  Include 
# HI field calculated in simulation to override the Cloudy estimation of HI.  
# Save LightRay to ray.h5 and use it locally as ray object.
ray = tri.make_simple_ray(fn, start_position=ray_start,
                          end_position=ray_end, data_filename='ray.h5',
                          fields=['density', 'temperature', 'metallicity',
                          'H_p0_number_density'])

# Now use the ray object to actually generate an absorption spectrum
# Use the settings (spectral range, LSF, and spectral resolution) for COS
# And save it as an output hdf5 file and plot it to an image.
sg = tri.SpectrumGenerator('COS')
sg.make_spectrum(ray)
sg.save_spectrum('spec_raw.h5')
sg.plot_spectrum('spec_raw.png')

# "Final" spectrum with added quasar, MW background, and gaussian noise 
# (SNR=30)
sg.add_qso_spectrum()
sg.add_milky_way_foreground()
sg.apply_lsf()
sg.add_gaussian_noise(30)
sg.save_spectrum('spec_final.h5')
sg.plot_spectrum('spec_final.png')
