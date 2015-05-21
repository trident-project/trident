# Example of a currently working script using trident to generate a COS 
# spectrum from a non-cosmological dataset

import trident as tri

# Load the dataset and define the coordinates of the start/end of the ray
fn = 'snap_r16n128vzw15p_098.bin'
bbox = [[-0.5, -0.5, -0.5],[0.5, 0.5, 0.5]]

#ds = yt.load(fn)
ray_start = [0,0,0]
ray_end = [1,1,1]

# Make a yt LightRay object for the temperature, density, and metallicity 
# fields of our dataset at redshift_start = redshift_end = 0.0.  Include true 
# HI to trump over Cloudy estimation of HI.  Save LightRay to ray.h5 and 
# ray.txt and use it internally as "ray"
lr = tri.LightRay(fn, load_kwargs={'bounding_box':bbox})
ray = lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                  fields=[('gas', 'temperature'), 
                          ('gas', 'density'), 
                          ('gas', 'H_nuclei_density'), 
                          ('gas', 'metallicity')],
                  solution_filename="ray.txt", data_filename="ray.h5",
                  get_los_velocity=True)

ray = tri.add_ion_fields_to_ray(ray, save_filename="mod_ray.h5")

# Now use the mod_ray h5 file to actually generate an absorption spectrum
# Set dlambda consistent with COS (need to do this with .set_instrument() 
# function based on spectral resolution R
lambda_min = 1200
lambda_max = 1400
dlambda = 0.01
n_lambda = (lambda_max - lambda_min) / dlambda + 1
sg = tri.SpectrumGenerator(lambda_min, lambda_max, n_lambda)
sg.make_spectrum("mod_ray.h5", output_file="spec.h5",
                line_list_file="line_list.txt")

# "Final" spectrum with added quasar, MW background, and gaussian noise (SNR=30)
sg.add_qso_spectrum(redshift=1.0)
sg.add_milky_way_foreground()
sg.apply_lsf('COS')
sg.add_gaussian_noise(10)
#sg.flux_field.clip(0, sg.flux_field.max(), out=sg.flux_field)
tri.plot_spectrum(sg.lambda_bins, sg.flux_field, 'spec.png')
