# Example of a currently working script using trident to generate a COS
# spectrum from a Gizmo dataset.  Because the sightline passes through the
# center of the galaxy, it should produce a DLA with lots of absorption lines.
# You can get this dataset at: # http://yt-project.org/data/
#
# wget http://yt-project.org/data/FIRE_M12i_ref11.tar.gz

import yt
import trident

# Set the dataset filename, load it into yt and define the trajectory
# of the LightRay.  This uses the maximum density location as the one end of
# the ray.  Define desired spectral features to include # all H, C, N, O, and
# Mg lines.
fn = 'FIRE_M12i_ref11/snapshot_600.hdf5'
ds = yt.load(fn)
_, c = ds.find_max(('gas', 'density'))
ray_start = c
ray_end = ds.domain_right_edge
line_list = ['H', 'C', 'N', 'O', 'Mg']

# Make a LightRay object including all necessary fields so you can add
# all H, C, N, O, and Mg fields to the resulting spectrum from your dataset.
# Save LightRay to ray.h5 and use it locally as ray object.
# Note: We use PartType0, the gas particle field type as our ftype!
ray = trident.make_simple_ray(ds, start_position=ray_start,
                              end_position=ray_end, data_filename='ray.h5',
                              lines=line_list, ftype='PartType0')

# Create a projection of the dataset in density along the x axis,
# overplot the trajectory of the ray, and save it.
# Note that this is using the 'gas' field type, which is the field deposited
# to the grid.
p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'))
p.annotate_ray(ray, arrow=True)
p.save('projection.png')

# Now use the ray object to actually generate an absorption spectrum
# Use the settings (spectral range, LSF, and spectral resolution) for COS
# And save it as an output text file and plot it to an image.
sg = trident.SpectrumGenerator('COS')
sg.make_spectrum(ray, lines=line_list)
sg.save_spectrum('spec_raw.txt')
sg.plot_spectrum('spec_raw.png')

# "Final" spectrum with added quasar, MW background, applied line-spread
# function, and added gaussian noise (SNR=30)
sg.add_qso_spectrum()
sg.add_milky_way_foreground()
sg.apply_lsf()
sg.add_gaussian_noise(30)
sg.save_spectrum('spec_final.txt')
sg.plot_spectrum('spec_final.png')
