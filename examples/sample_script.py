# Example of how one would use trident to generate a COS spectrum
# from a non-cosmological dataset

import yt
import trident as tri

# Load the dataset and define the coordinates of the start/end of the ray
ds = yt.load('DD0100/DD0100')
ray_start = [0,0,0]
ray_end = [1,1,1]

# tri.make_spectrum
spec = tri.make_spectrum(ds, ray_start, ray_end)
spec.add_qso_background(qso_redshift=2)
spec.add_mw_foreground()
spec.set_spectrograph('COS-G140L')
spec.add_noise(10)

# Display the spectrum in an iPython window
#spec.show()    

# Save the spectrum
spec.save()
