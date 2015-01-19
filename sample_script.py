import yt
import trident as tri

ds = yt.load('DD0100/DD0100')
ray_start = [0,0,0]
ray_end = [1,1,1]
spec = tri.make_spectrum(ds, ray_start, ray_end)
spec.add_qso_background(qso_redshift=2)
spec.add_mw_foreground()
spec.set_spectrograph('COS')
spec.add_noise(10)
spec.plot(lambda_min=1000, lambda_max=1800)
#spec.show()
spec.save()
