from matplotlib import pyplot
import numpy as np

from syn_spec_yt import SpectrumGenerator

if __name__ == "__main__":
    lambda_min = 800 # angstroms
    lambda_max = 2000
    n_lambda   = 10000
    sg = SpectrumGenerator(lambda_min, lambda_max, n_lambda)

    # if no filename given, will use a built-in file
    sg.load_line_list(filename=None)

    sg.make_spectrum("../light_rays/mod_lightray_123456800.h5",
                     output_file="../spectrum.h5",
                     line_list_file="../line_list.out")
    
    sg.add_qso_spectrum(redshift=0.5)
    sg.add_gaussian_noise(30)
    sg.flux_field.clip(0, sg.flux_field.max(), out=sg.flux_field)

    qso_spectrum = np.ones(n_lambda)
    sg.add_qso_spectrum(flux_field=qso_spectrum, redshift=0.5)
    pyplot.plot(sg.lambda_bins, sg.flux_field)
    pyplot.plot(sg.lambda_bins, qso_spectrum)
    pyplot.xlabel("wavelength [ang]")
    pyplot.ylabel("flux")
    pyplot.savefig("../spectrum.pdf")
