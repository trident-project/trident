from matplotlib import pyplot
import numpy as np
import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'

from syn_spec_yt import SpectrumGenerator
from spectrum_plot import plot_spectrum

if __name__ == "__main__":
    lambda_min = 1150 # angstroms
    lambda_max = 1850
    dlambda    = 0.05
    n_lambda   = (lambda_max - lambda_min) / dlambda + 1
    sg = SpectrumGenerator(lambda_min, lambda_max, n_lambda)

    # if no filename given, will use a built-in file
    sg.load_line_list(filename=None)

    sg.make_spectrum("../camerons_lines/sightlines/mod_sightline_rho_0020_0000.h5",
                     output_file="../cam_spectrum.h5",
                     line_list_file="../cam_line_list.out")

    plot_spectrum(sg.lambda_bins, sg.flux_field, "../cam_spectrum_pure.png")
    
    # pyplot.plot(sg.lambda_bins, sg.flux_field)
    # pyplot.xlabel("wavelength [ang]")
    # pyplot.ylabel("flux")
    # pyplot.savefig("../cam_spectrum_pure.pdf")

    # Add a composite qso spectrum at z = 0.5
    sg.add_qso_spectrum(redshift=0.5)

    # Add Milky Way foreground
    sg.add_milky_way_foreground()

    plot_spectrum(sg.lambda_bins, sg.flux_field, "../cam_spectrum_nonoise.png")

    # Add gaussian noise
    # SNR = 30
    sg.add_gaussian_noise(30)
    
    sg.flux_field.clip(0, sg.flux_field.max(), out=sg.flux_field)

    #qso_spectrum = np.ones(n_lambda)
    #sg.add_qso_spectrum(flux_field=qso_spectrum, redshift=3.0)
    #pyplot.plot(sg.lambda_bins, qso_spectrum)

    plot_spectrum(sg.lambda_bins, sg.flux_field, "../cam_spectrum_qso.png")

    # pyplot.plot(sg.lambda_bins, sg.flux_field)
    # pyplot.xlabel("wavelength [ang]")
    # pyplot.ylabel("flux")
    # pyplot.savefig("../cam_spectrum_qso.pdf")
