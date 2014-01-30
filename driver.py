import numpy as np

from syn_spec_yt import SpectrumGenerator
from spectrum_plot import plot_spectrum

if __name__ == "__main__":
    lambda_min = 1150 # angstroms
    lambda_max = 1800
    dlambda    = 0.01
    n_lambda   = (lambda_max - lambda_min) / dlambda + 1
    sg = SpectrumGenerator(lambda_min, lambda_max, n_lambda)

    # if no filename given, will use a built-in file
    sg.load_line_list()

    sg.make_spectrum("../light_rays/mod_lightray_123456789.h5",
                     output_file="../spectrum.h5",
                     line_list_file="../line_list.out")

    # Add a composite qso spectrum at z = 0.5
    # sg.add_qso_spectrum(redshift=0.5)

    # Add gaussian noise
    # SNR = 30
    # sg.add_gaussian_noise(SNR)

    plot_spectrum(sg.lambda_bins, sg.flux_field, "../spectrum.pdf")
