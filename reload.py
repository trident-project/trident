import numpy as np

from syn_spec_yt import SpectrumGenerator
from spectrum_plot import plot_spectrum

if __name__ == "__main__":
    lambda_min = 1150 # angstroms
    lambda_max = 1800
    dlambda    = 0.01
    n_lambda   = (lambda_max - lambda_min) / dlambda + 1
    sg = SpectrumGenerator(lambda_min, lambda_max, n_lambda)

    sg.load_spectrum("../spectrum.h5")
