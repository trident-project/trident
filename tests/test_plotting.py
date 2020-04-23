"""
Tests for Plotting functionality

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from trident.plotting import plot_spectrum
from trident.utilities import make_onezone_ray
from trident.spectrum_generator import \
    SpectrumGenerator, \
    load_spectrum
import tempfile
import shutil
import os

def test_plot_multiple_spectra():
    """
    Test that multiple spectra can be plotted on top of each other. Example
    from plot_spectrum docstrings.
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(column_densities={'H_p0_number_density':1e21},
                                   filename=filename)
    sg_final = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg_final.make_spectrum(ray, lines=['Ly a'])
    sg_final.save_spectrum(os.path.join(dirpath, 'spec_raw.h5'))
    sg_final.add_gaussian_noise(10)
    sg_raw = load_spectrum(os.path.join(dirpath, 'spec_raw.h5'))
    plot_spectrum([sg_raw.lambda_field, sg_final.lambda_field],
                  [sg_raw.flux_field, sg_final.flux_field],
                  stagger=0, step=[False, True], label=['Raw', 'Noisy'],
                  filename=os.path.join(dirpath, 'raw_and_noise.png'))
    shutil.rmtree(dirpath)
    assert True
