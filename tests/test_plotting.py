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
from trident.spectrum_generator import SpectrumGenerator
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
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines=['Ly a'])
    sg.save_spectrum(os.path.join(dirpath, 'spec.h5'))
    sg.add_gaussian_noise(10)
    sg.save_spectrum(os.path.join(dirpath, 'noise.h5'))
    sg1 = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg2 = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg1.load_spectrum(os.path.join(dirpath, 'spec.h5'))
    sg2.load_spectrum(os.path.join(dirpath, 'noise.h5'))
    plot_spectrum([sg1.lambda_field, sg2.lambda_field], 
                  [sg1.flux_field, sg2.flux_field], 
                  stagger=0, step=[False, True], 
                  filename=os.path.join(dirpath, 'raw_and_noise.png'))
    shutil.rmtree(dirpath)
    assert True
