"""
Tests for Spectrum Generation

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from trident.utilities import make_onezone_ray
from trident.spectrum_generator import \
    SpectrumGenerator, \
    load_spectrum
import tempfile
import shutil
import os

def test_save_load_spectrum_ascii():
    """
    Test that we can save and load spectra in the ASCII format
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(column_densities={'H_p0_number_density':1e21},
                                   filename=filename)
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines=['Ly a'])
    sg.save_spectrum(os.path.join(dirpath, 'spec.txt'))
    del(sg)
    sg = load_spectrum(os.path.join(dirpath, 'spec.txt'))
    sg.plot_spectrum(filename=os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
    assert True

def test_save_load_spectrum_hdf5():
    """
    Test that we can save and load spectra in the HDF5 format
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(column_densities={'H_p0_number_density':1e21},
                                   filename=filename)
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines=['Ly a'])
    sg.save_spectrum(os.path.join(dirpath, 'spec.h5'))
    del(sg)
    sg = load_spectrum(os.path.join(dirpath, 'spec.h5'))
    sg.plot_spectrum(filename=os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
    assert True

def test_save_load_spectrum_fits():
    """
    Test that we can save and load spectra in the FITS format
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(column_densities={'H_p0_number_density':1e21},
                                   filename=filename)
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines=['Ly a'])
    sg.save_spectrum(os.path.join(dirpath, 'spec.fits'))
    del(sg)
    sg = load_spectrum(os.path.join(dirpath, 'spec.fits'))
    sg.plot_spectrum(filename=os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
    assert True

def test_create_spectrum_all_lines():
    """
    Test that we can create a basic spectrum with all available lines
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(filename=filename)
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines='all')
    sg.plot_spectrum(os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
    assert True

def test_create_spectrum_H_O_lines():
    """
    Test that we can create a basic spectrum with hydrogen and oxygen lines
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(filename=filename)
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines=['H', 'O'])
    sg.plot_spectrum(os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
    assert True

def test_create_spectrum_H_lines_no_continuum():
    """
    Test that we can create a basic spectrum with H I lines but no Lyman
    continuum
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_onezone_ray(filename=filename)
    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg.make_spectrum(ray, lines=['H I'], ly_continuum=False)
    sg.plot_spectrum(os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
    assert True
