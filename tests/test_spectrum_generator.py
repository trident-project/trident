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

from trident.plotting import plot_spectrum
from trident.utilities import make_onezone_ray
from trident.spectrum_generator import \
    SpectrumGenerator, \
    load_spectrum
import tempfile
import shutil
import os
from yt.loaders import \
    load
from trident import \
    make_simple_ray
from trident.testing import \
    answer_test_data_dir
import numpy as np

GIZMO_SINGLE = os.path.join(answer_test_data_dir,
                            "FIRE_M12i_ref11/snapshot_600.hdf5")
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

def test_input_types():
    """
    Test that spectra can be generated from ray file, dataset, or
    data container. Edited to use gizmo to generate ray. Will not work
    with onezone_ray framework.
    """

    dirpath = tempfile.mkdtemp()
    ds = load(GIZMO_SINGLE)
    ray_start = ds.domain_left_edge
    ray_end = ds.domain_right_edge
    filename = os.path.join(dirpath, 'ray.h5')
    ray = make_simple_ray(ds, start_position=ray_start, end_position=ray_end, data_filename=filename)

    sg = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    spectra = []

    # from file
    sg.make_spectrum(filename, lines=['H I'], ly_continuum=False)
    spectra.append(sg.flux_field[:])

    # from dataset
    sg.make_spectrum(ray, lines=['H I'], ly_continuum=False)
    spectra.append(sg.flux_field[:])
    assert (spectra[0] == spectra[1]).all()

    # from data container
    sg.make_spectrum(ray.all_data(), lines=['H I'], ly_continuum=False)
    spectra.append(sg.flux_field[:])
    assert (spectra[0] == spectra[2]).all()

    plot_spectrum(sg.lambda_field, spectra,
                  filename=os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)

def test_empty_ray_spectrum():

    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ds = load(GIZMO_SINGLE)
    ray_start = ds.domain_right_edge
    ray_end = ds.domain_right_edge - ds.arr([1, 1, 1], 'kpc')
    line_list = ['H', 'C', 'N', 'O', 'Mg']

    ray = make_simple_ray(ds, start_position=ray_start,
                          end_position=ray_end, data_filename=filename,
                          lines=line_list, fail_empty=False)

    sg = SpectrumGenerator('COS')
    sg.make_spectrum(ray, lines=line_list)
    # Should be a flat spectrum
    assert (sg.flux_field[:] == np.ones_like(sg.flux_field)).all()
    shutil.rmtree(dirpath)
