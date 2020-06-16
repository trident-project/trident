"""
Tests for Utilities Code

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
import tempfile
import filecmp
from trident.utilities import \
    ensure_directory, \
    gunzip_file, \
    gzip_file, \
    make_onezone_dataset, \
    make_onezone_ray
from trident.spectrum_generator import \
    SpectrumGenerator
from trident.ray_generator import \
    make_simple_ray

def test_ensure_directory():
    """
    Tests ensure_directory code, which ensures a path exists and if not, 
    it creates such a path
    """
    tempdir = tempfile.mkdtemp()
    ensure_directory(tempdir)
    assert os.path.exists(tempdir)

    # is new random path created?
    subdir = os.path.join(tempdir, 'random')
    assert not os.path.exists(subdir)
    ensure_directory(subdir)
    assert os.path.exists(subdir)

def test_gzip_unzip_file():
    """
    Adds a test to ensure gzip and gunzip code is functionity properly by
    creating a random file, zipping it, unzipping it, and checking if the
    output is the same as the original file.
    """
    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()

    # Create a temporary file
    fp = tempfile.NamedTemporaryFile()
    fp.write(b'Hello world!')
    filepath = fp.name
    directory, filename = os.path.split(filepath)
    gzip_filename = os.path.join(tmpdir, filename+'.gz')
    gunzip_filename = os.path.join(tmpdir, filename)

    # Gzip it
    gzip_file(filepath, out_filename=gzip_filename, cleanup=False)

    # Gunzip it
    gunzip_file(gzip_filename, out_filename=gunzip_filename, cleanup=False)

    # Ensure contents of gunzipped file are the same as original
    filecmp.cmp(filepath, gunzip_filename)

def test_make_onezone_ray():
    """
    Tests the make_onezone_ray infrastructure by creating a ray and making a
    spectrum from it.
    """
    dirpath = tempfile.mkdtemp()
    ray_filename = os.path.join(dirpath, 'ray.h5')
    image_filename = os.path.join(dirpath, 'spec.png')
    ray = make_onezone_ray(column_densities={'H_p0_number_density':1e21},
                                   filename=ray_filename)    
    sg_final = SpectrumGenerator(lambda_min=1200, lambda_max=1300, dlambda=0.5)
    sg_final.make_spectrum(ray, lines=['Ly a'])
    sg_final.plot_spectrum(image_filename)

def test_make_onezone_dataset():
    """
    Tests the make_onezone_dataset infrastructure by generating a one_zone 
    dataset and then creating a ray and spectrum from it.
    """
    dirpath = tempfile.mkdtemp()
    ray_filename = os.path.join(dirpath, 'ray.h5')
    image_filename = os.path.join(dirpath, 'spec.png')
    ds = make_onezone_dataset()
    ray = make_simple_ray(ds, start_position=ds.domain_left_edge, 
                          end_position=ds.domain_right_edge, 
                          fields=['density', 'temperature', 'metallicity'],
                          data_filename=ray_filename)
    sg = SpectrumGenerator('COS') 
    sg.make_spectrum(ray)
    sg.plot_spectrum(image_filename)
