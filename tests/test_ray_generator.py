"""
Tests for ion balance code

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import trident as tri
import tempfile
import shutil
import os
from yt.loaders import \
    load
from trident import \
    LightRay, \
    make_simple_ray
from trident.testing import \
    answer_test_data_dir

COSMO_PLUS = os.path.join(answer_test_data_dir,
                          "enzo_cosmology_plus/AMRCosmology.enzo")
COSMO_PLUS_SINGLE = os.path.join(answer_test_data_dir,
                                 "enzo_cosmology_plus/RD0009/RD0009")
GIZMO_COSMO_SINGLE = os.path.join(answer_test_data_dir,
                                 "gizmo_cosmology_plus/snap_N128L16_150.hdf5")

def test_create_simple_grid_ray_with_lines():
    """
    Test to create a simple ray from a grid dataset using the lines kwargs
    Then generate a spectrum from it with those lines included.
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ds = tri.make_onezone_dataset()
    ray_start = ds.arr([0., 0., 0.], 'unitary')
    ray_end = ds.arr([1., 1., 1.], 'unitary')
    ray = tri.make_simple_ray(ds, start_position=ray_start, end_position=ray_end,
                              data_filename=filename,
                              lines=['H', 'C IV'], ftype='gas')
    sg = tri.SpectrumGenerator(lambda_min=1000, lambda_max=1400., dlambda=0.1)
    sg.make_spectrum(ray, lines=['H', 'C IV'])
    sg.plot_spectrum(os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)

def test_make_simple_ray_part():
    """
    Test the creation of a simple ray with default parameters for particle-
    based dataset
    """

    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ds = load(GIZMO_COSMO_SINGLE)
    ray = tri.make_simple_ray(ds, start_position=ds.domain_left_edge,
                              end_position=ds.domain_right_edge,
                              data_filename=filename,
                              lines=['Mg'])

def test_make_simple_ray_grid():
    """
    Test the creation of a simple ray with default parameters for grid-based
    dataset
    """

    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ds = load(COSMO_PLUS_SINGLE)
    ray = tri.make_simple_ray(ds, start_position=ds.domain_left_edge,
                              end_position=ds.domain_right_edge,
                              data_filename=filename,
                              lines=['Mg'])

def test_make_compound_ray_grid():
    """
    Test the creation of a compound ray with default parameters for grid-based
    dataset
    """

    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ray = tri.make_compound_ray(COSMO_PLUS, 'Enzo', 0.0, 0.03, lines=['O'],
                                fields=[('gas', 'temperature'),
                                        ('gas', 'metallicity')],
                                data_filename=filename)
