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

from __future__ import absolute_import
import trident as tri
import yt
from yt.testing import \
    fake_random_ds, \
    fake_amr_ds, \
    fake_particle_ds
import tempfile
import shutil
import os
import numpy as np

def test_create_simple_grid_ray_with_lines():
    """
    Test to create a simple ray from a grid dataset using the lines kwargs
    Then generate a spectrum from it with those lines included.
    """
    dirpath = tempfile.mkdtemp()
    filename = os.path.join(dirpath, 'ray.h5')
    ds = tri.make_onezone_dataset()
    ftype = 'stream'
    ray_start = ds.arr([0., 0., 0.], 'unitary')
    ray_end = ds.arr([1., 1., 1.], 'unitary')
    ray = tri.make_simple_ray(ds, start_position=ray_start, end_position=ray_end,
                              data_filename=filename, 
                              lines=['H', 'C IV'], ftype='gas')
    sg = tri.SpectrumGenerator(lambda_min=1000, lambda_max=1400., dlambda=0.1)
    sg.make_spectrum(ray, lines=['H', 'C IV'])
    sg.plot_spectrum(os.path.join(dirpath, 'spec.png'))
    shutil.rmtree(dirpath)
