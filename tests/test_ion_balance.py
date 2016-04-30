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
import yt
from yt.testing import fake_random_ds
import numpy as np

def test_add_all_ion_fields():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(32, fields=("density", "velocity_x", "velocity_y",
                                    "velocity_z", "temperature", "metallicity"),
                            units= ('g/cm**3', 'cm/s', 'cm/s',
                                    'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    tri.add_ion_fields(ds, ftype='stream')
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
