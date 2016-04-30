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
from yt.testing import \
    fake_random_ds, \
    fake_amr_ds, \
    fake_particle_ds

import numpy as np
# Make sure save occurs in tmpdir

def test_add_all_ion_fields_to_grid_ds():
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
    # Assure that a sampling of fields are added and can be sliced
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save()

def test_add_all_ion_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'stream'
    ad = ds.all_data()
    tri.add_ion_fields(ds, ftype='stream')
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save()

def test_add_all_ion_fields_to_particle_ds():
    """
    Test to add various ion fields
    """
    ds = fake_particle_ds(fields=('particle_mass',
                                  'particle_position_x',
                                  'particle_position_y',
                                  'particle_position_z',
                                  "particle_velocity_x",
                                  "particle_velocity_y",
                                  "particle_velocity_z",
                                  "density",
                                  "temperature",
                                  "metallicity",
                                  "smoothing_length"),
                           units=('g',
                                  'cm',
                                  'cm',
                                  'cm',
                                  'cm/s',
                                  'cm/s',
                                  'cm/s',
                                  'g/cm**3',
                                  'K',
                                  '',
                                  'cm'),
                           negative=(False,
                                     False,
                                     False,
                                     False,
                                     True,
                                     True,
                                     True,
                                     False,
                                     False,
                                     False,
                                     False))
    ftype = 'io',
    ad = ds.all_data()
    tri.add_ion_fields(ds, ftype='io')
    #len(ad[('gas', 'S_p3_ion_fraction')])
    #len(ad[('gas', 'H_mass')])
    #import pdb; pdb.set_trace()
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    #for field in fields:
    #    field = (ftype, field)
    #    assert field in ds.derived_field_list
    #    assert isinstance(ad[field], np.ndarray)
