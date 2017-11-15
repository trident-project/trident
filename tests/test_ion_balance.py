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
from trident.ion_balance import \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field
import yt
from yt.testing import \
    fake_random_ds, \
    fake_amr_ds, \
    fake_particle_ds
import tempfile
import shutil

import numpy as np

def test_add_ion_fraction_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_fraction_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_ion_fraction')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_number_density_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_number_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_density_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_mass_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_mass')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fraction_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_fraction_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_ion_fraction')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_number_density_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_number_density_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_number_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_density_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_density_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_mass_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'stream'
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds, ftype='stream')
    field = ('stream', 'O_p5_mass')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fields_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    ions = ['H', 'O', 'N V']
    tri.add_ion_fields(ds, ions, ftype='stream')
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_all_ion_fields_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    tri.add_ion_fields(ds, 'all', ftype='stream')
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_all_ion_fields_to_grid_ds_from_file():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'stream'
    ad = ds.all_data()
    tri.add_ion_fields(ds, 'all', ftype='stream', line_database='lines.txt')
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_all_ion_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'stream'
    ad = ds.all_data()
    ions = ['H', 'O', 'N V']
    tri.add_ion_fields(ds, ions, ftype='stream')
    fields = ['H_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fields_to_enzo():
    """
    Test to add various ion fields to Enzo dataset and slice on them
    """
    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    tri.add_ion_fields(ds, ['H', 'O VI'], ftype='gas')
    ad = ds.all_data()
    fields = ['H_p0_number_density', 'O_p5_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = ('gas', field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fields_to_gizmo():
    """
    Test to add various ion fields to gizmo dataset and slice on them
    """
    ds = yt.load('FIRE_M12i_ref11/snapshot_600.hdf5')
    tri.add_ion_fields(ds, ['H', 'O VI'], ftype='PartType0')
    ad = ds.all_data()
    fields = ['H_ion_fraction', 'O_p5_mass']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = ('gas', field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        yt.SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)
