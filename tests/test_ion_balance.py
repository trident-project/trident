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

from trident.ion_balance import \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field, \
    add_ion_fields
from yt import \
    load, \
    SlicePlot
from yt.testing import \
    fake_random_ds, \
    fake_amr_ds
import tempfile
import shutil
from trident.testing import \
    answer_test_data_dir, \
    assert_array_rel_equal
import os

import numpy as np

ISO_GALAXY = os.path.join(answer_test_data_dir,
                'IsolatedGalaxy/galaxy0030/galaxy0030')
FIRE_SIM = os.path.join(answer_test_data_dir,
                'FIRE_M12i_ref11/snapshot_600.hdf5')

def test_add_ion_fraction_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ad = ds.all_data()
    add_ion_fraction_field('O', 6, ds)
    field = ('gas', 'O_p5_ion_fraction')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_number_density_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds)
    field = ('gas', 'O_p5_number_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_density_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds)
    field = ('gas', 'O_p5_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_mass_field_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds, ftype='gas')
    field = ('gas', 'O_p5_mass')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fraction_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ad = ds.all_data()
    add_ion_fraction_field('O', 6, ds)
    field = ('gas', 'O_p5_ion_fraction')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_number_density_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ad = ds.all_data()
    add_ion_number_density_field('O', 6, ds)
    field = ('gas', 'O_p5_number_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_density_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ad = ds.all_data()
    add_ion_density_field('O', 6, ds)
    field = ('gas', 'O_p5_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_mass_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ad = ds.all_data()
    add_ion_mass_field('O', 6, ds)
    field = ('gas', 'O_p5_mass')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    dirpath = tempfile.mkdtemp()
    SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fields_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'gas'
    ad = ds.all_data()
    ions = ['H', 'O', 'N V']
    add_ion_fields(ds, ions)
    fields = ['H_p0_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_all_ion_fields_to_grid_ds():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'gas'
    ad = ds.all_data()
    add_ion_fields(ds, 'all')
    fields = ['H_p0_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_all_ion_fields_to_grid_ds_from_file():
    """
    Test to add various ion fields
    """
    ds = fake_random_ds(8, fields=("density", "velocity_x", "velocity_y",
                                   "velocity_z", "temperature", "metallicity"),
                           units= ('g/cm**3', 'cm/s', 'cm/s',
                                   'cm/s', 'K', ''))
    ftype = 'gas'
    ad = ds.all_data()
    add_ion_fields(ds, 'all', ftype=ftype, line_database='lines.txt')
    fields = ['H_p0_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_all_ion_fields_to_amr_ds():
    """
    Test to add various ion fields
    """
    ds = fake_amr_ds(fields=("density", "velocity_x", "velocity_y",
                             "velocity_z", "temperature", "metallicity"))
    ftype = 'gas'
    ad = ds.all_data()
    ions = ['H', 'O', 'N V']
    add_ion_fields(ds, ions, ftype=ftype)
    fields = ['H_p0_ion_fraction', 'H_p0_number_density', 'O_p5_mass', 'N_p4_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = (ftype, field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fields_to_enzo():
    """
    Test to add various ion fields to Enzo dataset and slice on them
    """
    ds = load(ISO_GALAXY)
    add_ion_fields(ds, ['H', 'O VI'], ftype='gas')
    ad = ds.all_data()
    fields = ['H_p0_number_density', 'O_p5_density']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = ('gas', field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_add_ion_fields_to_gizmo():
    """
    Test to add various ion fields to gizmo dataset and slice on them
    """
    ds = load(FIRE_SIM)
    add_ion_fields(ds, ['H', 'O VI'], ftype='PartType0')
    ad = ds.all_data()
    fields = ['H_p0_ion_fraction', 'O_p5_mass']
    # Assure that a sampling of fields are added and can be sliced
    dirpath = tempfile.mkdtemp()
    for field in fields:
        field = ('gas', field)
        assert field in ds.derived_field_list
        assert isinstance(ad[field], np.ndarray)
        SlicePlot(ds, 'x', field).save(dirpath)
    shutil.rmtree(dirpath)

def test_ion_fraction_field_is_from_on_disk_fields():
    """
    Test to add various ion fields to Enzo dataset and slice on them
    """
    ds = load(ISO_GALAXY)
    add_ion_fields(ds, ['H'], ftype='gas')
    ad = ds.all_data()
    # Assure that a sampling of fields are added and can be sliced
    arr1 = ad['H_p0_ion_fraction']
    arr2 = ad['H_p0_number_density'] / ad['H_nuclei_density']
    assert_array_rel_equal(arr1, arr2, decimals=15)

def test_to_not_overwrite_fields_for_grid():
    """
    Test to not overwrite an existing ion field
    """
    ds = load(ISO_GALAXY)
    val_before = ds.r['H_p0_number_density'][0]
    add_ion_fields(ds, ['H'], ftype='gas')
    val_after = ds.r['H_p0_number_density'][0]
    assert val_before == val_after

def test_to_not_overwrite_fields_for_particle():
    """
    Test to not overwrite an existing ion field
    """
    ds = load(FIRE_SIM)
    val_sph_before = ds.r[('PartType0', 'H_p0_number_density')][0]
    val_gas_before = ds.r[('gas', 'H_p0_number_density')][0]
    add_ion_fields(ds, ['H'], ftype='PartType0')
    val_sph_after = ds.r[('PartType0', 'H_p0_number_density')][0]
    val_gas_after = ds.r[('gas', 'H_p0_number_density')][0]
    assert val_sph_before == val_sph_after
    assert val_gas_before == val_gas_after

