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
    add_ion_fields, \
    calculate_ion_fraction, \
    atomic_mass, \
    mh
from yt import \
    load, \
    SlicePlot
from yt.testing import \
    fake_random_ds, \
    fake_amr_ds
from yt.utilities.physical_constants import mh
import tempfile
import shutil
from trident.testing import \
    answer_test_data_dir, \
    assert_array_rel_equal
import os
import os.path as path

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
                             "velocity_z", "temperature", "metallicity"),
                     units=('g/cm**3', 'cm/s', 'cm/s', 'cm/s', 'K', ''))
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
                             "velocity_z", "temperature", "metallicity"),
                     units=('g/cm**3', 'cm/s', 'cm/s', 'cm/s', 'K', ''))
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
                             "velocity_z", "temperature", "metallicity"),
                     units=('g/cm**3', 'cm/s', 'cm/s', 'cm/s', 'K', ''))
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
                             "velocity_z", "temperature", "metallicity"),
                     units=('g/cm**3', 'cm/s', 'cm/s', 'cm/s', 'K', ''))
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
                             "velocity_z", "temperature", "metallicity"),
                     units=('g/cm**3', 'cm/s', 'cm/s', 'cm/s', 'K', ''))
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

def test_add_ion_fields_to_enzo_with_nonsolar_abundances():
    """
    Test adding a non-tracked metal to an Enzo dataset with nonsolar abundance
    """
    abun = {"O": 5e-3}

    ds = load(ISO_GALAXY)
    ad = ds.all_data()
    add_ion_number_density_field('O', 6, ds, abundance_dict=abun)
    field = ('gas', 'O_p5_number_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    # Test values of added field
    num_dens = ds.quan(abun["O"], "1/Zsun") \
            * ad[("gas","metallicity")] \
            * ad[("gas","O_p5_ion_fraction")] \
            * ad[("gas","H_nuclei_density")]
    assert np.allclose(num_dens, ad[field])

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

def test_add_ion_fields_to_gizmo_with_nonsolar_abundances():
    """
    Test adding a non-tracked metal to an Enzo dataset with nonsolar abundance
    """
    abun = {"Na": 2e-5}

    ds = load(FIRE_SIM)
    ad = ds.all_data()
    add_ion_number_density_field('Na', 2, ds, abundance_dict=abun)
    field = ('gas', 'Na_p1_number_density')
    assert field in ds.derived_field_list
    assert isinstance(ad[field], np.ndarray)

    # Test values of added field
    num_dens = ds.quan(abun["Na"], "1.0/Zsun") \
            * ad[("gas","metallicity")] \
            * ad[("gas","Na_p1_ion_fraction")] \
            * ad[("gas","H_nuclei_density")]
    assert np.allclose(num_dens, ad[field])

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

def test_calculate_ion_fraction():
    """
    Test to ensure calculate_ion_fraction() works
    """
    temp = np.logspace(4, 7, 100)
    dens = np.ones_like(temp) * 1e-2
    reds = np.ones_like(temp) * 0.25
    # Does it work for arrays?
    calculate_ion_fraction('H I', dens, temp, reds)

    # Does it work for single values?
    calculate_ion_fraction('Mg II', dens[0], temp[0], reds[0])

    # Does it work for other ions?
    calculate_ion_fraction('O VI', dens[0], temp[0], reds[0])

    # Does it return all hydrogen being ionized at 1e7 K?
    assert calculate_ion_fraction('H II', 1e-2, 1e7, 0) == 1

    # Can it swap ionization tables?
    ion_filepath = os.environ.get("TRIDENT_ION_DATA", 
                                  path.join(path.expanduser("~"), ".trident"))
    ionfile_hm12 = path.join(ion_filepath, "hm2012_lr.h5")
    ionfile_fg09 = path.join(ion_filepath, "fg2009_lr.h5")

    # O VI shows the biggest difference between "model families";
    # Taira et al. 2025 https://ui.adsabs.harvard.edu/abs/2025ApJ...991..221T/abstract
    frac_hm12 = calculate_ion_fraction("O VI", dens, temp, reds,
                                       ionfile_hm12)
    frac_fg09 = calculate_ion_fraction("O VI", dens, temp, reds,
                                       ionfile_fg09)

    assert not np.allclose(frac_hm12, frac_fg09)

def test_species_fraction_field_is_used_for_ion_mass_and_number_density():
    """
    Ensure Gadget/AREPO-style X_fraction fields are used in preference to the
    bulk metallicity field when X_metallicity is absent.
    """
    ds = fake_random_ds(
        8,
        fields=("density", "temperature", "metallicity", "Mg_fraction"),
        units=('g/cm**3', 'K', '', '')
    )

    ftype = ds.default_fluid_type

    def _fixed_metallicity(field, data):
        return data.ds.arr(
            np.ones(data[(ftype, "density")].shape) * 1.0e-6, ""
        )

    def _fixed_mg_fraction(field, data):
        return data.ds.arr(
            np.ones(data[(ftype, "density")].shape) * 2.0e-2, ""
        )

    ds.add_field(
        (ftype, "metallicity"),
        function=_fixed_metallicity,
        sampling_type="cell",
        units="",
        force_override=True,
    )
    ds.add_field(
        (ftype, "Mg_fraction"),
        function=_fixed_mg_fraction,
        sampling_type="cell",
        units="",
        force_override=True,
    )

    ad = ds.all_data()

    # Add the fields under test on the dataset's actual fluid type.
    add_ion_fields(ds, ["Mg II"], ftype=ftype)

    ion_frac = ad[(ftype, "Mg_p1_ion_fraction")]
    ion_mass = ad[(ftype, "Mg_p1_mass")]
    ion_ndens = ad[(ftype, "Mg_p1_number_density")]

    density = ad[(ftype, "density")]
    mass = ad[(ftype, "mass")]
    mg_fraction = ad[(ftype, "Mg_fraction")]
    bulk_metallicity = ad[(ftype, "metallicity")]

    # The fixed code should use Mg_fraction.
    expected_mass = ion_frac * mass * mg_fraction
    expected_ndens = ion_frac * density * mg_fraction / (atomic_mass["Mg"] * mh)

    assert_array_rel_equal(ion_mass, expected_mass, decimals=15)
    assert_array_rel_equal(ion_ndens, expected_ndens, decimals=15)

    # And specifically should not be using the bulk metallicity fallback.
    wrong_mass = ion_frac * mass * bulk_metallicity
    wrong_ndens = ion_frac * density * bulk_metallicity / (atomic_mass["Mg"] * mh)

    assert not np.allclose(ion_mass, wrong_mass)
    assert not np.allclose(ion_ndens, wrong_ndens)
