"""
Unit test for the AbsorptionSpectrum analysis module
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014-2017, yt Development Team.
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
from yt.convenience import load
from yt.testing import \
    assert_allclose_units, \
    assert_almost_equal

from trident.absorption_spectrum.absorption_line import \
    voigt_old, voigt_scipy
from trident.absorption_spectrum.absorption_spectrum import \
    AbsorptionSpectrum
from trident.light_ray import \
    LightRay
from trident.testing import \
    answer_test_data_dir, \
    assert_array_rel_equal, \
    h5_answer_test, \
    TempDirTest

COSMO_PLUS = os.path.join(answer_test_data_dir,
                          "enzo_cosmology_plus/AMRCosmology.enzo")
COSMO_PLUS_SINGLE = os.path.join(answer_test_data_dir,
                                 "enzo_cosmology_plus/RD0009/RD0009")
GIZMO_PLUS = os.path.join(answer_test_data_dir,
                          "gizmo_cosmology_plus/N128L16.param")
GIZMO_PLUS_SINGLE = os.path.join(answer_test_data_dir,
                                 "gizmo_cosmology_plus/snap_N128L16_151.hdf5")
ISO_GALAXY = os.path.join(answer_test_data_dir,
                          "IsolatedGalaxy/galaxy0030/galaxy0030")
FIRE = os.path.join(answer_test_data_dir,
                    "FIRE_M12i_ref11/snapshot_600.hdf5")

def test_voigt_profiles():
    a = 1.7e-4
    x = np.linspace(5.0, -3.6, 60)
    assert_allclose_units, voigt_old(a, x), voigt_scipy(a, x), 1e-8

class AbsorptionSpectrumTest(TempDirTest):

    @h5_answer_test(assert_array_rel_equal, decimals=13)
    def test_absorption_spectrum_cosmo(self):
        """
        This test generates an absorption spectrum from a compound light ray on a
        grid dataset
        """

        lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.03)

        lr.make_light_ray(seed=1234567,
                          fields=['temperature', 'density', 'H_number_density'],
                          data_filename='lightray.h5')

        sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 912.323660  # Angstroms
        normalization = 1.6e17
        index = 3.0

        sp.add_continuum(my_label, field, wavelength, normalization, index)

        filename = "spectrum.h5"
        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file=filename,
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=True)
        return filename

    @h5_answer_test(assert_array_rel_equal, decimals=15)
    def test_absorption_spectrum_non_cosmo(self):
        """
        This test generates an absorption spectrum from a simple light ray on a
        grid dataset
        """

        lr = LightRay(COSMO_PLUS_SINGLE)

        ray_start = [0,0,0]
        ray_end = [1,1,1]
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_number_density'],
                          data_filename='lightray.h5')

        sp = AbsorptionSpectrum(1200.0, 1300.0, 10001)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        filename = "spectrum.h5"
        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file=filename,
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=True)
        return filename

    @h5_answer_test(assert_array_rel_equal, decimals=15)
    def test_absorption_spectrum_non_cosmo_novpec(self):
        """
        This test generates an absorption spectrum from a simple light ray on a
        grid dataset
        """

        lr = LightRay(COSMO_PLUS_SINGLE)

        ray_start = [0,0,0]
        ray_end = [1,1,1]
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_number_density'],
                          data_filename='lightray.h5', use_peculiar_velocity=False)

        sp = AbsorptionSpectrum(1200.0, 1300.0, 10001)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        filename = "spectrum.h5"
        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file=filename,
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=False)
        return filename

    def test_equivalent_width_conserved(self):
        """
        This tests that the equivalent width of the optical depth is conserved
        regardless of the bin width employed in wavelength space.
        Unresolved lines should still deposit optical depth into the spectrum.
        """

        lr = LightRay(COSMO_PLUS_SINGLE)

        ray_start = [0,0,0]
        ray_end = [1,1,1]
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_number_density'],
                          data_filename='lightray.h5')

        my_label = 'HI Lya'
        field = 'H_number_density'
        wave = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        lambda_min= 1200
        lambda_max= 1300
        lambda_bin_widths = [1e-3, 1e-2, 1e-1, 1e0, 1e1]
        total_tau = []

        for lambda_bin_width in lambda_bin_widths:
            n_lambda = ((lambda_max - lambda_min)/ lambda_bin_width) + 1
            sp = AbsorptionSpectrum(lambda_min=lambda_min, lambda_max=lambda_max,
                                    n_lambda=n_lambda)
            sp.add_line(my_label, field, wave, f_value, gamma, mass)
            wavelength, flux = sp.make_spectrum('lightray.h5')
            total_tau.append((lambda_bin_width * sp.tau_field).sum())

        # assure that the total tau values are all within 1e-3 of each other
        for tau in total_tau:
            assert_almost_equal(tau, total_tau[0], 3)

    def test_absorption_spectrum_fits(self):
        """
        This test generates an absorption spectrum and saves it as a fits file.
        """

        lr = LightRay(COSMO_PLUS_SINGLE)

        ray_start = [0,0,0]
        ray_end = [1,1,1]
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_number_density'],
                          data_filename='lightray.h5')

        sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 912.323660  # Angstroms
        normalization = 1.6e17
        index = 3.0

        sp.add_continuum(my_label, field, wavelength, normalization, index)

        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file='spectrum.fits',
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=True)

    @h5_answer_test(assert_array_rel_equal, decimals=12)
    def test_absorption_spectrum_cosmo_sph(self):
        """
        This test generates an absorption spectrum from a compound light ray on a
        particle dataset
        """

        lr = LightRay(GIZMO_PLUS, 'Gadget', 0.0, 0.01)

        lr.make_light_ray(seed=1234567,
                          fields=[('gas', 'temperature'),
                                  ('gas', 'H_number_density')],
                          data_filename='lightray.h5')

        sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

        my_label = 'HI Lya'
        field = ('gas', 'H_number_density')
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        my_label = 'HI Lya'
        field = ('gas', 'H_number_density')
        wavelength = 912.323660  # Angstroms
        normalization = 1.6e17
        index = 3.0

        sp.add_continuum(my_label, field, wavelength, normalization, index)

        filename = "spectrum.h5"
        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file=filename,
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=True)
        return filename

    @h5_answer_test(assert_array_rel_equal, decimals=16)
    def test_absorption_spectrum_non_cosmo_sph(self):
        """
        This test generates an absorption spectrum from a simple light ray on a
        particle dataset
        """

        ds = load(GIZMO_PLUS_SINGLE)
        lr = LightRay(ds)
        ray_start = ds.domain_left_edge
        ray_end = ds.domain_right_edge
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=[('gas', 'temperature'),
                                  ('gas', 'H_number_density')],
                          data_filename='lightray.h5')

        sp = AbsorptionSpectrum(1200.0, 1300.0, 10001)

        my_label = 'HI Lya'
        field = ('gas', 'H_number_density')
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        filename = "spectrum.h5"
        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file=filename,
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=True)
        return filename

    @h5_answer_test(assert_array_rel_equal, decimals=15)
    def test_absorption_spectrum_with_continuum(self):
        """
        This test generates an absorption spectrum from a simple light ray on a
        grid dataset and adds Lyman alpha and Lyman continuum to it
        """

        ds = load(ISO_GALAXY)
        lr = LightRay(ds)

        ray_start = ds.domain_left_edge
        ray_end = ds.domain_right_edge
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_number_density'],
                          data_filename='lightray.h5')

        sp = AbsorptionSpectrum(800.0, 1300.0, 5001)

        my_label = 'HI Lya'
        field = 'H_number_density'
        wavelength = 1215.6700  # Angstromss
        f_value = 4.164E-01
        gamma = 6.265e+08
        mass = 1.00794

        sp.add_line(my_label, field, wavelength, f_value,
                    gamma, mass, label_threshold=1.e10)

        my_label = 'Ly C'
        field = 'H_number_density'
        wavelength = 912.323660  # Angstroms
        normalization = 1.6e17
        index = 3.0

        sp.add_continuum(my_label, field, wavelength, normalization, index)

        filename = "spectrum.h5"
        wavelength, flux = sp.make_spectrum('lightray.h5',
                                            output_file=filename,
                                            line_list_file='lines.txt',
                                            use_peculiar_velocity=True)
        return filename

    def test_absorption_spectrum_with_zero_field(self):
        """
        This test generates an absorption spectrum with some
        particle dataset
        """

        ds = load(FIRE)
        lr = LightRay(ds)

        # Define species and associated parameters to add to continuum
        # Parameters used for both adding the transition to the spectrum
        # and for fitting
        # Note that for single species that produce multiple lines
        # (as in the OVI doublet), 'numLines' will be equal to the number
        # of lines, and f,gamma, and wavelength will have multiple values.

        HI_parameters = {
            'name': 'HI',
            'field': 'H_number_density',
            'f': [.4164],
            'Gamma': [6.265E8],
            'wavelength': [1215.67],
            'mass': 1.00794,
            'numLines': 1,
            'maxN': 1E22, 'minN': 1E11,
            'maxb': 300, 'minb': 1,
            'maxz': 6, 'minz': 0,
            'init_b': 30,
            'init_N': 1E14
        }

        species_dicts = {'HI': HI_parameters}


        # Get all fields that need to be added to the light ray
        fields = [('gas','temperature')]
        for s, params in species_dicts.items():
            fields.append(params['field'])

        # With a single dataset, a start_position and
        # end_position or trajectory must be given.
        # Trajectory should be given as (r, theta, phi)
        lr.make_light_ray(
            start_position=ds.arr([0., 0., 0.], 'unitary'),
            end_position=ds.arr([1., 1., 1.], 'unitary'),
            solution_filename='test_lightraysolution.txt',
            data_filename='test_lightray.h5',
            fields=fields)

        # Create an AbsorptionSpectrum object extending from
        # lambda = 900 to lambda = 1800, with 10000 pixels
        sp = AbsorptionSpectrum(900.0, 1400.0, 50000)

        # Iterate over species
        for s, params in species_dicts.items():
            # Iterate over transitions for a single species
            for i in range(params['numLines']):
                # Add the lines to the spectrum
                sp.add_line(
                    s, params['field'],
                    params['wavelength'][i], params['f'][i],
                    params['Gamma'][i], params['mass'],
                    label_threshold=1.e10)


        # Make and save spectrum
        wavelength, flux = sp.make_spectrum(
            'test_lightray.h5',
            output_file='test_spectrum.h5',
            line_list_file='test_lines.txt',
            use_peculiar_velocity=True)
