"""
tests for velocity bin space
"""

#-----------------------------------------------------------------------------
# Copyright (c) Trident Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
from yt.convenience import load
from yt.testing import \
    assert_allclose

from trident import \
    make_simple_ray, \
    SpectrumGenerator
from trident.testing import \
    answer_test_data_dir, \
    compare_spectra, \
    TempDirTest

COSMO_PLUS_SINGLE = os.path.join(answer_test_data_dir,
                                 "enzo_cosmology_plus/RD0009/RD0009")

class VelocitySpaceTest(TempDirTest):

    def setUp(self):
        super(VelocitySpaceTest, self).setUp()
        ds = load(COSMO_PLUS_SINGLE)
        line_list = ['H I 1216', 'H I 1026']

        make_simple_ray(ds, start_position=ds.domain_left_edge,
                        end_position=ds.domain_right_edge,
                        data_filename='ray.h5',
                        lines=line_list, ftype='gas')
        self.line_list = line_list

    def test_dlambda(self):
        """
        Compare auto velocity against setting same min, max, and dlambda
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=1.0, bin_space='velocity')
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        sg_comp = SpectrumGenerator(
            lambda_min=sg_auto.lambda_field[0],
            lambda_max=sg_auto.lambda_field[-1],
            dlambda=sg_auto.bin_width,
            bin_space='velocity')
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)
        compare_spectra(sg_auto, sg_comp, 'manually setting dlambda')

    def test_n_lambda(self):
        """
        Compare auto velocity against setting same min, max, and n_lambda
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=1.0, bin_space='velocity')
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        sg_comp = SpectrumGenerator(
            lambda_min=sg_auto.lambda_field[0],
            lambda_max=sg_auto.lambda_field[-1],
            n_lambda=sg_auto.lambda_field.size,
            bin_space='velocity')
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)
        compare_spectra(sg_auto, sg_comp, 'manually setting n_lambda')

    def test_dlambda_extra_wide(self):
        """
        Compare auto velocity against setting same dlambda with
        min and max set extra wide.
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=1.0, bin_space='velocity')
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        sg_comp = SpectrumGenerator(
            lambda_min=sg_auto.lambda_field[0],
            lambda_max=sg_auto.lambda_field[-1],
            n_lambda=sg_auto.lambda_field.size,
            bin_space='velocity')
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        assert_allclose(
            sg_auto.tau_field.sum(),
            sg_comp.tau_field.sum(), rtol=1e-7,
            err_msg='Total tau disagrees with extra wide spectrum.')

    def test_empty_spectrum(self):
        """
        Test we can make an empty spectrum without crashing.
        """

        sg = SpectrumGenerator(lambda_min=3000,
                               lambda_max='auto',
                               dlambda=1.0,
                               bin_space='velocity')
        sg.make_spectrum("ray.h5", lines=['H I 1216'],
                         output_file='blank.h5',
                         output_absorbers_file='blank.txt',
                         store_observables=True, ly_continuum=True)
        sg.plot_spectrum('spec_blank.png')
        assert sg.lambda_field is None
        assert sg.tau_field is None
        assert sg.flux_field is None

    def test_setting_lambda_min(self):
        """
        Test setting an arbitrary lambda_min with auto-lambda.
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=1.0, bin_space='velocity')
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        sg_comp = SpectrumGenerator(
            lambda_min=-10000., lambda_max='auto',
            dlambda=1.0, bin_space='velocity')
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        comp_lambda = sg_auto.lambda_field >= sg_comp.lambda_field[0]

        assert_allclose(
            sg_auto.tau_field[comp_lambda].sum(),
            sg_comp.tau_field.sum())

    def test_setting_lambda_max(self):
        """
        Test setting an arbitrary lambda_max with auto-lambda.
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=1.0, bin_space='velocity')
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        sg_comp = SpectrumGenerator(
            lambda_min='auto', lambda_max=-10000.,
            dlambda=1.0, bin_space='velocity')
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        comp_lambda = sg_auto.lambda_field <= sg_comp.lambda_field[-1]

        assert_allclose(
            sg_auto.tau_field[comp_lambda].sum(),
            sg_comp.tau_field.sum())
