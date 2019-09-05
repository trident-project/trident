"""
tests for auto-lambda feature
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
    TempDirTest

COSMO_PLUS_SINGLE = os.path.join(answer_test_data_dir,
                                 "enzo_cosmology_plus/RD0009/RD0009")

def compare_spectra(sg1, sg2, comp_key):
    assert_allclose(
        sg1.tau_field, sg2.tau_field, rtol=1e-6,
        err_msg='tau arrays for auto and %s disagree!' % comp_key)
    assert_allclose(
        sg1.lambda_field, sg2.lambda_field, rtol=1e-10,
        err_msg='lambda arrays for auto and %s disagree!' % comp_key)

    for my_line in sg1.line_observables_dict:
        for key in sg1.line_observables_dict[my_line]:
            assert_allclose(
                sg1.line_observables_dict[my_line][key],
                sg2.line_observables_dict[my_line][key], rtol=1e-7,
                err_msg='%s field for auto and %s disagree!' % (key, comp_key))

class AutoLambdaTest(TempDirTest):

    def setUp(self):
        super(AutoLambdaTest, self).setUp()
        ds = load(COSMO_PLUS_SINGLE)
        line_list = ['H I 1216', 'H I 1026']

        make_simple_ray(ds, start_position=ds.domain_left_edge,
                        end_position=ds.domain_right_edge,
                        data_filename='ray.h5',
                        lines=line_list, ftype='gas')
        self.line_list = line_list

    def test_dlambda(self):
        """
        Compare auto wavelength against setting same min, max, and dlambda
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=0.01)
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        sg_comp = SpectrumGenerator(
            lambda_min=sg_auto.lambda_field[0].d,
            lambda_max=sg_auto.lambda_field[-1].d,
            dlambda=sg_auto.bin_width.d)
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)
        compare_spectra(sg_auto, sg_comp, 'manually setting dlambda')

    def test_n_lambda(self):
        """
        Compare auto wavelength against setting same min, max, and n_lambda
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=0.01)
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        sg_comp = SpectrumGenerator(
            lambda_min=sg_auto.lambda_field[0].d,
            lambda_max=sg_auto.lambda_field[-1].d,
            n_lambda=sg_auto.lambda_field.size)
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)
        compare_spectra(sg_auto, sg_comp, 'manually setting n_lambda')

    def test_dlambda_extra_wide(self):
        """
        Compare auto wavelength against setting same dlambda with
        min and max set extra wide.
        """

        sg_auto = SpectrumGenerator(
            lambda_min='auto', lambda_max='auto',
            dlambda=0.01)
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              store_observables=True)

        sg_comp = SpectrumGenerator(
            lambda_min=sg_auto.lambda_field[0].d,
            lambda_max=sg_auto.lambda_field[-1].d,
            n_lambda=sg_auto.lambda_field.size)
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

        sg = SpectrumGenerator(lambda_max=1100,
                               lambda_min='auto',
                               dlambda=0.01)
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
            dlambda=0.01)
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        sg_comp = SpectrumGenerator(
            lambda_min=1100, lambda_max='auto',
            dlambda=0.01)
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
            dlambda=0.01)
        sg_auto.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        sg_comp = SpectrumGenerator(
            lambda_min='auto', lambda_max=1100,
            dlambda=0.01)
        sg_comp.make_spectrum("ray.h5", lines=self.line_list,
                              ly_continuum=False)

        comp_lambda = sg_auto.lambda_field <= sg_comp.lambda_field[-1]

        assert_allclose(
            sg_auto.tau_field[comp_lambda].sum(),
            sg_comp.tau_field.sum())
