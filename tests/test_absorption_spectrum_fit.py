"""
Unit test for the AbsorptionSpectrum analysis module
"""

#-----------------------------------------------------------------------------
# Copyright (c) Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
from yt.convenience import load

from trident import \
    make_simple_ray, \
    SpectrumGenerator
    
from trident.absorption_spectrum.absorption_spectrum_fit import \
    generate_total_fit
from trident.testing import \
    answer_test_data_dir, \
    h5_answer_test, \
    TempDirTest

ISO_GALAXY = os.path.join(answer_test_data_dir,
                          "IsolatedGalaxy/galaxy0030/galaxy0030")

class AbsorptionSpectrumFitTest(TempDirTest):

    def test_absorption_spectrum_with_continuum(self):
        """
        This test generates an absorption spectrum then does a default fit.
        Uses the example in the docs as a template:
        https://trident.readthedocs.io/en/latest/absorption_spectrum_fit.html
        """

        ds = load(ISO_GALAXY)
        line_list = ['O VI']
        ray = make_simple_ray(ds, start_position=ds.domain_left_edge,
                                      end_position=ds.domain_right_edge,
                                      data_filename='ray.h5',
                                      lines=line_list, ftype='gas')
        sg = SpectrumGenerator('COS')
        sg.make_spectrum(ray, lines=line_list)
        wavelength = sg.lambda_field
        flux = sg.flux_field

        OVI_parameters = {'name':'OVI',
                          'f':[.1325,.06580],
                          'Gamma':[4.148E8,4.076E8],
                          'wavelength':[1031.9261,1037.6167],
                          'numLines':2,
                          'maxN':1E17,'minN':1E11,
                          'maxb':300, 'minb':1,
                          'maxz':6, 'minz':0,
                          'init_b':20,
                          'init_N':1E12}

        speciesDicts = {'OVI':OVI_parameters}

        orderFits = ['OVI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength,
                                                       flux, 
                                                       orderFits, 
                                                       speciesDicts)
