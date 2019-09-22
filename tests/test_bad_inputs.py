"""
tests for giving bad inputs
"""

#-----------------------------------------------------------------------------
# Copyright (c) Trident Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from numpy.testing import \
    assert_raises

from trident.absorption_spectrum.absorption_spectrum import \
    AbsorptionSpectrum

def test_no_n_lambda_or_dlambda():
    with assert_raises(RuntimeError):
        AbsorptionSpectrum(
            lambda_min=1100, lambda_max=1200)

def test_n_lambda_with_auto():
    with assert_raises(RuntimeError):
        AbsorptionSpectrum(
            lambda_min='auto', lambda_max='auto',
            n_lambda=1000)

def test_min_greater_than_max():
    with assert_raises(RuntimeError):
        AbsorptionSpectrum(
            lambda_min=1200, lambda_max=1100,
            n_lambda=1000)
