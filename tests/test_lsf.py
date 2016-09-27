"""
Tests for Line Spread Function class and functions

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import
from __future__ import print_function
import pytest
from trident.lsf import LSF
import numpy as np
import astropy

def test_lsf_errors():
    """
    Test that LSF raises errors correctly
    """
    with pytest.raises(RuntimeError):
        lsf = LSF()
    with pytest.raises(RuntimeError):
        lsf = LSF(filename='NOTAFILE.txt')

def test_lsf_from_filename():
    """
    Create an LSF object from filename and assure that it works OK
    """
    lsf = LSF(filename='avg_COS.txt')
    assert isinstance(lsf.kernel, np.ndarray)
    assert len(lsf.kernel) > 0
    assert lsf.width == len(lsf.kernel)
    print(lsf)

def test_boxcar_lsf():
    """
    Create an LSF object with a boxcar function and assure it works OK
    """
    lsf = LSF(function='boxcar', width=5)
    assert isinstance(lsf.kernel, np.ndarray)
    assert len(lsf.kernel) == 5
    assert lsf.width == 5
    print(lsf)

def test_gaussian_lsf():
    """
    Create an LSF object with a gaussian function and assure it works OK
    """
    lsf = LSF(function='gaussian', width=5)
    assert isinstance(lsf.kernel, astropy.convolution.kernels.Gaussian1DKernel)
    assert len(lsf.kernel.array) == 41
    assert lsf.width == 5
    print(lsf)

