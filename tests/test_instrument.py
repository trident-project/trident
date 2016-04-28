"""
Tests for Instrument Class and Functions

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from trident.instrument import Instrument

def test_instrument():
    """
    Create a couple of Instrument objects and assure that their values
    are stored appropriately and that the __repr__ class works OK.
    """
    inst = Instrument(1000, 2000, n_lambda=101)
    assert inst.lambda_min == 1000
    assert inst.lambda_max == 2000
    assert inst.n_lambda == 101
    assert inst.dlambda == 10
    print inst

    inst = Instrument(100, 20000, dlambda=10, lsf_kernel='COS', name='test')
    assert inst.lambda_min == 100
    assert inst.lambda_max == 20000
    assert inst.n_lambda == 1991
    assert inst.dlambda == 10
    assert inst.lsf_kernel == 'COS'
    assert inst.name == 'test'
    print inst
