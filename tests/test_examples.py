"""
Tests to assure the code runs on examples

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from trident import verify

def test_verify():
    """
    Test that the verify function works and the code can create a dataset,
    generate a LightRay and then construct a spectrum from it.
    """
    verify()
    assert True