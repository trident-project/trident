import h5py
from functools import wraps
from numpy.testing import \
    assert_array_equal
import os
import shutil
import tempfile

"""
Testing utilities for Trident

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

def in_tmpdir(func, *args, **kwargs):
    """
    Make a temp dir, cd into it, run operation,
    return to original location, remove temp dir.
    """

    @wraps(func)
    def do_in_tmpdir(*args, **kwargs):
        tmpdir = tempfile.mkdtemp()
        curdir = os.getcwd()
        os.chdir(tmpdir)
        func(*args, **kwargs)
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return do_in_tmpdir

def h5_dataset_compare(fn1, fn2):
    """
    Compare all datasets between two hdf5 files.
    """

    fh1 = h5py.File(fn1, "r")
    fh2 = h5py.File(fn2, "r")
    assert list(fh1.keys()) == list(fh2.keys()), \
      "Files have different datasets!"
    for key in fh1.keys():
        assert_array_equal(fh1[key].value, fh2[key].value,
                           err_msg="%s field not equal!" % key)
