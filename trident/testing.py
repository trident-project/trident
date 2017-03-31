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

import h5py
from functools import wraps
from numpy.testing import \
    assert_array_equal
import os
import shutil
import tempfile
from unittest import \
    TestCase
from yt.funcs import \
    ensure_dir
from trident.utilities import \
    parse_config

# If GENERATE_TEST_RESULTS=1, just generate test results.
generate_results = int(os.environ.get("GENERATE_TEST_RESULTS", 0)) == 1
answer_test_data_dir = \
  os.path.abspath(os.path.expanduser(parse_config('answer_test_data_dir')))
test_results_dir = ensure_dir(
  os.path.join(answer_test_data_dir, "test_results"))

class TestInTempDir(TestCase):
    """
    A test class that runs in a temporary directory and removes it afterward.
    """

    def setUp(self):
        self.curdir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

def h5_answer_test(func):
    """
    HDF5 answer test decorator.

    Put this decorator above testing functions that return a
    filename for a file generated within that file.

    If the environment variable, GENERATE_TEST_RESULTS is 1,
    then the filename will be renamed to the name of the test
    function and stored.

    If the GENERATE_TEST_RESULTS is 0, then check to make sure
    the comparison file exists, then compare this file with the
    output of the test function.
    """

    def do_h5_answer_test(*args, **kwargs):
        # name the file after the function
        filename = "%s.h5" % func.__name__
        result_filename = os.path.join(test_results_dir, filename)

        if not generate_results:
            assert os.path.exists(result_filename), \
              "Result file, %s, not found!" % result_filename

        output_filename = func(*args, **kwargs)
        if generate_results:
            os.rename(output_filename, result_filename)
        else:
            h5_dataset_compare(output_filename, result_filename)

    return do_h5_answer_test

def h5_dataset_compare(fn1, fn2, compare=None):
    """
    Compare all datasets between two hdf5 files.
    """

    if compare is None:
        compare = assert_array_equal
    fh1 = h5py.File(fn1, "r")
    fh2 = h5py.File(fn2, "r")
    assert list(fh1.keys()) == list(fh2.keys()), \
      "Files have different datasets!"
    for key in fh1.keys():
        compare(fh1[key].value, fh2[key].value,
                err_msg="%s field not equal!" % key)
