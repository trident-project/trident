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
from numpy.testing import \
    assert_allclose, \
    assert_array_equal
import os
import shutil
import tempfile
from unittest import \
    TestCase
from yt.funcs import \
    ensure_dir
from yt.testing import \
    assert_rel_equal
from trident.config import \
    parse_config

# If TRIDENT_GENERATE_TEST_RESULTS=1, just generate test results.
generate_results = int(os.environ.get("TRIDENT_GENERATE_TEST_RESULTS", 0)) == 1
answer_test_data_dir = ensure_dir(
  os.path.abspath(os.path.expanduser(parse_config('answer_test_data_dir'))))
test_results_dir = os.path.join(answer_test_data_dir, "test_results")
if generate_results:
    ensure_dir(test_results_dir)

class TempDirTest(TestCase):
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

def h5_answer_test(compare, **kwargs):
    """
    HDF5 answer test decorator.

    Put this decorator above testing functions that return a
    filename for a file generated within that file.

    If the environment variable, TRIDENT_GENERATE_TEST_RESULTS is 1,
    then the filename will be renamed to the name of the test
    function and stored.

    If the TRIDENT_GENERATE_TEST_RESULTS is 0, then check to make sure
    the comparison file exists, then compare this file with the
    output of the test function.
    """

    def real_h5_answer_test(func):
        def wrapper(*args):
            # name the file after the function
            filename = "%s.h5" % func.__name__
            result_filename = os.path.join(test_results_dir, filename)

            # check that answers exist
            if not generate_results:
                assert os.path.exists(result_filename), \
                  "Result file, %s, not found!" % result_filename

            output_filename = func(*args)
            # if generating, move files to results dir
            if generate_results:
                os.rename(output_filename, result_filename)
            # if comparing, run the comparison
            else:
                h5_dataset_compare(output_filename, result_filename,
                                   compare=compare, **kwargs)
        return wrapper
    return real_h5_answer_test

def h5_dataset_compare(fn1, fn2, compare=None, **kwargs):
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
        compare(fh1[key][()], fh2[key][()],
                err_msg="%s field not equal!" % key, **kwargs)

def assert_array_rel_equal(a1, a2, decimals=16, **kwargs):
    """
    Wraps assert_rel_equal with, but decimals is a keyword arg.
    """
    assert_rel_equal(a1, a2, decimals, **kwargs)


def compare_spectra(sg1, sg2, comp_key):
    """
    Compare two spectra.
    """
    assert_allclose(
        sg1.tau_field, sg2.tau_field, rtol=1e-7,
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
