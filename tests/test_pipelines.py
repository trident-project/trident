"""
Tests to assure the code runs on example pipelines

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import os
from yt.convenience import load
from yt.testing import \
    assert_almost_equal

from trident import \
    verify, \
    make_simple_ray, \
    make_compound_ray, \
    SpectrumGenerator
from trident.testing import \
    answer_test_data_dir, \
    generate_results, \
    test_results_dir, \
    TempDirTest

# Datasets required for running some answer tests
enzo_small = os.path.join(answer_test_data_dir, 'enzo_small')
err_precision = 8

def test_verify():
    """
    Test that the verify function works and the code can create a dataset,
    generate a LightRay and then construct a spectrum from it.
    """
    verify()
    assert True

class PipelineTest(TempDirTest):

    def test_enzo_small_simple(self):
        """
        This is an answer test, which compares the results of this test
        against answers generated from a previous version of the code.

        This test generates a COS spectrum from a single Enzo dataset 
        using a simple ray and compare the ray and spectral output data
        against a known answer.

        """

        # Set the dataset filename, load it into yt and define the trajectory
        # of the LightRay to cross the box from one corner to the other.
        ds = load(os.path.join(enzo_small, 'RD0009/RD0009'))
        ray_start = ds.domain_left_edge
        ray_end = ds.domain_right_edge

        # Make a LightRay object including all necessary fields so you can add
        # all H, C, N, O, and Mg fields to the resulting spectrum from your dataset.
        # Save LightRay to ray.h5 and use it locally as ray object.
        ray_fn = 'enzo_small_simple_ray.h5'
        ray = make_simple_ray(ds, start_position=ray_start,
                              end_position=ray_end, data_filename=ray_fn,
                                      lines=['H', 'C', 'N', 'O', 'Mg'], ftype='gas')

        # Now use the ray object to actually generate an absorption spectrum
        # Use the settings (spectral range, LSF, and spectral resolution) for COS
        # And save it as an output hdf5 file and plot it to an image.
        sg = SpectrumGenerator('COS')
        sg.make_spectrum(ray, lines=['H', 'C', 'N', 'O', 'Mg'])
        raw_file = 'enzo_small_simple_spec_raw.h5'
        raw_file_compare = os.path.join(test_results_dir, raw_file)
        sg.save_spectrum(raw_file)
        sg.plot_spectrum('enzo_small_simple_spec_raw.png')

        # "Final" spectrum with added quasar, MW background, applied line-spread
        # function, and added gaussian noise (SNR=30)
        # Save as a text file and plot it to an image.
        sg.add_qso_spectrum()
        sg.add_milky_way_foreground()
        sg.apply_lsf()
        sg.add_gaussian_noise(30, seed=1)
        final_file = 'enzo_small_simple_spec_final.h5'
        final_file_compare = os.path.join(test_results_dir, final_file)
        sg.save_spectrum(final_file)
        sg.plot_spectrum('enzo_small_simple_spec_final.png')

        if generate_results:
            os.rename(raw_file, raw_file_compare)
            os.rename(final_file, final_file_compare)

        else:
            old_spec = h5py.File(raw_file_compare, 'r')
            new_spec = h5py.File(raw_file, 'r')
            for key in old_spec.keys():
                assert_almost_equal(new_spec[key].value, old_spec[key].value, \
                                    decimal=err_precision,
                                    err_msg='Raw spectrum array does not match '+\
                                    'for enzo_small_simple answer test')
            old_spec.close()
            new_spec.close()

            old_spec = h5py.File(final_file_compare, 'r')
            new_spec = h5py.File(final_file, 'r')
            for key in old_spec.keys():
                assert_almost_equal(new_spec[key].value, old_spec[key].value, \
                                    decimal=err_precision,
                                    err_msg='Final spectrum array does not match '+\
                                    'for enzo_small_simple answer test')
            old_spec.close()
            new_spec.close()

    def test_enzo_small_compound(self):
        """
        This is an answer test, which compares the results of this test
        against answers generated from a previous version of the code.

        This test generates a COS spectrum from a single Enzo dataset 
        using a compound ray and compare the ray and spectral output data
        against a known answer.
        """

        # Set the dataset filename, the redshift range to be covered, and the
        # ions in which we're interested, and give the random number generator
        # a seed (for determining the ray trajectory).  Use these to create a
        # compound ray spanning multiple datasets to cover the desired redshift
        # range.
        fn = os.path.join(enzo_small, 'AMRCosmology.enzo')
        redshift_start = 0.00
        redshift_end = 0.015
        lines = ['H', 'C', 'N', 'O', 'Mg']
        seed = 1
        ray_fn = 'enzo_small_compound_ray.h5'
        ray = make_compound_ray(fn, 'Enzo', redshift_start, redshift_end,
                                lines=lines, seed=seed, ftype='gas',
                                data_filename=ray_fn)

        # Now use the ray object to actually generate an absorption spectrum
        # Use the settings (spectral range, LSF, and spectral resolution) for COS
        # And save it as an output hdf5 file and plot it to an image.
        sg = SpectrumGenerator('COS')
        sg.make_spectrum(ray, lines=lines)
        raw_file = 'enzo_small_compound_spec_raw.h5'
        raw_file_compare = os.path.join(test_results_dir, raw_file)
        sg.save_spectrum(raw_file)
        sg.plot_spectrum('enzo_small_compound_spec_raw.png')

        # "Final" spectrum with added quasar, MW background, applied line-spread
        # function, and added gaussian noise (SNR=30)
        # Save as a text file and plot it to an image.
        sg.add_qso_spectrum()
        sg.add_milky_way_foreground()
        sg.apply_lsf()
        sg.add_gaussian_noise(30, seed=1)
        final_file = 'enzo_small_compound_spec_final.h5'
        sg.save_spectrum(final_file)
        final_file_compare = os.path.join(final_file)
        sg.plot_spectrum('enzo_small_compound_spec_final.png')

        if generate_results:
            os.rename(raw_file, raw_file_compare)
            os.rename(final_file, final_file_compare)

        else:
            old_spec = h5py.File(raw_file_compare, 'r')
            new_spec = h5py.File(raw_file, 'r')
            for key in old_spec.keys():
                assert_almost_equal(new_spec[key].value, old_spec[key].value, \
                                    decimal=err_precision,
                                    err_msg='Raw spectrum array does not match '+\
                                    'for enzo_small_compound answer test')
            old_spec.close()
            new_spec.close()

            old_spec = h5py.File(final_file_compare, 'r')
            new_spec = h5py.File(final_file, 'r')
            for key in old_spec.keys():
                assert_almost_equal(new_spec[key].value, old_spec[key].value, \
                                    decimal=err_precision,
                                    err_msg='Final spectrum array does not match '+\
                                    'for enzo_small_compound answer test')
            old_spec.close()
            new_spec.close()
