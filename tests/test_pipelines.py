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

from __future__ import absolute_import
from yt import load
from trident import \
    verify, \
    make_simple_ray, \
    make_compound_ray, \
    SpectrumGenerator, \
    parse_config
from yt.testing import \
    assert_equal, \
    requires_file
import tempfile
import shutil
import os
import h5py as h5

# Datasets required for running some answer tests
answer_test_data_dir = parse_config('answer_test_data_dir')
enzo_small = os.path.join(answer_test_data_dir, 'enzo_small')

def test_verify():
    """
    Test that the verify function works and the code can create a dataset,
    generate a LightRay and then construct a spectrum from it.
    """
    verify()
    assert True

@requires_file(enzo_small)
def test_enzo_small_simple(generate_answers=False):
    """
    This is an answer test, which compares the results of this test
    against answers generated from a previous version of the code.

    ANSWER TESTING REQUIRES SMALL SAMPLE DATASET AND ANSWER DATASET

    First, get the dataset and answer for this dataset.
    Then, unzip and untar the dataset somewhere in your filesystem.
    Finally, modify your .trident/config.tri file so trident knows where you
    put the unzipped dataset.

    All of this can be accomplished with the commands:

    $ wget http://trident-project.org/data/tests/enzo_small.tar.gz
    $ tar -zxvf enzo_small.tar.gz
    $ mv enzo_small <wherever/you/want>
    $ echo answer_test_data_dir = <wherever/you/want> >> ~/.trident/config.tri

    Finally, this test generates a COS spectrum from a single Enzo dataset 
    using a simple ray and compare the ray and spectral output data
    against a known answer.

    If you want to generate the answer tests because a substantial change
    has occurred in the code necessitating it.  Import this function 
    explicitly in a python shell and run it with generate_answers=True.  Copy
    the h5 files to the same directory as the enzo_small directory.  You 
    may want to contact the developers about making this replace the 
    standard test files available for download on http://trident-project.org
    """
    # Figure out where sample datafile is placed
    answer_test_data_dir = parse_config('answer_test_data_dir')

    # If not generating new gold standard answers, 
    # Set up temp directory to store output data so we can delete everything 
    # after having saved it
    if not generate_answers:
        tempdir = tempfile.mkdtemp()
    else:
        tempdir = '.'

    # Set the dataset filename, load it into yt and define the trajectory
    # of the LightRay to cross the box from one corner to the other.
    ds = load(os.path.join(enzo_small, 'RD0009/RD0009'))
    ray_start = ds.domain_left_edge
    ray_end = ds.domain_right_edge
    ray = os.path.join(tempdir, 'enzo_small_simple_ray.h5')

    # Make a LightRay object including all necessary fields so you can add
    # all H, C, N, O, and Mg fields to the resulting spectrum from your dataset.
    # Save LightRay to ray.h5 and use it locally as ray object.
    ray = make_simple_ray(ds, start_position=ray_start,
                          end_position=ray_end, data_filename=ray,
                                  lines=['H', 'C', 'N', 'O', 'Mg'], ftype='gas')

    # Now use the ray object to actually generate an absorption spectrum
    # Use the settings (spectral range, LSF, and spectral resolution) for COS
    # And save it as an output hdf5 file and plot it to an image.
    sg = SpectrumGenerator('COS')
    sg.make_spectrum(ray, lines=['H', 'C', 'N', 'O', 'Mg'])
    sg.save_spectrum(os.path.join(tempdir, 'enzo_small_simple_spec_raw.h5'))
    sg.plot_spectrum(os.path.join(tempdir, 'enzo_small_simple_spec_raw.png'))

    # "Final" spectrum with added quasar, MW background, applied line-spread
    # function, and added gaussian noise (SNR=30)
    # Save as a text file and plot it to an image.
    sg.add_qso_spectrum()
    sg.add_milky_way_foreground()
    sg.apply_lsf()
    sg.add_gaussian_noise(30, seed=1)
    sg.save_spectrum(os.path.join(tempdir, 'enzo_small_simple_spec_final.h5'))
    sg.plot_spectrum(os.path.join(tempdir, 'enzo_small_simple_spec_final.png'))

    # If not generating answer test files, 
    # Compare output spectra against gold standard files
    if not generate_answers:
        old_spec_fn = os.path.join(enzo_small, \
                                   'enzo_small_simple_spec_raw.h5')
        old_spec = h5.File(old_spec_fn, 'r')
        new_spec_fn = os.path.join(tempdir, 'enzo_small_simple_spec_raw.h5')
        new_spec = h5.File(new_spec_fn, 'r')
        for key in old_spec.keys():
            assert_equal(new_spec[key].value, old_spec[key].value, \
                        'Raw spectrum array does not match for '+\
                        'enzo_small_simple answer test')
        old_spec.close()
        new_spec.close()

        old_spec_fn = os.path.join(enzo_small, \
                                   'enzo_small_simple_spec_final.h5')
        old_spec = h5.File(old_spec_fn, 'r')
        new_spec_fn = os.path.join(tempdir, 'enzo_small_simple_spec_final.h5')
        new_spec = h5.File(new_spec_fn, 'r')
        for key in old_spec.keys():
            assert_equal(new_spec[key].value, old_spec[key].value, \
                        'Final spectrum array does not match for '+\
                        'enzo_small_simple answer test')
        old_spec.close()
        new_spec.close()

    # Cleaning up if not generating the answer gold standards
    if not generate_answers:
        shutil.rmtree(tempdir)

@requires_file(enzo_small)
def test_enzo_small_compound(generate_answers=False):
    """
    This is an answer test, which compares the results of this test
    against answers generated from a previous version of the code.

    ANSWER TESTING REQUIRES SMALL SAMPLE DATASET AND ANSWER DATASET

    First, get the dataset and answer for this dataset.
    Then, unzip and untar the dataset somewhere in your filesystem.
    Finally, modify your .trident/config.tri file so trident knows where you
    put the unzipped dataset.

    All of this can be accomplished with the commands:

    $ wget http://trident-project.org/data/tests/enzo_small.tar.gz
    $ tar -zxvf enzo_small.tar.gz
    $ mv enzo_small <wherever/you/want>
    $ echo answer_test_data_dir = <wherever/you/want> >> ~/.trident/config.tri

    Finally, this test generates a COS spectrum from a single Enzo dataset 
    using a compound ray and compare the ray and spectral output data
    against a known answer.

    If you want to generate the answer tests because a substantial change
    has occurred in the code necessitating it.  Import this function 
    explicitly in a python shell and run it with generate_answers=True.  Copy
    the h5 files to the same directory as the enzo_small directory.  You 
    may want to contact the developers about making this replace the 
    standard test files available for download on http://trident-project.org
    """
    # Figure out where sample datafile is placed
    answer_test_data_dir = parse_config('answer_test_data_dir')

    # If not generating new gold standard answers, 
    # Set up temp directory to store output data so we can delete everything 
    # after having saved it
    if not generate_answers:
        tempdir = tempfile.mkdtemp()
    else:
        tempdir = '.'

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
    ray_fn = os.path.join(tempdir, 'enzo_small_compound_ray.h5')
    ray = make_compound_ray(fn, 'Enzo', redshift_start, redshift_end,
                            lines=lines, seed=seed, ftype='gas',
                            data_filename=ray_fn)

    # Now use the ray object to actually generate an absorption spectrum
    # Use the settings (spectral range, LSF, and spectral resolution) for COS
    # And save it as an output hdf5 file and plot it to an image.
    sg = SpectrumGenerator('COS')
    sg.make_spectrum(ray, lines=lines)
    sg.save_spectrum(os.path.join(tempdir, 'enzo_small_compound_spec_raw.h5'))
    sg.plot_spectrum(os.path.join(tempdir, 'enzo_small_compound_spec_raw.png'))
    
    # "Final" spectrum with added quasar, MW background, applied line-spread
    # function, and added gaussian noise (SNR=30)
    # Save as a text file and plot it to an image.
    sg.add_qso_spectrum()
    sg.add_milky_way_foreground()
    sg.apply_lsf()
    sg.add_gaussian_noise(30, seed=1)
    sg.save_spectrum(os.path.join(tempdir, 'enzo_small_compound_spec_final.h5'))
    sg.plot_spectrum(os.path.join(tempdir, 'enzo_small_compound_spec_final.png'))
    
    # If not generating answer test files, 
    # Compare output spectra against gold standard files
    if not generate_answers:
        old_spec_fn = os.path.join(enzo_small, \
                                   'enzo_small_compound_spec_raw.h5')
        old_spec = h5.File(old_spec_fn, 'r')
        new_spec_fn = os.path.join(tempdir, 'enzo_small_compound_spec_raw.h5')
        new_spec = h5.File(new_spec_fn, 'r')
        for key in old_spec.keys():
            assert_equal(new_spec[key].value, old_spec[key].value, \
                        'Raw spectrum array does not match for '+\
                        'enzo_small_compound answer test')
        old_spec.close()
        new_spec.close()

        old_spec_fn = os.path.join(enzo_small, \
                                   'enzo_small_compound_spec_final.h5')
        old_spec = h5.File(old_spec_fn, 'r')
        new_spec_fn = os.path.join(tempdir, 'enzo_small_compound_spec_final.h5')
        new_spec = h5.File(new_spec_fn, 'r')
        for key in old_spec.keys():
            assert_equal(new_spec[key].value, old_spec[key].value, \
                        'Final spectrum array does not match for '+\
                        'enzo_small_compound answer test')
        old_spec.close()
        new_spec.close()

    # Cleaning up if not generating the answer gold standards
    if not generate_answers:
        shutil.rmtree(tempdir)
