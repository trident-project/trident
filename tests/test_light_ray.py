"""
Unit test for the light_ray analysis module
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016-2017, yt Development Team.
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.convenience import \
    load
from yt.testing import \
    assert_array_equal, \
    assert_almost_equal
from trident import \
    LightRay, \
    make_simple_ray
from trident.testing import \
    answer_test_data_dir, \
    TempDirTest
import os

COSMO_PLUS = os.path.join(answer_test_data_dir,
                          "enzo_cosmology_plus/AMRCosmology.enzo")
COSMO_PLUS_SINGLE = os.path.join(answer_test_data_dir,
                                 "enzo_cosmology_plus/RD0009/RD0009")
GIZMO_COSMO_SINGLE = os.path.join(answer_test_data_dir,
                                 "gizmo_cosmology_plus/snap_N128L16_150.hdf5")

def compare_light_ray_solutions(lr1, lr2):
    assert len(lr1.light_ray_solution) == len(lr2.light_ray_solution)
    if len(lr1.light_ray_solution) == 0:
        return
    for s1, s2 in zip(lr1.light_ray_solution, lr2.light_ray_solution):
        for field in s1:
            if field in ["next", "previous"]:
                continue
            if isinstance(s1[field], np.ndarray):
                assert_array_equal(s1[field], s2[field])
            else:
                assert s1[field] == s2[field]

class LightRayTest(TempDirTest):

    def test_light_ray_cosmo(self):
        """
        This test generates a cosmological light ray
        """
        lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.03)

        lr.make_light_ray(seed=1234567,
                          fields=['temperature', 'density', 'H_p0_number_density'],
                          data_filename='lightray.h5')

        ds = load('lightray.h5')
        compare_light_ray_solutions(lr, ds)

    def test_light_ray_cosmo_nested(self):
        """
        This test generates a cosmological light ray confing the ray to a subvolume
        """
        left = np.ones(3) * 0.25
        right = np.ones(3) * 0.75

        lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.03)

        lr.make_light_ray(seed=1234567, left_edge=left, right_edge=right,
                          fields=['temperature', 'density', 'H_p0_number_density'],
                          data_filename='lightray.h5')

        ds = load('lightray.h5')
        compare_light_ray_solutions(lr, ds)

    def test_light_ray_cosmo_nonperiodic(self):
        """
        This test generates a cosmological light ray using non-periodic segments
        """
        lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.03)

        lr.make_light_ray(seed=1234567, periodic=False,
                          fields=['temperature', 'density', 'H_p0_number_density'],
                          data_filename='lightray.h5')

        ds = load('lightray.h5')
        compare_light_ray_solutions(lr, ds)

    def test_light_ray_non_cosmo(self):
        """
        This test generates a non-cosmological light ray
        """
        lr = LightRay(COSMO_PLUS_SINGLE)

        ray_start = [0,0,0]
        ray_end = [1,1,1]
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_p0_number_density'],
                          data_filename='lightray.h5')

        ds = load('lightray.h5')
        compare_light_ray_solutions(lr, ds)

    def test_light_ray_redshift_coverage_grid(self):
        """
        Tests to assure a light ray covers the full redshift range appropriate
        for that comoving line of sight distance.  Was not always so!
        """
        ds = load(COSMO_PLUS_SINGLE)
        ray = make_simple_ray(ds, start_position=ds.domain_left_edge, end_position=ds.domain_right_edge, lines=['H'])
        assert_almost_equal(ray.r['redshift'][0], 6.99900695e-03, decimal=8)
        assert_almost_equal(ray.r['redshift'][-1], -1.08961751e-02, decimal=8)

    def test_light_ray_non_cosmo_from_dataset(self):
        """
        This test generates a non-cosmological light ray created from an already
        loaded dataset
        """
        ds = load(COSMO_PLUS_SINGLE)
        lr = LightRay(ds)

        ray_start = [0,0,0]
        ray_end = [1,1,1]
        lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                          fields=['temperature', 'density', 'H_p0_number_density'],
                          data_filename='lightray.h5')

        ds = load('lightray.h5')
        compare_light_ray_solutions(lr, ds)

    def test_light_ray_redshift_coverage(self):
        """
        Tests to assure a light ray covers the full redshift range appropriate
        for that comoving line of sight distance.  Was not always so!
        """
        ds = load(GIZMO_COSMO_SINGLE)
        ray = make_simple_ray(ds, start_position=ds.domain_left_edge, end_position=ds.domain_right_edge, lines=['H'])
        assert_almost_equal(ray.r['redshift'][0], 0.00489571, decimal=8)
        assert_almost_equal(ray.r['redshift'][-1], -0.00416831, decimal=8)

    def test_light_ray_redshift_monotonic(self):
        """
        Tests to assure a light ray redshift decreases monotonically
        when ray extends outside the domain.
        """
        ds = load(COSMO_PLUS_SINGLE)
        ray = make_simple_ray(ds, start_position=ds.domain_center,
                              end_position=ds.domain_center+ds.domain_width)
        assert((np.diff(ray.data['redshift']) < 0).all())
