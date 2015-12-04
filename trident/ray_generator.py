"""
SpectrumGenerator class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os

from yt.analysis_modules.cosmological_observation.api import \
    LightRay
from yt.convenience import \
    load
from yt.funcs import \
    mylog, \
    YTArray

def make_simple_ray(parameter_filename, start_position, end_position,
                    fields=None, solution_filename=None, data_filename=None, 
                    trajectory=None, redshift=None, 
                    setup_function=None, load_kwargs=None):
    """
    A helper function to yt's LightRay class simplified for use with single 
    dataset outputs.

    **Parameters**

    """
    if fields is None:
        fields = ['density', 'temperature', 'metallicity']
    lr = LightRay(parameter_filename, load_kwargs=load_kwargs)
    return lr.make_light_ray(start_position=start_position,
                             end_position=end_position,
                             trajectory=trajectory,
                             fields=fields, 
                             setup_function=setup_function,
                             solution_filename=solution_filename,
                             data_filename=data_filename,
                             redshift=redshift)

def make_compound_ray(parameter_filename, simulation_type,
                      start_redshift, end_redshift,
                      fields=None, solution_filename=None, data_filename=None, 
                      use_minimum_datasets=True, deltaz_min=0.0,
                      minimum_coherent_box_fraction=0.0, seed=None, 
                      setup_function=None, load_kwargs=None):

    """
    A helper function to yt's LightRay class simplified for use with multiple
    dataset outputs from a single simulation.

    **Parameters**

    """
    if fields is None:
        fields = ['density', 'temperature', 'metallicity']
    lr = LightRay(parameter_filename, 
                  simulation_type=simulation_type, 
                  near_redshift=start_redshift, 
                  far_redshift=end_redshift,
                  use_minimum_datasets=use_minimum_datasets,
                  deltaz_min=deltaz_min,
                  minimum_coherent_box_fraction=minimum_coherent_box_fraction,
                  load_kwargs=load_kwargs)

    return lr.make_light_ray(seed=seed, 
                             fields=fields, 
                             setup_function=setup_function,
                             solution_filename=solution_filename,
                             data_filename=data_filename,
                             redshift=None, njobs=-1)
