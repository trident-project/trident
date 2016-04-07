"""
SpectrumGenerator class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
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
                    fields=None, solution_filename=None, 
                    data_filename=None, trajectory=None, redshift=None, 
                    setup_function=None, load_kwargs=None):
    """
    A helper function to create a yt LightRay object for single dataset 
    outputs.  Returns a LightRay object for later use with Trident's
    SpectrumGenerator.


    **Parameters**

    parameter_filename : string
        The path to the simulation parameter file or dataset

    start_position, end_position : list of floats
        The coordinates of the starting and ending position of the ray in 
        code_length units.

    fields : optional, list of strings
        The list of which fields to store in the output LightRay.  If none
        selected, defaults to ['density', 'temperature', 'metallicity'].
        Default: None

    solution_filename : optional, string
        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    data_filename : optional, string
        Output filename for ray data stored as an HDF5 file.  Note that 
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    trajectory : optional, list of floats
        The (r, theta, phi) direction of the LightRay.  Use either end_position
        or trajectory, but not both.
        Default: None

    redshift : optional, float
        Sets the highest cosmological redshift of the ray.  By default, it will
        use the cosmological redshift of the dataset, if set, and if not set,
        it will use a redshift of 0.
        Default: None

    setup_function : optional, function
        A function that will be called on the dataset as it is loaded but 
        before the LightRay is generated.  Very useful for adding derived 
        fields and other manipulations of the dataset prior to LightRay 
        creation.
        Default: None

    load_kwargs : optional, dict
        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None
    """
    if fields is None:
        fields = ['density', 'temperature', 'metallicity']
    if data_filename is None:
        data_filename = 'ray.h5'
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
                      near_redshift, far_redshift,
                      fields=None, solution_filename=None, data_filename=None, 
                      use_minimum_datasets=True, deltaz_min=0.0,
                      minimum_coherent_box_fraction=0.0, seed=None, 
                      setup_function=None, load_kwargs=None):

    """
    A helper function to create yt LightRay objects spanning multiple dataset 
    outputs from a single simulation.  Returns a LightRay object for later use 
    with Trident's SpectrumGenerator.  Note that this functionality only 
    currently works with Enzo and Gadget yt frontends.  For inclusion of 
    additional simulation frontends, please contact the Trident developers.
    Returns a yt LightRay object.


    **Parameters**

    parameter_filename : string
        The path to the simulation parameter file or dataset

    simulation_type : string
        The simulation type of the parameter file.  At present, this
        functionality only works with "Enzo" and "Gadget" yt frontends.

    near_redshift, far_redshift : floats
        The near and far redshift bounds of the LightRay through the 
        simulation datasets.

    fields : optional, list of strings
        The list of which fields to store in the output LightRay.  If none
        selected, defaults to ['density', 'temperature', 'metallicity'].
        Default: None

    solution_filename : optional, string
        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    data_filename : optional, string
        Output filename for ray data stored as an HDF5 file.  Note that 
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    use_minimum_datasets : optional, bool
        Use the minimum number of datasets to make the ray continuous 
        through the supplied datasets from the near_redshift to the 
        far_redshift.  If false, the LightRay solution will contain as many
        datasets as possible to enable the light ray to traverse the 
        desired redshift interval.
        Default: True

    deltaz_min : optional, float
        The minimum delta-redshift value between consecutive datasets used
        in the LightRay solution.
        Default: 0.0

    minimum_coherent_box_fraction : optional, float
        When use_minimum_datasets is set to False, this parameter specifies
        the fraction of the total box width to be traversed before 
        rerandomizing the ray location and trajectory.
        Default: 0.0

    seed : optional, int
        Sets the seed for the random number generator used to determine the
        location and trajectory of the LightRay as it traverses the 
        simulation datasets.  For consistent results between LightRays, 
        use the same seed value.
        Default: None

    setup_function : optional, function
        A function that will be called on the dataset as it is loaded but 
        before the LightRay is generated.  Very useful for adding derived 
        fields and other manipulations of the dataset prior to LightRay 
        creation.
        Default: None

    load_kwargs : optional, dict
        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    **Parameters**

    """
    if fields is None:
        fields = ['density', 'temperature', 'metallicity']
    lr = LightRay(parameter_filename, 
                  simulation_type=simulation_type, 
                  near_redshift=near_redshift, 
                  far_redshift=far_redshift,
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
