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

def make_simple_ray(dataset_file, start_position, end_position,
                    fields=None, solution_filename=None, 
                    data_filename=None, trajectory=None, redshift=None, 
                    setup_function=None, load_kwargs=None):
    """
    Create a yt LightRay object for a single dataset (eg CGM).  This is a 
    wrapper function around yt's LightRay interface to reduce some of the 
    complexity there.  
    
    A simple ray is a straight line passing through a single dataset
    where each gas cell intersected by the line is sampled for the desired
    fields and stored.  Several additional fields are created and stored
    including ``dl`` and ``dredshift`` to represent the path length in space and 
    redshift for each element in the ray, ``v_los`` to represent the line of
    sight velocity along the ray, and ``redshift``, ``redshift_dopp``, and 
    ``redshift_eff`` to represent the cosmological redshift, doppler redshift
    and effective redshift (combined doppler and cosmological) for each
    element of the ray.

    A simple ray is typically specified by its start and end positions in the 
    dataset volume.  Because a simple ray only probes a single output, it
    lacks foreground absorbers between the observer at z=0 and the redshift
    of the dataset that one would naturally encounter.  Thus it is usually
    only appropriate for studying the circumgalactic medium rather than
    the intergalactic medium.
    
    This function can accept a yt dataset already loaded in memory, 
    or it can load a dataset if you pass it the dataset's filename and 
    optionally any load_kwargs or setup_function necessary to load/process it
    properly before generating the LightRay object.
    
    **Parameters**

    :dataset_file: string or yt Dataset object
    
        Either a yt dataset or the filename of a dataset on disk.  If you are
        passing it a filename, consider usage of the ``load_kwargs`` and
        ``setup_function`` kwargs.

    :start_position, end_position: list of floats or YTArray object

        The coordinates of the starting and ending position of the desired 
        ray.  If providing a raw list, coordinates are assumed to be in 
        code length units, but if providing a YTArray, any units can be
        specified.

    :fields: optional, list of strings

        The list of which fields to store in the output LightRay.  If none
        selected, defaults to ['density', 'temperature', 'metallicity'].
        Default: None

    :solution_filename: optional, string

        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    :data_filename: optional, string
    
        Output filename for ray data stored as an HDF5 file.  Note that 
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    :trajectory: optional, list of floats

        The (r, theta, phi) direction of the LightRay.  Use either end_position
        or trajectory, but not both.
        Default: None

    :redshift: optional, float

        Sets the highest cosmological redshift of the ray.  By default, it will
        use the cosmological redshift of the dataset, if set, and if not set,
        it will use a redshift of 0.
        Default: None

    :setup_function: optional, function

        A function that will be called on the dataset as it is loaded but 
        before the LightRay is generated.  Very useful for adding derived 
        fields and other manipulations of the dataset prior to LightRay 
        creation.
        Default: None

    :load_kwargs: optional, dict

        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    **Example**

    Generate a simple ray passing from the lower left corner to the upper 
    right corner through some dataset

    >>> import trident
    >>> import yt
    >>> ds = yt.load('path/to/dataset')
    >>> ray = trident.make_simple_ray(ds, 
    ... start_position=ds.domain_left_edge, end_position=ds.domain_right_edge,
    ... fields=['density', 'temperature', 'metallicity'])
    """
    if fields is None:
        fields = ['density', 'temperature', 'metallicity']
    if data_filename is None:
        data_filename = 'ray.h5'
    lr = LightRay(dataset_file, load_kwargs=load_kwargs)
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
    Create a yt LightRay object for multiple consecutive datasets (eg IGM).  
    This is a wrapper function around yt's LightRay interface to reduce some 
    of the complexity there.  
    
    .. note::
        
        The compound ray functionality has only been implemented for the
        Enzo and Gadget code.  If you would like to help us implement
        this functionality for your simulation code, please contact us
        about this on the mailing list.

    A compound ray is a series of straight lines passing through multiple
    consecutive outputs from a single cosmological simulation to approximate
    a continuous line of sight to high redshift.

    Because a single continuous ray traversing a simulated volume can only 
    cover a small range in redshift space (e.g. 100 Mpc only covers the 
    redshift range from z=0 to z=0.023), the compound ray passes rays through  
    multiple consecutive outputs from the same simulation to approximate the
    path of a single line of sight to high redshift.  By probing all of the
    foreground material out to any given redshift, the compound ray is
    appropriate for studies of the intergalactic medium and circumgalactic 
    medium.

    By default, it selects a random starting location and trajectory in 
    each dataset it traverses, to assure that the same cosmological structures
    are not being probed multiple times from the same direction.  In doing 
    this, the ray becomes discontinuous across each dataset.

    The compound ray requires the parameter_filename of the simulation run.
    This is *not* the dataset filename from a single output, but the parameter
    file that was used to run the simulation itself.  It is in this parameter
    file that the output frequency, simulation volume, and cosmological 
    parameters are described to assure full redshift coverage can be achieved
    for a compound ray.  It also requires the simulation_type of the simulation.

    Unlike the simple ray, which is specified by its start and end positions 
    in the dataset volume, the compound ray requires the near_redshift and
    far_redshift to determine which datasets to use to get full coverage
    in redshift space as the ray propagates from near_redshift to far_redshift.
    
    Like the simple ray produced by :class:`~trident.make_simple_ray`, 
    each gas cell intersected by the LightRay is sampled for the desired
    fields and stored.  Several additional fields are created and stored
    including ``dl`` and ``dredshift`` to represent the path length in space and 
    redshift for each element in the ray, ``v_los`` to represent the line of
    sight velocity along the ray, and ``redshift``, ``redshift_dopp``, and 
    ``redshift_eff`` to represent the cosmological redshift, doppler redshift
    and effective redshift (combined doppler and cosmological) for each
    element of the ray.

    **Parameters**

    :parameter_filename: string

        The simulation parameter file *not* the dataset filename

    :simulation_type: string

        The simulation type of the parameter file.  At present, this
        functionality only works with "Enzo" and "Gadget" yt frontends.

    :near_redshift, far_redshift: floats

        The near and far redshift bounds of the LightRay through the 
        simulation datasets.

    :fields: optional, list of strings

        The list of which fields to store in the output LightRay.  If none
        selected, defaults to ['density', 'temperature', 'metallicity'].
        Default: None

    :solution_filename: optional, string

        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    :data_filename: optional, string

        Output filename for ray data stored as an HDF5 file.  Note that 
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    :use_minimum_datasets: optional, bool

        Use the minimum number of datasets to make the ray continuous 
        through the supplied datasets from the near_redshift to the 
        far_redshift.  If false, the LightRay solution will contain as many
        datasets as possible to enable the light ray to traverse the 
        desired redshift interval.
        Default: True

    :deltaz_min: optional, float

        The minimum delta-redshift value between consecutive datasets used
        in the LightRay solution.
        Default: 0.0

    :minimum_coherent_box_fraction: optional, float

        When use_minimum_datasets is set to False, this parameter specifies
        the fraction of the total box width to be traversed before 
        rerandomizing the ray location and trajectory.
        Default: 0.0

    :seed: optional, int

        Sets the seed for the random number generator used to determine the
        location and trajectory of the LightRay as it traverses the 
        simulation datasets.  For consistent results between LightRays, 
        use the same seed value.
        Default: None

    :setup_function: optional, function

        A function that will be called on the dataset as it is loaded but 
        before the LightRay is generated.  Very useful for adding derived 
        fields and other manipulations of the dataset prior to LightRay 
        creation.
        Default: None

    :load_kwargs: optional, dict

        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    **Example**

    Generate a compound ray passing from the redshift 0 to redshift 0.05
    through a multi-output simulation.

    >>> import trident
    >>> import yt
    >>> ds = yt.load('path/to/simulation/parameter/file')
    >>> ray = trident.make_compound_ray(fn, simulation_type='Enzo',
    ... near_redshift=0.0, far_redshift=0.05,
    ... fields=['density', 'temperature', 'metallicity'])
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
