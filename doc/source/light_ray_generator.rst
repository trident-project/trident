.. _light-ray-generator:

Light Ray Generator
===================

The :class:`~trident.light_ray.LightRay` is the one-dimensional object representing
the pencil beam of light traveling from the source to the observer. Light
rays can stack multiple datasets together to span a redshift interval
larger than the simulation box.

.. image:: _images/lightray.png

The preferred manner for generating rays uses the
:func:`~trident.ray_generator.make_simple_ray` for
:class:`~trident.light_ray.LightRay` 's spanning a single dataset
and
:func:`~trident.ray_generator.make_compound_ray` for
:class:`~trident.light_ray.LightRay` 's spanning multiple datasets.

Simple Rays
-----------

For a simple ray, you specify the dataset to use, the start and end coordinates
of your 1D line, and then optionally any additional fields you want stored on the
:class:`~trident.light_ray.LightRay` or optionally any ionic species you will want to
use with this ray::

    import yt
    import trident
    ds = yt.load('FIRE_M12i_ref11/snapshot_600.hdf5')
    ray = trident.make_simple_ray(ds,
                                  start_position=[0, 0, 0],
                                  end_position=[60000, 60000, 60000],
                                  lines=['H', 'Mg', 'O'],
                                  fields=[('gas', 'temperature'), ('gas', 'metallicity')])

Compound Rays
-------------

For a compound ray, you specify the simulation parameter filename, the simulation code,
the start and end redshifts of the :class:`~trident.light_ray.LightRay`, and optionally
any additional fields you want stored or any ionic species you will want to use with
this ray::

    import trident
    fn = 'enzo_cosmology_plus/AMRCosmology.enzo'
    ray = trident.make_compound_ray(fn,
                                    simulation_type='Enzo',
                                    near_redshift=0.0,
                                    far_redshift=0.1,
                                    lines=['H', 'Mg', 'O'],
                                    fields=[('gas', 'temperature'), ('gas', 'metallicity')])

Ray Fields
----------

The resulting ``ray`` is a :class:`~trident.light_ray.LightRay` object, consisting of a series
of arrays representing the different fields it probes in the original dataset along
its length.  Each element in the arrays represents a different resolution element
along the path of the ray.  The ray also possesses some special fields not originally
present in the original dataset:

    * ``('gas', 'l')`` Location along the LightRay length from 0 to 1.
    * ``('gas', 'dl')`` Pathlength of resolution element (not a *true* pathlength for particle-based codes)
    * ``('gas', 'redshift')`` Cosmological redshift of resolution element
    * ``('gas', 'redshift_dopp')`` Doppler redshift of resolution element
    * ``('gas', 'redshift_eff')`` Effective redshift (combined cosmological and Doppler)

Like any dataset, you can see what fields are present on the ray by examining its
``derived_field_list`` (e.g., ``print(ds.derived_field_list``).  If you want more ions
present on this ray than are currently available, you can add them with
:class:`~trident.ion_balance.add_ion_fields` (see: :ref:`ion-balance`).

This ``ray`` object is also saved to disk as an HDF5 file, which can later be loaded
into ``yt`` as a stand-alone dataset.  By default it is saved as ``ray.h5``, but you
can specify other filenames when you create it.  To later access this file and
load it into yt, load it like any other dataset: ``ds = yt.load('ray.h5')``.

Calculating Column Densities
----------------------------

Perhaps we wish to know the total column density of a particular ion present along
this :class:`~trident.light_ray.LightRay`. This can easily be done by multiplying the desired
ion number density field by the pathlength field, ``dl``, to yield an array of
column densities for each resolution element, and then summing them together::

    column_density_HI = ray.r[('gas', 'H_p0_number_density')] * ray.r[('gas', 'dl')]
    print('HI Column Density = %g' % column_density_HI.sum())

Examining LightRay Solution Data
--------------------------------

When a :class:`~trident.light_ray.LightRay` is created, it saves the source information
from the dataset that produced it in a dictionary, including its filename, its start
and end points in the original dataset, etc.  This is all accessible when
you load up the :class:`~trident.light_ray.LightRay` again through the
``light_ray_solution``::

    import yt
    ds = yt.load('ray.h5')
    print(ds.light_ray_solution)

    [{'end': unyt_array([1., 1., 1.], 'unitary'),
    'filename': 'snapshot_600.hdf5',
    'redshift': 0.05,
    'start': unyt_array([0.48810148, 0.51748806, 0.54316002], 'unitary'),
    'traversal_box_fraction': unyt_quantity(0.83878521, 'unitary'),
    'unique_identifier': '1436307563512020127'}]

Useful Tips for Making Compound LightRays
-----------------------------------------

Below are some tips that may come in handy for creating proper LightRays.  For full
use of these, you may have to create the :class:`~trident.light_ray.LightRay`
by hand instead of using the :func:`~trident.ray_generator.make_compound_ray` helper
function.

How many snapshots do I need for a compound ray?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The number of snapshots required to traverse some redshift interval depends
on the simulation box size and cosmological parameters.  Before running an
expensive simulation only to find out that you don't have enough outputs
to span the redshift interval you want, have a look at the guide
`Planning Simulations for LightCones or LightRays
<https://yt-astro-analysis.readthedocs.io/en/latest/planning_cosmology_simulations.html>`__.
The functionality described there will allow you to calculate the precise
number of snapshots and specific redshifts at which they should be written.

My snapshots are too far apart!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``max_box_fraction`` keyword, provided when creating the Lightray,
allows the user to control how long a ray segment can be for an
individual dataset.  Be default, the LightRay generator will try to
make segments no longer than the size of the box to avoid sampling the
same structures more than once.  However, this can be increased in the
case that the redshift interval between datasets is longer than the
box size.  Increasing this value should be done with caution as longer
ray segments run a greater risk of coming back to somewhere near their
original position.

What if I have a zoom-in simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A zoom-in simulation has a high resolution region embedded within a
larger, low resolution volume.  In this type of simulation, it is likely
that you will want the ray segments to stay within the high resolution
region.  To do this, you must first specify the size of the high
resolution region when creating the LightRay using the
``max_box_fraction`` keyword.  This will make sure that
the calculation of the spacing of the segment datasets only takes into
account the high resolution region and not the full box size.  If your
high resolution region is not a perfect cube, specify the smallest side.
Then, in the call to
:func:`~trident.light_ray.LightRay.make_light_ray`,
use the ``left_edge`` and ``right_edge`` keyword arguments to specify the
precise location of the high resolution region.

Technically speaking, the ray segments should no longer be periodic
since the high resolution region is only a sub-volume within the
larger domain.  To make the ray segments non-periodic, set the
``periodic`` keyword to False.  The LightRay generator will continue
to generate randomly oriented segments until it finds one that fits
entirely within the high resolution region.  If you have a high
resolution region that can move and change shape slightly as structure
forms, use the `min_level` keyword to mandate that the ray segment only
pass through cells that are refined to at least some minimum level.

If the size of the high resolution region is not large enough to
span the required redshift interval, the `LightRay` generator can
be configured to treat the high resolution region as if it were
periodic simply by setting the ``periodic`` keyword to True.  This
option should be used with caution as it will lead to the creation
of disconnected ray segments within a single dataset.

I want a continous trajectory over the entire ray.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set the ``minimum_coherent_box_fraction`` keyword argument to a very
large number, like infinity (``numpy.inf``).
