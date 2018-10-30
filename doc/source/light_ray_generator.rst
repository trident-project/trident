.. _light-ray-generator:

Light Ray Generator
===================

The :class:`~trident.LightRay` is the one-dimensional object representing
the pencil beam of light traveling from the source to the observer. Light
rays can stack multiple datasets together to span a redshift interval
larger than the simulation box.

.. image:: _images/lightray.png

A ray segment records the information of all grid cells intersected by the
ray as well as the path length, ``dl``, of the ray through the cell.  Column
densities can be calculated by multiplying physical densities by the path
length.

.. _multi-ray:

Configuring the Light Ray Generator
-----------------------------------

Below follows the creation of a light ray from multiple datasets stacked
together.  The primary Trident interface to this is covered in
:ref:`compound-ray`.  A light ray can also be made from a single dataset.
For information on this, see :ref:`single-ray`.

The arguments required to instantiate a :class:`~trident.LightRay` are
the simulation parameter file, the simulation type, the nearest redshift,
and the furthest redshift.

.. code-block:: python

  from trident import LightRay
  lr = LightRay("enzo_tiny_cosmology/32Mpc_32.enzo",
                simulation_type="Enzo",
                near_redshift=0.0, far_redshift=0.1)

Making Light Ray Data
---------------------

Once the LightRay object has been instantiated, the
:func:`~trident.light_ray.LightRay.make_light_ray`
function will trace out the rays in each dataset and collect information for all the
fields requested.  The output file will be an yt-loadable dataset containing all the
cell field values for all the cells that were intersected by the ray.  A
single LightRay object can be used over and over to make multiple
randomizations, simply by changing the value of the random seed with the
``seed`` keyword.

.. code-block:: python

  lr.make_light_ray(seed=8675309,
                    fields=['temperature', 'density'],
                    use_peculiar_velocity=True)

  # Optionally, we can now overplot the part of this ray that intersects
  # one output from the source dataset in a ProjectionPlot
  ds = yt.load('enzo_tiny_cosmology/RD0004/RD0004')
  p = yt.ProjectionPlot(ds, 'z', 'density')
  p.annotate_ray(lr)
  p.save()

.. _single-ray:

Light Rays Through Single Datasets
----------------------------------

LightRays can also be configured for use with single datasets. In this case,
one must specify the ray's trajectory explicitly.  The main Trident interface
to this functionality is covered in :ref:`simple-ray`.

.. code-block:: python

  from trident import LightRay
  import yt

  ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
  lr = LightRay(ds)

  lr.make_light_ray(start_position=ds.domain_left_edge,
                    end_position=ds.domain_right_edge,
                    solution_filename='lightraysolution.txt',
                    data_filename='lightray.h5',
                    fields=['temperature', 'density'])

  # Overplot the ray on a projection.
  p = yt.ProjectionPlot(ds, 'z', 'density')
  p.annotate_ray(lr)
  p.save()

Alternately, the ``trajectory`` keyword can be used in place of `end_position`
to specify the (r, theta, phi) direction of the ray.

Useful Tips for Making LightRays
--------------------------------

Below are some tips that may come in handy for creating proper LightRays.

How many snapshots do I need?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
