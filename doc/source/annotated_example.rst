.. _annotated-example:

Annotated Example
=================

Making spectra with Trident is intended to require a minimal amount of effort 
by the user.  This section will walk you through the steps necessary to 
produce a synthetic spectrum based on simulation data and to view its path
through the parent dataset.

The basic process requires three main steps:

    1. Generate a :class:`~trident.LightRay` from the simulation data 
       representing a sightline through the data.
    2. Define the desired spectrum features and use the ``LightRay`` to 
       create a corresponding synthetic spectrum.
    3. Create a projected image and overplot the path of the ``LightRay``.

.. note::

    The actual code contained in the snippets below can be found online 
    `here <https://bitbucket.org/trident-project/trident/src/tip/examples/working_script.py>`_,
    or locally in ``trident.path/examples/working_script.py``

LightRay Generation
--------------------

A :class:`~trident.LightRay` is a 1D object representing the path a ray of
light takes through a simulation volume on its way from some bright background
object to the observer.  It records all of the gas fields it intersects along
the way for use in construction of a spectrum.  

In order to generate a ``LightRay`` from your data, you need to first make sure 
that you've imported both the yt and Trident packages::

   import yt
   import trident

Then, you need to specify from which dataset to extract the light ray.  For 
this example, we'll be using a `sample Enzo dataset 
<http://yt-project.org/data/>`_ but the process should work for other data 
formats and simulation types::

   fn = 'enzo_cosmology_plus/RD0009/RD0009'

We need to decide the trajectory that the :class:`~trident.LightRay` will take
through our simulation volume.  This arbitrary trajectory is specified with
coordinates in code length units (e.g. [x_start, y_start, z_start] to 
[x_end, y_end, z_end]). Probably the simplest trajectory is cutting
diagonally from the origin of the simulation volume to its outermost corner
using the yt ``domain_left_edge`` and ``domain_right_edge`` attributes.  Here
we load the dataset into yt to get access to these attributes::

    ds = yt.load(fn)
    ray_start = ds.domain_left_edge
    ray_end = ds.domain_right_edge

Let's define what lines or species we want to be added to our final spectrum.
In this case, we want to deposit all hydrogen, carbon, nitrogen, oxygen,
and magnesium lines to the resulting spectrum from the dataset::

    line_list = [‘H’, 'C', 'N', "O', 'Mg']

We can now generate the light ray using the :class:`~trident.make_simple_ray`
function by passing the dataset and the trajectory endpoints to it as well
as telling trident to save the resulting ray dataset to an HDF5 file. We
explicitly instruct trident to pull all necessary fields from the dataset
in order to be able to add the lines from our ``line_list``.
Lastly, we set the ``ftype`` keyword as the field type of the fields
where Trident will look to find density, temperature, and metallicity for
building the required ion fields::

    ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,
                                  data_filename="ray.h5",
                                  lines=line_list,
                                  ftype='gas')

.. warning::
    It is imperative that you set the `ftype` keyword properly for your dataset.
    An ``ftype`` of 'gas' is adequate for grid-based codes, but not particle.
    Particle-based datasets must set ``ftype`` to the field type
    of their gas particles (e.g. ``PartType0``) to assure that Trident builds 
    the ion fields on the particles themselves before smoothing these fields 
    to the grid.  By not setting this correctly, you risk bad ion values by
    building from smoothed gas fields.

Overlaying a LightRay's Trajectory on a Projection
--------------------------------------------------

Here we create a projection of the density field along the x axis of the 
dataset, and then overplot the path the LightRay takes through the simulation,
before saving it to disk::

    p = yt.ProjectionPlot(ds, ‘x’, ‘density’)
    p.annotate_ray(ray, arrow=True)
    p.save(‘projection.png’)

.. image:: http://trident-project.org/data/doc_images/annotated_example/projection.png

Spectrum Generation
-------------------

Now that we have our ``LightRay``, we can use it to generate a spectrum.
To create a spectrum, we need to make a :class:`~trident.SpectrumGenerator`
object defining our desired wavelength range and bin size.  You can do this
by manually setting these features, or just using one of the presets for 
an instrument.  Currently, we have three pre-defined instruments, the G130M,
G160M, and G140L observing modes for the Cosmic Origins Spectrograph aboard
the Hubble Space Telescope: ``COS-G130M``, ``COS-G160M``, and ``COS-G140L``.
Notably, instrument ``COS`` aliases to ``COS-G130M``.

We then use this :class:`~trident.SpectrumGenerator` object to make a *raw* 
spectrum according to the intersecting fields it encountered in the 
corresponding :class:`~trident.LightRay`.  We save this spectrum to disk, and
plot it::

    sg = trident.SpectrumGenerator('COS-G130M')
    sg.make_spectrum(ray, lines=line_list)
    sg.save_spectrum('spec_raw.txt')
    sg.plot_spectrum('spec_raw.png')

.. image:: http://trident-project.org/data/doc_images/annotated_example/spec_raw.png
   :width: 700

From here we can do some post-processing to the spectrum to include 
additional features that would be present in an actual observed spectrum.
We add a background quasar spectrum, a Milky Way foreground, apply the
COS line spread function, and add gaussian noise with SNR=30::

    sg.add_qso_spectrum()
    sg.add_milky_way_foreground()
    sg.apply_lsf()
    sg.add_gaussian_noise(30)

Finally, we use plot and save the resulting spectrum to disk::

    sg.save_spectrum('spec_final.txt')
    sg.plot_spectrum('spec_final.png')

which produces:

.. image:: http://trident-project.org/data/doc_images/annotated_example/spec_final.png
   :width: 700

To create more complex or ion-specific spectra, refer to :ref:`advanced-spectra`.
