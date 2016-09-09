.. _making-spectra:

Making Spectra
==============

Making spectra with Trident is intended to require a minimal amount of effort 
by the user.  This section will walk you through the steps necessary to 
produce a synthetic spectrum based on simulation data.

The basic process requires two main steps:

    1. Generate a :class:`~trident.LightRay` from the simulation data.
    2. Define the desired spectrum features and use the light ray to 
       create the corresponding synthetic spectrum.

.. note::

    The actual code contained in the snippets below can be found online 
    `here <https://bitbucket.org/trident-project/trident/src/tip/examples/working_script.py>`_,
    or locally in ``trident.path/examples/working_script.py``

Light Ray Generation
--------------------

A :class:`~trident.LightRay` is a 1D object representing the path a ray of
light takes through a simulation volume on its way from some bright background
object to the observer.  It records all of the gas fields it intersects along
the way for use in construction of a spectrum.  

In order to generate a light ray from your data, you need to first make sure 
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
[x_end, y_end, z_end]), but probably the simplest trajectory is to go
diagonally from the origin of the simulation volume to its outermost corner
using the yt ``domain_left_edge`` and ``domain_right_edge`` attributes.  Here
we load the dataset into yt to get access to these attributes::

    ds = yt.load(fn)
    ray_start = ds.domain_left_edge
    ray_end = ds.domain_right_edge

We can now generate the light ray using the :class:`~trident.make_simple_ray`
function by passing the dataset and the trajectory endpoints to it as well
as telling trident to save the resulting ray dataset to an HDF5 file. We
explicitly instruct trident to include all necessary fields to the ray
so that we can use the ray to deposit certain lines.  In this case, we want
to deposit all hydrogen, carbon, nitrogen, oxygen, and magnesium lines.  
Lastly, we set the ftype keyword to represent the field type of the fields
to search where we might find these ion fields::

    ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,
                                  data_filename="ray.h5",
                                  lines=['H', 'C', 'N', 'O', 'Mg'],
                                  ftype='gas')

.. note::
    it is also possible to generate a :class:`~trident.LightRay` spanning 
    multiple consecutive datasets to approximate an IGM sightline.  Additional
    documentation will be added soon, but in the meanwhile, take a look at 
    :class:`~trident.make_compound_ray`.

Spectrum Generation
-------------------

Now that we have our light ray, we can use it to generate a spectrum.
To create a spectrum, we need to make a :class:`~trident.SpectrumGenerator`
object defining our desired wavelength range and bin size.  You can do this
by manually setting these features, or just using one of the presets for 
an instrument.  Currently, the only working instrument we have is for COS,
the Cosmic Origins Spectrograph aboard the Hubble Space Telescope, but this
will be supplemented in the future.  We then use this 
:class:`~trident.SpectrumGenerator` to make a *raw* spectrum according to the
intersecting fields it encountered in the corresponding 
:class:`~trident.LightRay`::

    sg = trident.SpectrumGenerator('COS')
    sg.make_spectrum(ray, lines=['H', 'C', 'N', 'O', 'Mg'])

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

.. image:: _images/spec.png
   :width: 700

To create more complex or ion-specific spectra, refer to :ref:`advanced-spectra`.
