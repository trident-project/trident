.. _advanced-spectra:

Advanced Spectrum Generation
============================

In addition to generating a basic spectrum as demonstrated in 
the :ref:`annotated example <annotated-example>`, the user can also 
customize the generated spectrum in a variety of ways.  One can choose which
spectral lines to deposit or choose different settings for the characteristics
of the spectrograph, and more.  The following code goes through the process of 
setting these properties and shows what impact it has on resulting spectra.

For this demonstation, we'll be using a light ray passing through a very dense
disk of gas, taken from the initial output from the AGORA isolated box 
simulation using ART-II in `Kim et al. (2016)
<http://adsabs.harvard.edu/abs/2016ApJ...833..202K>`_.
If you'd like to try to reproduce the spectra included below you can get 
the :class:`~trident.LightRay` file from the Trident sample data using the 
command:

.. highlight:: none

::

   $ wget http://trident-project.org/data/sample_data/ART-II_ray.h5

Now, we'll load up the ray using yt::

   import yt
   import trident
   ray = yt.load('ART-II_ray.h5')

Setting the spectrograph
------------------------

Let's set the characteristics of the spectrograph we will use to create
this spectrum.  We can either choose the wavelength range and resolution
and line spread function explicitly, or we can choose one of the preset
instruments that come with Trident.  To list the presets and their respective
values, use this command::

    print(trident.valid_instruments)
    
Currently, we have `three settings for the Cosmic Origins Spectrograph 
<http://www.stsci.edu/hst/cos/design/gratings/>`_ available:
``COS-G130M``, ``COS-G140L``, and ``COS-G160M``, but we plan to add more
instruments soon.  To use one of them, we just use the name string in the 
SpectrumGenerator class::

   sg = trident.SpectrumGenerator('COS-G130M')

But instead, let's just set our wavelength range manually
from 1150 angstroms to 1250 angstroms with a resolution of 0.01 angstroms::

   sg = trident.SpectrumGenerator(lambda_min=1150, lambda_max=1250, dlambda=0.01)

From here, we can pass the ray to the :class:`~trident.SpectrumGenerator` object 
to use in the construction of a spectrum. 

Choosing what absorption features to include
--------------------------------------------

There is a :class:`~trident.LineDatabase` class that controls which spectral 
lines you can add to your spectrum.  Trident provides you with a default
:class:`~trident.LineDatabase` with 213 spectral lines commonly used in CGM 
and IGM studies, but you can create your own :class:`~trident.LineDatabase` 
with different lines.  To see a list of all the lines included in the default 
line list::

    ldb = trident.LineDatabase('lines.txt')
    print(ldb)

which is reading lines from the 'lines.txt' file present in the 
data directory (see :ref:`where is Trident installed? <where-installed>`)
We can specify any subset of these spectral lines to use when creating the 
spectrum from our master line list.  So if you're interested in just looking 
at neutral hydrogen lines in your spectrum, you can see what lines will be 
included with the command::

    print(ldb.parse_subset('H I'))

As a first pass, we'll create a spectrum that just include lines produced 
by hydrogen::

    sg.make_spectrum(ray, lines=['H'])
    sg.plot_spectrum('spec_H.png')

The resulting spectrum contains a nice, big Lyman-alpha feature.

.. image:: http://trident-project.org/data/doc_images/spectra/spec_H.png

If, instead, we want to shows the lines that would be in our spectral range 
due to carbon, nitrogen, and oxygen, we can do the following::

    sg.make_spectrum(ray, lines=['C', 'N', 'O'])
    sg.plot_spectrum('spec_CNO.png')

And now we have:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_CNO.png

We can see how these two spectra combined when we include all of the same 
lines::

    sg.make_spectrum(ray, lines=['H', 'C', 'N', 'O'])
    sg.plot_spectrum('spec_HCNO.png')

which gives:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_HCNO.png

We can get even more specific, by generating a spectrum that only contains 
lines due to a single ion species.  For example, we might just want the 
lines from four-times-ionized nitrogen, N V::

    sg.make_spectrum(ray, lines=['N V'])
    sg.plot_spectrum('spec_NV.png')

This spectrum only shows a couple of small lines on the right hand side.

.. image:: http://trident-project.org/data/doc_images/spectra/spec_NV.png

But if that level of specificity isn't enough, we can request individual lines::

    sg.make_spectrum(ray, lines=['C I 1193', 'C I 1194'])
    sg.plot_spectrum('spec_CI_1193_1194.png')

And we end up with:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_CI_1193_1194.png

Or we can just include all of the available lines in our 
:class:`~trident.LineDatabase` with::

    sg.make_spectrum(ray, lines='all')
    sg.plot_spectrum('spec_all.png')

Giving us:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_all.png

To understand how to further customize your spectra, look at the documentation 
for the :class:`~trident.SpectrumGenerator` and :class:`~trident.LineDatabase`
classes and other :ref:`API <api-reference>` documentation.
