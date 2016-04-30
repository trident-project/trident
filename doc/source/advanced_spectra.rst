.. _advanced-spectra:

Advanced Spectrum Generation
============================

Outside of generating a basic "observed" spectrum with the standard features that might present in a realistic asborption spectrum, the user can also customize the generated spectrum to have an arbitrary set of desired properties, e.g. wavelength range, spectral resolution, and included absorption lines. The following code goes through the process of setting such properties and displays the resulting spectra.

For this demostation, we'll be using a light ray that we've previously generated from a RAMSES dataset.  If you'd like to try to reproduce the spectra included below you can get the ray from the Trident project space using:

.. highlight:: none

::

   $ wget http://trident-project.org/data/sample_rays/ramses_ray.h5

Now, we'll load up a ray using yt.

.. code-block:: python

   import yt
   import trident

   ray = yt.load("ramses_ray.h5")

We then generate a spectrum that goes from 1150 angstroms to 1250 angstroms with a resolution of 0.01 angstroms.

.. code-block:: python

   sg = tri.SpectrumGenerator(lambda_min=1150, lambda_max=1250, dlambda=0.01)

From here, we can pass the ray to the SpectrumGenerator object to use in the construction of a spectrum.  As a first pass, we'll create a spectrum that just include the lines produced by hydrogen.

.. code-block:: python

    sg.make_spectrum(ray, lines=["H"])
    sg.plot_spectrum('spec_H.png')

The resulting spectrum contains a nice, big Lyman-alpha feature.

.. image:: http://trident-project.org/data/doc_images/spectra/spec_H.png

.. note::
    If you do try recreating this spectrum yourself, be aware that since this deep Lyman-alpha feature is the combination of many combined absorption features, the generate of the spectrum may take several minutes.

If, instead, we want to shows the lines that would be in our spectral range due to carbon, nitrogen, and oxygen, we can do the following.

.. code-block:: python

    sg.make_spectrum(ray, lines=["C", "N", "O"])
    sg.plot_spectrum('spec_CNO.png')

And now we have:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_CNO.png

We can see how these two spectra combined when we include all of the same lines:

.. code-block:: python

    sg.make_spectrum(ray, lines=["H", "C", "N", "O"])
    sg.plot_spectrum('spec_HCNO.png')

which gives:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_HCNO.png

We can get even more specific, by generating a spectrum that only contains lines due to a single ion species.  For example, we might just want the lines from four-times-ionized nitrogen, N V:

.. code-block:: python

    sg.make_spectrum(ray, lines=["N V"])
    sg.plot_spectrum('spec_NV.png')

This spectrum only shows a couple of small lines on the right hand side.

.. image:: http://trident-project.org/data/doc_images/spectra/spec_NV.png

But if that level of specificity isn't enough, we can request individual lines!

.. code-block:: python

    sg.make_spectrum(ray, lines=["C I 1193", "C I 1194"])
    sg.plot_spectrum('spec_CI_1193_1194.png')

And we end up with:

.. image:: http://trident-project.org/data/doc_images/spectra/spec_CI_1193_1194.png

To understand how to further customize your spectra, refer to the class and function documentation provided in the API Reference section.
