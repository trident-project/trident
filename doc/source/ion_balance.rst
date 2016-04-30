.. _ion-balance:

Using Ion Balance
============================

In addition to being able to create absorption spectra for atomic line transitions of species that aren't explicitly tracked in the simulation, the ion_balance package contained within Trident can be used to generate entire species fields that can be analyzed using the standard yt analysis packages.  This page provides some examples as to how these fields can be generated and analyzed.

Generating species fields
_________________________

As always, we first need to import yt and Trident and then load up a dataset:

.. code-block:: python

   import yt
   import trident

   file_name = 'enzo_cosmology_plus/RD0009/RD0009'
   ds = yt.load(file_name)

We also need to specify the ionization look-up table that we wish to use for producing the species values, you should have acquired one when you first installed Trident, for this example, we'll plan on using the high resolution file.  We can then use the table to generate a new species field, e.g. OVI (five-times-ionized oxygen).

.. code-block:: python

   ionization_table = 'hm2012_hr.h5'

   tri.add_ion_number_density("O", 6, ionization_table, ds)

To show how one can use this newly generated field, we'll make a project of the OVI number density field to produce a OVI column density in the simulation.

.. code-block:: python

   proj = yt.ProjectionPlot(ds, "z", "O_p5_number_density")
   proj.save()

which produces:

.. image:: <link to image here>

We can also generate a field to represent the mass of OVI in the simulation and create a phase plot to show there the OVI mass lives as a function of density and temperature.

.. code-block:: python

   tri.add_ion_mass_field("O", 6, ionization_table, ds)

   # we need to extract a data source from the dataset to make a phase plot
   ad = ds.all_data()

   phase = yt.PhasePlot(ad, "density", "temperature",
                        ["O_p5_mass"], weight_field="O_p5_mass",
                        fractional=True)
   phase.save()

resulting in:

.. image:: <link to image here>
