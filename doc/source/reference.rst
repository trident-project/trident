.. _api-reference:

API Reference
=============

Generating Rays
---------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   ~trident.ray_generator.make_simple_ray
   ~trident.ray_generator.make_compound_ray
   ~trident.light_ray.LightRay

Generating Spectra
------------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   ~trident.spectrum_generator.SpectrumGenerator
   ~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum
   ~trident.instrument.Instrument
   ~trident.lsf.LSF
   ~trident.line_database.Line
   ~trident.line_database.LineDatabase

Plotting Spectra
----------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   ~trident.spectrum_generator.load_spectrum
   ~trident.spectrum_generator.plot_spectrum

Adding Ion Fields
-----------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   ~trident.ion_balance.add_ion_fields
   ~trident.ion_balance.add_ion_fraction_field
   ~trident.ion_balance.add_ion_number_density_field
   ~trident.ion_balance.add_ion_density_field
   ~trident.ion_balance.add_ion_mass_field

Miscellaneous Utilities
-----------------------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   ~trident.utilities.make_onezone_dataset
   ~trident.utilities.make_onezone_ray
   ~trident.roman.to_roman
   ~trident.roman.from_roman
   ~trident.config.trident_path
   ~trident.config.trident
   ~trident.config.verify
   ~trident.absorption_spectrum.absorption_spectrum_fit.generate_total_fit
