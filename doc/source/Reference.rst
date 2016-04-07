.. api-reference:

API Reference
=============

#.. automodapi::


Interface for Generating Rays
-----------------------------

.. autosummary::
    :toctree: generated/

.. autofunction:: trident.make_simple_ray
.. autofunction:: trident.make_compound_ray
.. autoclass:: trident.LightRay
    :members:

Interface for Generating Spectra
--------------------------------

.. autosummary::
    :toctree: generated/

.. autoclass:: trident.SpectrumGenerator
   :members:
.. autoclass:: trident.Instrument
   :members:
.. autoclass:: trident.Line
   :members:
.. autoclass:: trident.LineDatabase
   :members:
.. autoclass:: trident.LSF
   :members:

Interface for Plotting Spectra
------------------------------

.. autosummary::
    :toctree: generated/

.. autofunction:: trident.plot_spectrum

Interface for Adding Ion Fields
-------------------------------

.. autosummary::
    :toctree: generated/

.. autofunction:: trident.add_ion_density_field
.. autofunction:: trident.add_ion_number_density_field
.. autofunction:: trident.add_ion_fraction_field
.. autofunction:: trident.add_ion_mass_field

Miscellaneous Utilities
-----------------------

.. autosummary::
    :toctree: generated/

.. autofunction:: trident.toRoman
.. autofunction:: trident.fromRoman
.. autofunction:: trident.trident
