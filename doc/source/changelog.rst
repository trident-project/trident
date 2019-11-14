.. _changelog:

Changelog
=========

This document summarizes changes to the codebase with time for Trident.

Contributors
------------

The `CREDITS file <https://github.com/trident-project/trident/blob/master/CREDITS>`_
has an updated list of contributors to the codebase.

Version 1.2.2 (November 14, 2019)
---------------------------------

This is a bug fix release.

Bug Fixes
^^^^^^^^^

- Shift wavelength of velocity center to redshift from light ray solution
  (`PR 102 <https://github.com/trident-project/trident/pull/102>`__)

Version 1.2.1 (October 1, 2019)
-------------------------------

This is a bug fix release.

Bug Fixes
^^^^^^^^^

- Logging info doesn't use correct units
  (`PR 99 <https://github.com/trident-project/trident/pull/99>`__)

Version 1.2 (September 19, 2019)
--------------------------------

New Features
^^^^^^^^^^^^

- Add support for creating spectra in velocity space
  (`PR 90 <https://github.com/trident-project/trident/pull/90>`__)
- Add ability to set wavelength limits to 'auto'
  (`PR 87 <https://github.com/trident-project/trident/pull/87>`__)
- Add ability to make spectrum from data container
  (`PR 91 <https://github.com/trident-project/trident/pull/91>`__)

Bug Fixes
^^^^^^^^^

- off by one error in subgrid index calculation
  (`PR 85 <https://github.com/trident-project/trident/pull/85>`__)
- fixing error in _ion_mass
  (`PR 81 <https://github.com/trident-project/trident/pull/81>`__)
- H_p0_number_density to default, H_number_density to alias
  (`PR 78 <https://github.com/trident-project/trident/pull/78>`__)
- Implements a better way of calculating the redshift on the light ray
  (`PR 71 <https://github.com/trident-project/trident/pull/71>`__)
- Assures that ion_fraction field reflects on-disk fields
  (`PR 64 <https://github.com/trident-project/trident/pull/64>`__)
- Fix atomic data for Si II 1260 according to Morton 2003
  (`PR 43 <https://github.com/trident-project/trident/pull/43>`__)
- A check to avoid failure when no continuum absorbers were found in the ray
  (`PR 39 <https://github.com/trident-project/trident/pull/39>`__)
- Auto-addition of H_nuclei_density to a LightRay when present in the base dataset
  (`PR 39 <https://github.com/trident-project/trident/pull/39>`__)
- Adding max_box_fraction as a kwarg to make_compound_ray
  (`PR 37 <https://github.com/trident-project/trident/pull/37>`__)
- Updated trident_path() to be OS Independent
  (`PR 36 <https://github.com/trident-project/trident/pull/36>`__)
- simplify setting up ion fields using the "local" field type
  (`PR 30 <https://github.com/trident-project/trident/pull/30>`__)
- split and join filenames using os.sep instead of assuming unix
  (`PR 29 <https://github.com/trident-project/trident/pull/29>`__)
- updated oscillator strengths and gamma's for Si II 1206 and Si III 1260
  (`PR 25 <https://github.com/trident-project/trident/pull/25>`__)

Minor Enhancements
^^^^^^^^^^^^^^^^^^

- Calculating LOS velocity with relative velocities to account for bulk motion
  (`PR 93 <https://github.com/trident-project/trident/pull/93>`__)
- Enabling use of `output_absorbers_file` kwarg in SpectrumGenerator
  (`PR 58 <https://github.com/trident-project/trident/pull/58>`__)
- Switching imports from yt.analysis_modules to yt_astro_analysis
  (`PR 55 <https://github.com/trident-project/trident/pull/55>`__)
- Enable passing in-memory LineDatabase to SpectrumGenerator
  (`PR 42 <https://github.com/trident-project/trident/pull/42>`__)
- Added equivalent width calculation to line_observables_dict
  (`PR 40 <https://github.com/trident-project/trident/pull/40>`__)
- Numerous documentation updates
- Updates and fixes to testing

Version 1.1 (November 18, 2017)
-------------------------------

- Trident development has changed from mercurial to git, and the source has
  moved from bitbucket to github.  This was done in recognition that more
  people interact with git/github than do with hg/bitbucket, as well as to
  follow our major dependency yt in making the same transition.  All previous
  repository history (e.g., commits, versions, tags, etc.) is retained under
  this transition. For users operating on the development branch of
  Trident, you must re-install Trident in order to continue to get updates.
  The :ref:`installation instructions <installation>` were updated accordingly.
- We totally rebuilt the testing interface to Trident, which includes
  more coverage in unit tests and answer tests over both grid-based and
  particle-based datasets.  We now have continuous integration through
  `Travis <https://travis-ci.org/trident-project/trident>`_ that tests the code
  daily and with each new pull request to assure consistent code results and to
  minimize bugs.  For more information, see :ref:`testing`.
- Much of the original Trident codebase was developed in yt as the base classes
  :class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
  and :class:`~trident.LightRay`.  We have now stripped these classes out of
  yt and moved them entirely into Trident for more flexibility, stability, and
  autonomy moving forward.  This should not affect the user as these changes
  were behind the scenes.
- Added ``store_observables`` keyword to
  :func:`~trident.SpectrumGenerator.make_spectrum` to store a
  dictionary of observable properties (e.g., tau, column density, and thermal_b)
  for each cell along a line of sight for use in post-processing.  See source
  of :class:`~trident.SpectrumGenerator` for more information.
- Added an approximate ``flux_error`` field to output spectra, since many
  observational tools require its presence.  See
  :func:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum.error_func`
  for more details.
- Made ``min_tau`` a keyword to
  :func:`~trident.SpectrumGenerator.make_spectrum` to enable higher precision
  (although more time intensive) absorption line deposition.
- Added ability to specify an arbitrary noise vector with
  :func:`~trident.SpectrumGenerator.add_noise_vector`.
- A `bugfix <https://github.com/yt-project/yt/pull/1611>`_ was made
  in yt to the temperature field for Gadget-based code outputs.  The internal
  energy field was mistakenly being read in co-moving instead of physical units,
  which led to gas temperatures being low by a factor of (1+z).
  This is now resolved in yt dev and thus we recommend Trident users use
  yt dev until yt 3.5 stable is released.
- `Another bugfix <https://github.com/astropy/astropy/pull/5782>`_ was made
  in Trident dependency `astropy <https://github.com/astropy/astropy/>`_ to
  the convolve function, which is used in
  :func:`~trident.SpectrumGenerator.apply_lsf`.  This may cause slight
  backwards-incompatible changes when applying line spread functions to
  post-process spectra.
- Replaced internal instances of ``particle_type`` with ``sampling_type`` to
  match similar yt conversion.

Version 1.0 (November 16, 2017)
-------------------------------

Initial release.  See our :ref:`method paper <citation>` for details.

- Create absorption-line spectra for any trajectory through a simulated
  data set mimicking both background quasar and down-the-barrel configurations.
- Reproduce the spectral characteristics of common instruments like the
  Cosmic Origins Spectrograph.
- Operate across the ultraviolet, optical, and infrared using customizable
  absorption-line lists.
- Trace simulated physical structures directly to spectral features.
- Approximate the presence of ion species absent from the simulation outputs.
- Generate column density maps for any ion.
- Provide support for all major astrophysical hydrodynamical codes.
