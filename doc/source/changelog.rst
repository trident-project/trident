.. _changelog:

Changelog
=========

This document summarizes changes to the codebase with time for Trident.

Contributors
------------

The `CREDITS file <https://github.com/trident-project/trident/blob/master/CREDITS>`_ 
has an updated list of contributors to the codebase.

Version 1.1
-----------

- Trident development has changed from mercurial to git, and the source has 
  moved from bitbucket to github.  This was done in recognition that more
  people interact with git/github than do with hg/bitbucket, as well as to
  follow our major dependency yt in making the same transition.  All previous
  repository history (e.g., commits, versions, tags, etc.) is retained under
  this transition, but for users operating on the development branch of 
  Trident, you must re-install Trident in order to continue to get updates.
- We totally rebuilt the testing interface to Trident, which includes 
  more coverage in unit tests and answer tests over both grid-based and
  particle-based datasets.  We now have continuous integration through Travis
  that tests the code daily and with each new pull request to assure 
  consistent code results and to minimize bugs.  For more information, see
  :ref:`testing`.
- Much of the original Trident codebase was developed in yt as base the classes
  `AbsorptionSpectrum` and `LightRay`.  We have now stripped these classes 
  out of yt and moved them entirely into Trident for more flexibility, 
  stability, and autonomy moving forward.  This should not affect the user
  as these changes were behind the scenes.
- Added ability to specify an arbitrary noise vector with `SpectrumGenerator.add_noise_vector`

Version 1.0
-----------

Initial release.  See our :ref:`method paper <citation>` for details. 
- Create absorption-line spectra for any trajectory through a simulated 
  data set mimicking both background quasar and down-the-barrel configurations
- Reproduce the spectral characteristics of common instruments like the 
  Cosmic Origins Spectrograph
- Operate across the ultraviolet, optical, and infrared using customizable 
  absorption-line lists
- Trace simulated physical structures directly to spectral features
- Approximate the presence of ion species absent from the simulation outputs
- Generate column density maps for any ion
- Provide support for all major astrophysical hydrodynamical codes
