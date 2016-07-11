.. _faq:

Frequently Asked Questions
==========================

.. _astropy-problem:

Trident fails with error: ``IndexError: list index out of range``
-----------------------------------------------------------------

One of our dependencies, `astropy <http://astropy.readthedocs.io/en/stable/>`_, 
has problems when you try to import it while you're in another package's
directory.  If you're in the Trident home directory when you run python
and import Trident, which imports astropy, you'll get this error::

    $ cd <path/to/trident>
    $ python
    >>> import trident
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "trident/__init__.py", line 42, in <module>
        from trident.spectrum_generator import \
    File "trident/spectrum_generator.py", line 18, in <module>
        from yt.analysis_modules.absorption_spectrum.api import \
    File "/Users/me/src/yt/yt/analysis_modules/absorption_spectrum/api.py", line 16, in <module>
        from .absorption_spectrum import \
    File "/Users/me/src/yt/yt/analysis_modules/absorption_spectrum/absorption_spectrum.py", line 36, in <module>
        pyfits = _astropy.pyfits
    File "/Users/me/src/yt/yt/utilities/on_demand_imports.py", line 38, in pyfits
        import astropy.io.fits as pyfits
    File "/Users/me/src/miniconda2/lib/python2.7/site-packages/astropy/__init__.py", line 274, in <module>
        __bibtex__ = _get_bibtex()
    File "/Users/me/src/miniconda2/lib/python2.7/site-packages/astropy/__init__.py", line 268, in _get_bibtex
        refcontents = re.findall(r'\{[^()]*\}', citation.read())[0]
    IndexError: list index out of range

To avoid this error, just move into another directory outside the Trident home 
directory.

.. _what-version-am-i-running:

What version of Trident am I running?
-------------------------------------

To learn what version of Trident you're running, type::

    $ python
    >>> import trident
    >>> print trident.__version__

If you have a version ending in dev, it means you're on the development branch
and you should also figure out which particular changeset you're running.  You
can do this by::

    $ cd <path/to/trident>
    $ hg id

To figure out what version of yt you're running, type::

    $ yt version

If you're writing to the mailing list with a problem, be sure to include all
of the above with your bug report or question.

Where is Trident installed?  Where are its data files?
------------------------------------------------------

One can easily identify where Trident is installed::

    $ python
    >>> import trident
    >>> print trident.path

The data files are located in that path with an appended ``/data``.

.. _mailing-list:

How do I join the mailing list?
-------------------------------

You can join our mailing list for announcements, bugs reports, and changes
at:

https://groups.google.com/forum/#!forum/trident-project-users
