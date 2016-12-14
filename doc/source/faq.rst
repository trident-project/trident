.. _faq:

Frequently Asked Questions
==========================

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

How do I learn more about the algorithms used in Trident?
---------------------------------------------------------

We have a full description of all the methods used in Trident including
citations to previous related works in our `Trident method paper 
<https://arxiv.org/abs/1612.03935>`_.

How do I cite Trident in my research?
-------------------------------------

Please cite our `method paper <https://arxiv.org/abs/1612.03935>`_ and also 
include a footnote to our `website <http://trident-project.org>`_ so others 
know where to acquire the Trident software.


