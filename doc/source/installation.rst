.. _installation:

Installation
============

Follow these steps to successfully install Trident and its dependencies.

Step 1: Install yt  
------------------

`yt <http://yt-project.org>`_ is a python-based software package for the 
analysis and visualization of a variety of different datasets, including 
astrophysical hydrodynamical data.  Installation of Trident currently 
requires that you install the development version of yt as described below.  
If you already have a working installation of the development branch of yt, 
then skip to step 2.

There are several methods for installing yt, which are all discussed in 
detail in the `yt installation documentation 
<http://yt-project.org/docs/dev/installing.html>`_.  Note that Trident 
currently requires the development branch of yt to work properly.
If you're just starting out, we recommend the `anaconda installation method 
<http://yt-project.org/docs/dev/installing.html#installing-yt-using-anaconda>`_ 
as the least work for new users to get yt working.

.. _install-trident:

Step 2: Install Trident
-----------------------

There are two ways to install the Trident code itself.  The easiest 
method is to use pip to install the most recent stable version of Trident.  
Alternatively, you can install the development version of Trident--the version 
appropriate for users who want to hack on the code and get access 
to features not yet available in the stable release.  Don't worry, you can 
always switch between the two versions easily enough by following the directions
in :ref:`uninstallation`.

Installing the Stable Release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install the most recent stable release of Trident using pip::

    $ pip install trident

Installing the Development Branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the development version, you'll pull the source code from its 
repository using mercurial, which should be installed as part of your yt 
installation.  After that, you'll use pip to install the source directly.  
Go to your desired source code installation directory and run::

    $ hg clone http://bitbucket.org/trident-project/trident
    $ cd trident
    $ pip install -e .

Step 3: Get Ionization Table
----------------------------

In order to calculate the ionization fractions for various ions from 
density, temperature, metallicity fields, you will need an ionization table 
datafile and a configuration file.  Because this datafile can be large, it is
not packaged with the main source code.  The first time you import Trident into 
Python, the code will attempt to automatically set this all up for you with 
a series of interactive prompts.  **This step requires an internet connection.**

.. code-block:: bash

    $ python
    >>> import trident
    ...Series of Interactive Prompts...

If you cannot directly access the internet on this computer, or you lack write
access to your ``$HOME`` directory, or this step fails for any reason, please 
follow our documentation on :ref:`manual-config`.

Step 4: Verify Installation
---------------------------

Once you've installed yt, Trident, and downloaded your ion table data file, 
everything *should* be working, but it's good to verify your installation.
Trident provides a simple test function to verify that your install is 
able to create a simple one-zone dataset, generate a ray through it, and 
create a spectrum from that ray.  This should take about 5 seconds to run 
on a modern computer, and if it succeeds it demonstrates that your installation
has been totally successful::

    $ python
    >>> import trident
    >>> trident.verify()

Step 5: Science!
----------------

Congratulations, you're now ready to use Trident!  Please refer to the 
documentation for how to use it with your data or with one of our sample 
datasets.


.. _manual-config:

Manually Installing your Ionization Table
-----------------------------------------

If for some reason you are unable to install the config file and ionization
table data automatically, you must set it up manually.  When Trident runs,
it looks for a configuration file called ``config.tri`` in the 
``$HOME/.trident`` directory or alternatively in the current working 
directory (for users lacking write access to their ``$HOME`` directories).  
This configuration file is simple in that it tells Trident a few things about 
your install including the location and filename of your desired ionization 
table.  Manually create a text file called ``config.tri`` with contents 
following the form::

    [Trident]
    ion_table_dir = ~/.trident
    ion_table_file = hm2012_hr.h5

To manually obtain an ion table datafile, download and gunzip one from:
http://trident-project.org/data/ion_table .  While the ``config.tri`` file needs 
to exist in your ``$HOME/.trident`` directory or in the working directory
when you import trident, the ion_table datafile can exist anywhere on the 
file system.  Just assure that the config file points to the proper location 
and filename of the ion table datafile.

Now, to confirm everything is working properly, try to load Trident in a python 
session::

    $ python
    >>> import trident

If this fails or you have additional problems, please contact our mailing list.

.. _uninstallation:

Uninstallation or Switching Code Versions
-----------------------------------------

Uninstallation of the Trident source code is easy.  If you installed the 
stable version of the code via pip, just run::

    $ pip uninstall trident

If you installed the dev version of Trident, you'll have to delete the source
as well::

    $ pip uninstall trident
    $ rm -rf /path/to/trident/source

If you want to switch between the two stable and development versions, just
*uninstall* your version of the code as above, and then install the desired
version as described in :ref:`install-trident`

To fully remove the code from your system, remember to remove any ion table
datafiles you may have downloaded in your ``$HOME/.trident`` directory, 
and follow the instructions for how to `uninstall yt 
<http://yt-project.org/docs/dev/installing.html>`_.

.. _updating:

Updating to the Latest Version
------------------------------

If you want more recent features, you should periodically update your Trident
codebase.  If you installed the "stable" version of the code using pip, then 
simply run::

    $ pip install -U trident

If you installed the "development" version of the code, it's slightly more
involved::

    $ cd <path/to/trident>
    $ hg pull
    $ hg up
    $ pip install -e .
