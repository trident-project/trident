.. _installation:

Installation
============

Follow these steps to successfully install Trident and its dependencies.

.. _versions:

Versions of Trident
-------------------

Currently, there are two versions of Trident: the `stable version 
<http://trident.readthedocs.io/en/stable>`_, and the `developent version
<http://trident.readthedocs.io/en/latest>`_.  The stable version is tried,
tested, and stable.  The development version has new features under
current development, which may lead to the occasional error but is generally
very stable too.  The installation steps are slightly different between the two
(easier for stable than dev) so pay attention in the steps below.
Don't worry if you want to change later, you can always 
switch between the two versions easily enough by following the directions
in :ref:`uninstallation`.

.. _step-1:

Step 1: Install yt  
------------------

`yt <http://yt-project.org>`_ is a python-based software package for the 
analysis and visualization of a variety of different datasets, including 
astrophysical hydrodynamical data.  yt is a dependency of Trident, so you
must install it before Trident will work.  There are several methods for 
installing yt, which are all discussed in detail in the `yt installation 
documentation <http://yt-project.org/doc/installing.html>`_.  

Installing yt for the Stable Release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The stable release of Trident runs off of the latest stable release of yt.
You can install yt by following any of the methods on its `installation page
<http://yt-project.org/doc/installing.html>`_, but we find that the easiest
is to use the all-in-one install script::

    $ wget http://bitbucket.org/yt_analysis/yt/raw/stable/doc/install_script.sh
    $ bash install_script.py

If you already have conda installed, it's even faster and easier::

    $ conda install -c conda-forge yt

Finally, if you already have the development version of yt, you can skip the 
above steps and just update to the stable version::

    $ cd <path/to/yt/repo>
    $ hg up stable

Installing yt for the Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The development version of Trident runs off of the development version of yt,
since much of the functionality is co-developed between the two packages.
You can install yt by following any of the methods on its `installation page
<http://yt-project.org/docs/dev/installing.html>`_, but we find that the 
easiest is to use the all-in-one install script and build from source
(i.e. with settings ``INST_CONDA=1`` and ``INST_YT_SOURCE=1``)::

    $ wget http://bitbucket.org/yt_analysis/yt/raw/yt/doc/install_script.sh
    $ ... edit the install_script.sh to mark INST_YT_SOURCE=1 ...
    $ bash install_script.sh

Alternatively, if you already have conda installed, you can get the most
recent nightly build of the development version of yt with::

    $ conda install -c http://use.yt/with_conda yt

.. _install-trident:
.. _step-2:

Step 2: Install Trident
-----------------------

Installing the Stable Release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install the most recent stable release of Trident using pip::

    $ pip install trident

Installing the Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the development version, you'll pull the source code from its 
repository using mercurial, which should be installed as part of your yt 
installation.  After that, you'll use pip to install the source directly.  
Go to your desired source code installation directory and run::

    $ hg clone http://bitbucket.org/trident-project/trident
    $ cd trident
    $ pip install -e .

.. _step-3:

Step 3: Get Ionization Table and Verify Installation
----------------------------------------------------

In order to calculate the ionization fractions for various ions from 
density, temperature, metallicity fields, you will need an ionization table 
datafile and a configuration file.  Because this datafile can be large, it is
not packaged with the main source code.  The first time you try to do anything
that requires it, Trident will attempt to automatically set this all up for 
you with a series of interactive prompts.  **This step requires an internet 
connection the first time you run it.**

In addition, Trident provides a simple test function to verify that your 
install is functioning correctly.  This function not only tries to set up
your configuration and download your ion table file, but it will 
create a simple one-zone dataset, generate a ray through it, and 
create a spectrum from that ray.  This should execute very quickly, 
and if it succeeds it demonstrates that your installation has been totally 
successful::

    $ python
    >>> import trident
    >>> trident.verify()
    ...Series of Interactive Prompts...

If you cannot directly access the internet on this computer, or you lack write
access to your ``$HOME`` directory, or this step fails for any reason, please 
follow our documentation on :ref:`manual-config`.

.. _step-4:

Step 4: Science!
----------------

Congratulations, you're now ready to use Trident!  Please refer to the 
documentation for how to use it with your data or with one of our sample 
datasets.  Because Trident is in beta, the docs are not complete, and 
the API may still change in slight ways.  Please join our :ref:`mailing list 
<mailing-list>` for announcements about when the code is officially released.

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

Now, to confirm everything is working properly, verify your installation
following :ref:`step-3`.  If this fails or you have additional problems, 
please contact our mailing list.

.. _uninstallation:

Uninstallation or Switching Code Versions
-----------------------------------------

Uninstallation of the Trident source code is easy.  If you installed the 
stable version of the code via pip, just run::

    $ pip uninstall trident

If you installed the dev version of Trident, you'll have to delete the source
as well::

    $ pip uninstall trident
    $ rm -rf </path/to/trident/repo>

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
codebase.  

Updating to the Latest Stable Release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you installed the "stable" version of the code using pip, then 
you can easily update your trident and yt installations::

    $ pip install -U trident
    $ yt update

Updating to the Latest Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you installed the "development" version of the code, it's slightly more
involved::

    $ cd <path/to/trident/repo>
    $ hg pull
    $ hg up
    $ pip install -e .
    $ yt update

For more information on updating your yt installation, see the `yt update 
instructions 
<http://yt-project.org/docs/dev/installing.html#updating-yt-and-its-dependencies>`_.
