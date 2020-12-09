.. _installation:

Installation
============

Follow these steps to successfully install Trident and its dependencies.

.. _versions:

Versions of Trident
-------------------

There are currently two versions of Trident: a `stable version
<http://trident.readthedocs.io/en/stable>`_ and a `development version
<http://trident.readthedocs.io/en/latest>`_.  Make sure you are reading the
correct docs for the version you are using!

The stable version is tried
and tested and easy to install with pip.  The development version is actively
being updated with new features including superior support for particle-based
datasets (previously known as the demeshening).  Note that the stable version
of trident requires the stable version of yt, and the development version of
trident requires the development version of yt, due to some
backwards-incompatible changes regarding particle-support in yt/trident.

Thus, the installation steps are slightly different for stable and development,
so pay attention in the steps below.  Don't worry if you want to change later,
you can always switch between the two versions easily enough by following the
directions in :ref:`uninstallation`.

Trident's Major Dependency: yt
------------------------------

`yt <http://yt-project.org>`_ is a python-based software package for the
analysis and visualization of a different numerical datasets, including
astrophysical hydrodynamical data.  yt is the primary dependency of Trident,
so you must install it before Trident will work.  There are several methods
for installing yt, which are all discussed in detail in the `yt installation 
documentation <http://yt-project.org/doc/installing.html>`_.  Use the one
that is appropriate for you.  We find that using
`conda <https://docs.conda.io/en/latest/>`_ is the most streamlined and
reliable.

.. _stable-trident:

Installing the Stable Version of yt and Trident
-----------------------------------------------

Installation of the stable versions of yt and Trident is quite simple:

```
$ pip install yt
$ pip install trident
```

Now, you can try to run Trident for the first time, where it will download
some additional files.  See :ref:`step-3`, for more information.

```
$ python
>>> import trident
```

Follow the instructions to download the ion_balance table and then verify that
everything is working correctly.  You should now be ready to do some 
:ref:`step-4`

Installing the Development Version of yt and Trident
----------------------------------------------------

Step 0: Ensure Conda is Installed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conda is a package manager providing a clean, stand-alone installation of
python that is self-contained in its installation directory.  yt & trident
require a modern installation of python to work.  conda provides that 
installation.

You can see if conda is already installed by running:

```
$ conda -h
```

If conda is installed, move to the next step.  Otherwise install Mini-conda.

Use the appropriate conda install script for your architecture.  We recommend
getting the latest version of conda for Python3 for your architecture here:
https://repo.continuum.io/miniconda/

For modern macs:

```
$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
$ bash Miniconda3-latest-MacOSX-x86_64.sh
```

For modern linux machines:

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

At the end of the installation step, allow conda to add its installation directory to the $PATH.

.. _step-1:

Step 1: Install yt  
^^^^^^^^^^^^^^^^^^

First you need yt's major dependencies:

```
$ conda install numpy cython mpi4py git
```

Now you pull directly from the yt github repository to access
the up-to-date version of the source code and build it.

```
$ git clone https://github.com/yt-project/yt.git yt
$ cd yt
$ pip install -e .
$ cd ..
```

Note, you'll also need a separate library, yt_astro_analysis in order to 
have Trident work directly with yt.

```
$ git clone https://github.com/yt-project/yt_astro_analysis.git yt_astro_analysis
$ cd yt_astro_analysis
$ pip install -e .
$ cd ..
```

.. _install-trident:
.. _step-2:
.. _install-dev:

Step 2: Install Trident
^^^^^^^^^^^^^^^^^^^^^^^

Like yt, in order to get the development version of Trident, you must clone
and build the up-to-date source code from its repository.

```
$ git clone https://github.com/trident-project/trident.git trident
$ cd trident
$ pip install -e .
$ cd ..
```

.. _step-3:

Step 3: Get Ionization Table and Verify Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^

Congratulations, you're now ready to use Trident!  Please refer to the 
documentation for how to use it with your data or with one of our sample 
datasets.  A good place to start is the 
:ref:`annotated example <annotated-example>`, and the `example scripts found
in the source code 
<https://github.com/trident-project/trident/blob/master/examples>`_.

Please join our :ref:`mailing list 
<mailing-list>` or :ref:`slack channel <slack-channel>` for announcements
and updates when new features are added to the code.

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
    $ rm -rf <YOUR_PATH_TO_TRIDENT_REPO>

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

    $ cd <YOUR_PATH_TO_TRIDENT_REPO>
    $ git pull origin master
    $ pip install -e .
    $ yt update

For more information on updating your yt installation, see the `yt update 
instructions 
<http://yt-project.org/docs/dev/installing.html#updating-yt-and-its-dependencies>`_.
