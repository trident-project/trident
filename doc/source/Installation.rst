.. _installation:

Installation
============

To install Trident, please assure yt is installed (right now it only works on
the tip of the yt development branch, not with the stable version).  Generally,
the easiest way to install yt is using Anaconda as described here:
http://yt-project.org/docs/dev/installing.html#installing-yt-using-anaconda
Make sure that you also obtain the source code as instructed.
For notes on all other installation methods see:
http://yt-project.org/docs/dev/installing.html.

Once you have yt working, download the Trident repository at:
http://bitbucket.org/trident-project/trident .  With mercurial installed, you
can do this by:

.. highlight:: none

::

   $ hg clone https://bitbucket.org/trident-project/trident
   $ cd trident
   $ python setup.py develop

In addition, you will need a data file for running Trident.  Get
one of the following files and put it in the trident/data/ion_balance/
directory.  The best results are achieved using the high-resolution data
file (1.6GB), but the low-resolution data is also acceptable (205MB).

To make sure the data tables end up in the right place, go into the
ion_balance data directory:

.. highlight:: none

::

   $ cd <path/to/trident/repo>/data/ion_balance

Then, if you want the high resolution version (1.6GB):

.. highlight:: none

::

   $ wget http://trident-project.org/data/ion_balance/hm2012_hr.h5

Or, if you'd prefer to use the low resolution version (205MB):

.. highlight:: none

::

   $ wget http://trident-project.org/data/ion_balance/hm2012_lr.h5


Now you should be ready to go.  To test your installation, you can run our
working script.  You will require the enzo_cosmology_plus yt public dataset
in order for this example to work.  To get it and test it, exit the trident
source directory into a test working directory and run the following:

.. highlight:: none

::

    $ cd <sandbox/directory>
    $ wget http://yt-project.org/data/enzo_cosmology_plus.tar.gz
    $ tar -zxvf enzo_cosmology_plus.tar.gz
    $ cp <path/to/trident/repo>/examples/working_script.py .
    $ python working_script.py

Now it's time for :ref:`making-spectra`!
