.. _installation:

Installation
------------

The Trident source is most easily installed using pip.  But you can only
do this once you've installed yt and its dependencies.  If you already have
the dev version of yt installed, then skip to step 2.  


Step 1: Install yt  
^^^^^^^^^^^^^^^^^^

There are several methods for doing this all discussed in detail in the yt 
docs here: http://yt-project.org/docs/dev/installing.html .
Note that Trident currently requires the development branch of yt to work
properly.

If you're just starting out, we recommend the anaconda method as the least
work for new users to get yt installed as described here (be sure to get 
the development version of yt):
http://yt-project.org/docs/dev/installing.html#installing-yt-using-anaconda


Step 2: Install Trident
^^^^^^^^^^^^^^^^^^^^^^^

Now, you can install Trident using pip::

    $ pip install trident


Step 3: Get Ionization Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to calculate the ionization fractions for various ions from 
density, temperature, metallicity fields, you will need an ionization table 
datafile and a configuration file.  The first time you import Trident into 
Python, the code will attempt to automatically set this all up for you with 
a series of interactive prompts.  **This step requires an internet connection.**

.. code-block:: bash

    $ python
    >>> import trident

If you cannot directly access the internet on this computer, or you lack write
access to your ``$HOME`` directory, or this step fails for any reason, please 
follow our documentation on :ref:`manual-config`.

Step 4: Science!
^^^^^^^^^^^^^^^^

You're now ready to use Trident.  Please refer to the documentation for how
to use it with your data or with one of our sample datasets.


.. _manual-config:

Manually Installing your Ionization Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
