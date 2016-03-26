.. _installation:

Installation
============

The Trident source is most easily installed using pip.  But you can only
do this once you've installed yt and its dependencies.  If you already have
yt installed, then skip to step 2.


Step 1: Install yt  There are several methods for doing this all discussed in
detail in the yt docs here: http://yt-project.org/docs/dev/installing.html
Note that Trident currently requires the development branch of yt to work
properly.

If you're just starting out, we recommend the anaconda method as the least
work for new users to get yt installed as described here:
http://yt-project.org/docs/dev/installing.html#installing-yt-using-anaconda


Step 2: Install Trident using pip.::

    $ pip install trident


Step 3: Get the ``ion_balance`` datafile.  This datafile is needed to determine
the ionization fractions for ions from density, temperature, metallicity, and
redshift fields.  The first time you run Trident, it will attempt to
automatically download this file for you with some interactive prompts.

Alternatively, you can download it manually by downloading one of the
following files and either placing it in your $HOME/.trident directory, or
placing its path in the .trident/config file.  See the docs for more
information.

High-Res (1.6GB)::

    $ wget http://trident-project.org/data/ion_balance/hm2012_hr.h5

Low-Res (205MB)::

    $ wget http://trident-project.org/data/ion_balance/hm2012_lr.h5

The best results are achieved using the high-resolution data
file (1.6GB), but the low-resolution data is also acceptable (205MB).

You're now ready to use Trident.  
