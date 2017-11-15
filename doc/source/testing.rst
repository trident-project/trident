.. _testing:

Testing
=======

We maintain a series of tests in Trident to make sure the code gives consistent
results and to catch accidental breakages in our source code and dependencies.
These tests are run by `Travis <https://travis-ci.org/>`_ automatically and 
regularly to assure consistency in functionality, but you can run them locally
too (see below).  The tests consist of a mix of unit tests (tests to assure Trident 
functions don't actively fail) and answer tests (tests comparing newly 
generated results against some old established results to assure consistency).

Running the Test Suite
----------------------

Running the test suite requires a version of Trident installed from
source (see :ref:`install-dev`).

The tests are run using the ``pytest`` Python module.  This can be
installed with ``conda`` or ``pip``.

.. code-block:: bash

   $ conda install pytest

The test suite requires a number of datasets for testing functionality.
Trident comes with a helper script that will download all the datasets and 
untar them.  Before running this, make sure you have the 
``answer_test_data_dir`` variable set in your config file (see :ref:`step-3`).  
This variable should point to a directory where these datasets will be stored.  
The helper script is located in the ``tests`` directory of the Trident source.

.. code-block:: bash

   $ cd tests
   $ python download_test_data.py

If this is your first time running the tests, then you need to generate a
"gold standard" for the answer tests. See :ref:`generating-answer-tests` 
before continuing with running the tests, otherwise your answer tests will 
fail.

Make sure you're on the desired version of yt and trident that you want to 
test and use (usually the tip of the development branch).  

.. code-block:: bash

   $ cd /path/to/yt/
   $ git checkout master
   $ pip install -e .
   $ cd /path/to/trident
   $ git checkout master
   $ pip install -e .

The test suite is run by calling ``py.test`` from within the ``tests`` 
directory.

.. code-block:: bash

   $ cd tests
   $ py.test
   ============================= test session starts ==============================
   platform darwin -- Python 3.6.0, pytest-3.0.7, py-1.4.32, pluggy-0.4.0
   rootdir: /Users/britton/Documents/work/yt/extensions/trident/trident, inifile:
   collected 52 items

   test_absorption_spectrum.py ..........
   test_download.py .
   test_generate.py .
   test_instrument.py .
   test_ion_balance.py ............
   test_light_ray.py .....
   test_line_database.py .......
   test_lsf.py ....
   test_pipelines.py ...
   test_plotting.py .
   test_ray_generator.py .
   test_spectrum_generator.py ......

   ========================= 52 passed in 117.32 seconds ==========================

If a test fails for some reason, you will be given a detailed traceback and
reason for it failing.  You can use this to identify what is wrong with your
source or perhaps a change in the code of your dependencies.

.. _generating-answer-tests:

Generating Test Results
-----------------------

These are a set of answer tests created with an older stable version of trident 
and yt that we think is accurate.  Periodically, this gold standard must be 
updated as bugs are caught or new more accurate behavior is enabled by new 
development.

To generate the test results, you must roll back the Trident and yt source back
to the older "trusted" versions of the code.  You can find the tags for the 
most recent "trusted" versions of the code by running 
``gold_standard_versions.py`` and the rebuilding yt and Trident with these
versions of the code
Lastly, set the ``TRIDENT_GENERATE_TEST_RESULTS`` environment variable to 1 
and run the tests:

.. code-block:: bash

   $ cd tests
   $ python gold_standard_versions.py
   
   Latest Gold Standard Commit Tags
   yt = 38b79c094ca9
   Trident = gold-standard-v1

   To update to them, `git checkout <tag>` in appropriate repository

   $ cd /path/to/yt
   $ git checkout 38b79c094ca9
   $ pip install -e .
   $ cd /path/to/trident
   $ git checkout gold-standard-v1
   $ pip install -e .
   $ export TRIDENT_GENERATE_TEST_RESULTS=1
   $ py.test

The test results should now be stored in the ``answer_test_data_dir`` that
you specified in your Trident configuration file (see :ref:`step-3`).
