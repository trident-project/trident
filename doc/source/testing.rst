.. _testing:

Running the Test Suite
======================

Running the test suite requires a version of Trident installed from
source (see :ref:`install-dev`).

Trident maintains a suite of unit and answer tests to ensure that
development doesn't change code behavior in unexpected ways.  The
tests are run using the ``pytest`` Python module.  This can be
installed with ``conda`` or ``pip``.

.. code-block:: bash

   $ conda install pytest

The test suite requires a number of datasets as well as results
files for answer comparison.  Trident comes with a helper script
that will download all the datasets and untar them.  Before running
this, make sure you have the ``answer_test_data_dir`` variable set in
your config file (see :ref:`step-3`).  This variable should point to
a directory.  The helper script is located in the ``tests`` directory
of the Trident source.

.. code-block:: bash

   $ cd tests
   $ python download_test_data.py

Once the test data has been downloaded, the test suite is run by
calling ``py.test`` from within the ``tests`` directory.

.. code-block:: bash

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

Generating Test Results
=======================

If new tests have been added or the code's behavior has changed (in a good way)
such that the tests no longer pass, new results must be generated.  Before
generating new results, be sure to update the results version number in
``tests/test_results_version.txt``.  Then, set the ``TRIDENT_GENERATE_TEST_RESULTS``
environment variable to 1 and rerun the tests:

.. code-block:: bash

   $ cd tests
   $ export TRIDENT_GENERATE_TEST_RESULTS=1
   $ py.test

The results must then be tarred up and uploaded to the Trident website.
