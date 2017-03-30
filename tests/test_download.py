"""
Unit test for download_file function
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os

from download_test_data import \
    download_datasets
from trident.utilities import \
    download_file

def test_download_file():
    env = os.environ
    if int(env.get("RUN_DOWNLOAD_TEST", 0)) == 1:
        local_dir = "."
        urls = ["http://trident-project.org/data/ion_table/config.tri",
                "http://trident-project.org/data/ion_table/hm2012_lr.h5.gz"]
        for url in urls:
            download_file(url, local_directory=local_dir)
            assert os.path.exists(
                os.path.join(local_dir, os.path.basename(url)))
        download_datasets(local_dir=local_dir)
