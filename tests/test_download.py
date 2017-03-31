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

def test_download_file():
    env = os.environ
    if int(env.get("RUN_DOWNLOAD_TEST", 0)) == 1:
        download_datasets()
