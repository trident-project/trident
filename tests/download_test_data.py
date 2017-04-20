"""
Testing utilities for Trident

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
import tarfile
from trident.testing import \
    answer_test_data_dir
from trident.utilities import \
    download_file

def download_datasets(local_dir=None, progress_bar=True):
    if local_dir is None:
        local_dir = answer_test_data_dir
    urls = open("test_datasets.txt", "r").readlines()
    for url in urls:
        if url.strip().startswith("#"):
            continue
        url = url.strip()
        print ("Downloading %s to %s." % (url, local_dir))
        filename = os.path.join(local_dir, os.path.basename(url))
        download_file(url, local_directory=local_dir, progress_bar=True)
        assert os.path.exists(filename), \
          "Failed to download %s." % url
        print ("Untarring %s." % filename)
        tar = tarfile.open(filename)
        tar.extractall(path=local_dir)
        tar.close()
    print ("Data downloaded and untarred successfully. Tarfiles can be deleted")

if __name__ == "__main__":
    download_datasets()
