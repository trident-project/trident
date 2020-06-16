"""
Tests for Utilities Code

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
import tempfile
import filecmp
from trident.utilities import \
    ensure_directory, \
    gunzip_file, \
    gzip_file, \
    make_onezone_dataset, \
    make_onezone_ray

def test_ensure_directory():
    """
    Tests ensure_directory code, which ensures a path exists and if not, 
    it creates such a path
    """
    tempdir = tempfile.mkdtemp()
    ensure_directory(tempdir)
    assert os.path.exists(tempdir)

    # is new random path created?
    subdir = os.path.join(tempdir, 'random')
    assert not os.path.exists(subdir)
    ensure_directory(subdir)
    assert os.path.exists(subdir)

def test_gzip_unzip_file():
    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()

    # Create a temporary file
    fp = tempfile.NamedTemporaryFile()
    fp.write(b'Hello world!')
    filepath = fp.name
    directory, filename = os.path.split(filepath)
    gzip_filename = os.path.join(tmpdir, filename+'.gz')
    gunzip_filename = os.path.join(tmpdir, filename)

    # Gzip it
    gzip_file(filepath, out_filename=gzip_filename, cleanup=False)

    # Gunzip it
    gunzip_file(gzip_filename, out_filename=gunzip_filename, cleanup=False)

    # Ensure contents of gunzipped file are the same as original
    filecmp.cmp(filepath, gunzip_filename)

