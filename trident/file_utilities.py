"""
Utilities for downloading and storing ion_balance data files.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import urllib2
import gzip
import os

def download(url, local_directory=None, local_filename=None):
    """
    Downloads a file from the provided URL.  If progress_bar is set
    to true, generates a progress bar for the user (Default: True).  
    Optionally accepts a local directory and local filename in which to 
    place the file.  By default, it will place the using in the current 
    directory with its original filename.
    """

    # Following the base description on stack overflow:
    # http://stackoverflow.com/questions/22676/how-do-i-download-a-file-over-http-using-python

    filename = url.split('/')[-1]
    u = urllib2.urlopen(url)
    filehandle = open(filename, 'wb')
    meta = u.info()
    filesize = int(meta.getheaders("Content-Length")[0])
    print "Downloading: %s Bytes: %s" % (filename, filesize)

    filesize_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        filesize_dl += len(buffer)
        filehandle.write(buffer)
        status = r"%10d  [%3.2f%%]" % (filesize_dl, filesize_dl * 100. / filesize)
        status = status + chr(8)*(len(status)+1)
        print status,

    filehandle.close()

def gunzip_file(in_filename, out_filename=None, cleanup=True):
    """
    Uses gzip to unzip a file and save it to disk.  If out_filename is not
    provided, it simply saves the gunzipped file without the ".gz" suffix.
    if cleanup is True, removes the zipped version of the file after 
    gunzipping.
    """
    in_file = gzip.open(in_filename, 'rb')
    # if out_filename not set, defaults to stripping off the .gz of in_filename
    if out_filename is None: 
        out_filename = ".".join(in_filename.split('.')[:-1])
    out_file = open(out_filename, 'wb')
    out_file.write( in_file.read() )
    in_file.close()
    out_file.close()
    if cleanup:
        os.remove(in_filename)

def gzip_file(in_filename, out_filename=None, cleanup=True):
    """
    Uses gzip to zip a file and save it to disk.  If out_filename is not
    provided, it simply saves the gzipped file with a ".gz" suffix.  If 
    cleanup is True, removes the unzipped version of the file after gzipping.
    """
    in_file = open(in_filename, 'rb')
    # if out_filename not set, defaults to appending .gz to the in_filename
    if out_filename is None: 
        out_filename = in_filename + ".gz"
    out_file = gzip.open(out_filename, 'wb')
    out_file.write( in_file.read() )
    in_file.close()
    out_file.close()
    if cleanup:
        os.remove(in_filename)
