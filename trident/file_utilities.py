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
import sys
from os.path import expanduser
from ConfigParser import SafeConfigParser
import requests
import tempfile
import shutil

def ensure_directory(directory):
    """
    Ensures a directory exists by creating it if it does not.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def download_file(url, local_directory=None, local_filename=None):
    """
    Downloads a file from the provided URL.  If progress_bar is set
    to true, generates a progress bar for the user (Default: True).  
    Optionally accepts a local directory and local filename in which to 
    place the file.  By default, it will place the using in the current 
    directory with its original filename.
    """

    # Following the base description on stack overflow:
    # http://stackoverflow.com/questions/22676/how-do-i-download-a-file-over-http-using-python

    if local_filename is None:
        local_filename = url.split('/')[-1]
    if local_directory is None:
        local_directory = '.'
    ensure_directory(local_directory)
    filepath = os.path.join(local_directory, local_filename)
    filehandle = open(filepath, 'wb')

    u = urllib2.urlopen(url)
    meta = u.info()
    filesize = int(meta.getheaders("Content-Length")[0])
    print "Downloading file: %s     Size [bytes]: %s" % (local_filename, filesize)
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

#def _scan_remote_directory():

def parse_config():
    """
    This function runs every time Trident gets imported.  It assures that
    Trident knows where to look for ion_balance datafiles.  If something is
    missing, either the config file, or a valid ion_balance datafile, it
    tries to set it up for the user.
    """
    # Assure the ~/.trident directory exists, and read in the config file.
    home = expanduser("~")
    directory = os.path.join(home, '.trident')
    config_filename = os.path.join(directory, 'config')
    try:
        parser = SafeConfigParser()
        parser.read(config_filename)
        ion_balance_data_dir = parser.get('Trident', 'ion_balance_data_dir')
        ion_balance_data_file = parser.get('Trident', 'ion_balance_data_file')
    except:
        ion_balance_data_dir, ion_balance_data_file = create_config()
    return ion_balance_data_dir, ion_balance_data_file

def create_config():
    """
    """
    ensure_directory(expanduser('~/.trident'))
    datadir, datafile = get_datafiles()
    config = SafeConfigParser()
    config.add_section('Trident')
    config.set('Trident', 'ion_balance_data_dir', datadir)
    config.set('Trident', 'ion_balance_data_file', datafile)
    config_filename = expanduser('~/.trident/config')
    with open(config_filename, 'w') as configfile:
        config.write(configfile)
    return datadir, datafile

def get_datafiles(url=None):
    """
    """
    default_dir = expanduser('~/.trident')
    print ""
    print "To finalize your Trident installation, we will"
    print "  create a `~/.trident` directory"
    print "  create a config file in your `~.trident` directory"
    print "  download a data file for calculating ionization fractions"
    print ""
    print "Where would you like Trident to store the data file?"
    print "[%s]" % default_dir

    # First assure that the .trident directory is created for storing
    # the config file.
    ensure_directory(default_dir)

    datadir = raw_input().rstrip()
    if datadir == '':
        datadir = default_dir
    datadir = expanduser(datadir)
    ensure_directory(datadir)

    if url is None:
        url = 'http://trident-project.org/data/ion_balance/'
    page = str(requests.get(url).text)
    line_list = page.split('\n')
    i = 1
    filenames = []
    print "The following data files are available:"
    for line in line_list:
       if not line.startswith('<!--X'):
            continue
       line = line[6:-3]  # strips off the HTML comment chars
       filenames.append(line.split()[0])
       print "%d) %s" % (i, line)
       i += 1

    print ""
    print "Which number would you like to download and use? [1]"
    value = raw_input().rstrip()
    if value == '':
        value = '1'
    while 1:
        try:
            filename = filenames[int(value)-1] 
            break
        except IndexError:
            print "%d is not a valid option.  Please pick one of the listed numbers."

    fileurl = os.path.join(url, filename)
    tempdir = tempfile.mkdtemp()
    print ""
    download_file(fileurl, local_directory=tempdir)
    print " Unzipping file: %s" % filename
    gunzip_file(os.path.join(tempdir, filename), 
                out_filename=os.path.join(datadir, filename[:-3]))
    shutil.rmtree(tempdir)
    return datadir, filename[:-3]
