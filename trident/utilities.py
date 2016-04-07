"""
Miscellaneous Utilities for Trident

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import urllib2
import gzip
import os
import sys
from os.path import \
    expanduser
from ConfigParser import \
    SafeConfigParser
import requests
import tempfile
import shutil
from yt.funcs import \
    get_pbar

def ensure_directory(directory):
    """
    Ensures a directory exists by creating it if it does not.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def download_file(url, progress_bar=True, local_directory=None, 
                 local_filename=None):
    """
    Downloads a file from the provided URL.  If progress_bar is set
    to true, generates a progress bar for the user (Default: True).  
    Optionally accepts a local directory and local filename in which to 
    place the file.  By default, it will download the file into the current 
    directory with its original filename.
    """

    # Following the base description on stack overflow:
    # http://stackoverflow.com/questions/22676/how-do-i-download-a-file-over-http-using-python

    # Set defaults
    if local_filename is None:
        local_filename = url.split('/')[-1]
    if local_directory is None:
        local_directory = '.'
    ensure_directory(local_directory)

    # open local file handle
    filepath = os.path.join(local_directory, local_filename)
    filehandle = open(filepath, 'wb')

    # Get information about remote filesize
    u = urllib2.urlopen(url)
    meta = u.info()
    filesize = int(meta.getheaders("Content-Length")[0])/2**10 # in kB
    if progress_bar:
        pbar = get_pbar("Downloading file: %s" % local_filename, filesize)
    filesize_dl = 0
    block_sz = 8192

    # Download file in blocks and update statusbar until done with transfer
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        filesize_dl += int(len(buffer))/2**10
        filehandle.write(buffer)
        if progress_bar:
            pbar.update(filesize_dl)
    if progress_bar:
        pbar.finish()
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

def parse_config():
    """
    This function runs every time Trident gets imported.  It assures that
    Trident knows where to look for ion table datafiles.  If something is
    missing, either the config file, or a valid ion table datafile, it
    tries to set everything up for the user.  For more information about doing
    this manually see the installation documentation.
    """
    # Assure the ~/.trident directory exists, and read in the config file.
    home = expanduser("~")
    directory = os.path.join(home, '.trident')
    config_filename = os.path.join(directory, 'config')
    try:
        parser = SafeConfigParser()
        parser.read(config_filename)
        ion_table_dir = parser.get('Trident', 'ion_table_dir')
        ion_table_file = parser.get('Trident', 'ion_table_file')
    except:
        ion_table_dir, ion_table_file = create_config()

    if not os.path.exists(os.path.join(ion_table_dir, 
                                       ion_table_file)):
        print ""
        print "No ion table data file found in %s" % ion_table_dir
        ion_table_file = get_datafiles(ion_table_dir)
        parser.set('Trident', 'ion_table_file', ion_table_file)
        with open(config_filename, 'w') as configfile:
            parser.write(configfile)
    return ion_table_dir, ion_table_file

def create_config():
    """
    This function is run if it appears that the configuration has not yet
    been set up for the user.  It will attempt to create a configuration file
    and download an ion table datafile from the web.  It does this using
    user interaction from the python prompt.
    """
    default_dir = expanduser('~/.trident')
    trident()
    print "It appears that this is your first time using Trident.  To finalize your"
    print "Trident installation, you must:"
    print " * create a `~/.trident` directory"
    print " * create a config file in your `~.trident` directory"
    print " * download an ion table file for calculating ionization fractions"
    print ""
    print "You can do this manually by following the installation docs, or we can"
    print "do it automatically now if you have web access."
    print ""
    print "Would you like to do this automatically? ([y]/n)"
    value = raw_input().rstrip()
    if not value == '' and not value == 'y':
        sys.exit('Instructions at http://trident.readthedocs.org/en/latest/Installation.html')

    print ""
    print "Where would you like Trident to store the ion table file?"
    print "[%s]" % default_dir

    # First assure that the .trident directory is created for storing
    # the config file.
    ensure_directory(default_dir)

    datadir = raw_input().rstrip()
    if datadir == '':
        datadir = default_dir
    datadir = expanduser(datadir)

    # Try to create data directory if it doesn't exist
    try:
        ensure_directory(datadir)
        print "Using %s" % datadir
        print ""
    except:
        print "Cannot create directory %s" % datadir
        raise

    # Get the datafile from the web
    datafile = get_datafiles(datadir=datadir)

    # Create the config file and make it to the datadir and datafiles chosen
    config = SafeConfigParser()
    config.add_section('Trident')
    config.set('Trident', 'ion_table_dir', datadir)
    config.set('Trident', 'ion_table_file', datafile)
    config_filename = expanduser('~/.trident/config')
    with open(config_filename, 'w') as configfile:
        config.write(configfile)

    print ""
    print "Installation complete.  Now let's do some science!"

    # Return the ion_table_dir and ion_table_file so they
    # can be set as trident global variables for future use by ion_balance 
    # classes/functions
    return datadir, datafile

def get_datafiles(datadir=None, url=None):
    """
    If the user lacks an ion table datafile, this attempts to download one
    from the web using interaction with the user.
    """
    if datadir is None:
        datadir = expanduser('~/.trident')
    ensure_directory(datadir)

    # ion table datafiles are stored here
    if url is None:
        url = 'http://trident-project.org/data/ion_table/'
    
    # Try to figure out which datafiles are available on the remote server
    try:
        page = str(requests.get(url).text)
    except:
        print "Cannot seem to access %s; Do you have internet access?" % url
        raise

    line_list = page.split('\n')
    i = 1 # counter
    filenames = []
    print "The following data files are available:"
    for line in line_list:
       if not line.startswith('<!--X'): # identifies files by comments
            continue
       line = line[6:-3]  # strips off the HTML comment chars
       filenames.append(line.split()[0])
       print "%d) %s" % (i, line)
       i += 1

    # User chooses which file
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

    # Actually downloads the file to a temporary directory before gunzipping.
    fileurl = os.path.join(url, filename)
    tempdir = tempfile.mkdtemp()
    print ""
    download_file(fileurl, local_directory=tempdir)
    print "  Unzipping file: %s" % filename
    gunzip_file(os.path.join(tempdir, filename), 
                out_filename=os.path.join(datadir, filename[:-3]))
    shutil.rmtree(tempdir) 
    return filename[:-3]

def trident():
    """
    Prints a nice ASCII logo!
    """
    print """
MMMMMMMMMMMMMMMMMMM.............................................................
M...B....C....D...M.MMMMMMMMMM.......MM........MM.....................MM7.......
M...MM...M...MM...M.....MM.....................MM.....................MM7.......
M....MM..M..MM....M.....MM...MMMMMM..MM..MMMMMMMM..MMMMMMM..MMDMMMMM MMMMMMM....
M....MM..M..MM....M.....MM...MM ..MM.MM. MM....MM.,MM...MM+.MM....MM..MM7.......
M.....MMMMMMM.....M.....MM...MM......MM. MM....MM.8MMMMMMMD.MM....MM..MM7.......
M....... M........M.....MM...MM......MM..MM...8MM. MM...._..MM....MM..MMZ.MM ...
M........M........M.....MM...MM......MM..MMMMM.MM. MMMMMMM..MM....MM...MMMMM....
MMMMMMMMMMMMMMMMMMM.............................................................
"""

def trident_path():
    """
    A function returning the path of the trident installation directory.
    Useful for identifying where data files are (e.g. path/data).  Note that
    ion table datafiles are downloaded separate and placed in another
    location according to the ~/.trident/config file.
    """
    return '/'.join(os.path.dirname(__file__).split('/')[:-1])
