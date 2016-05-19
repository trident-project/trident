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
import h5py as h5
import numpy as np
from yt.units import \
    cm, \
    pc, \
    Zsun
from yt import \
    load_uniform_grid, \
    YTArray, \
    load

def ensure_directory(directory):
    """
    Ensures a directory exists by creating it if it does not.  Works for 
    both absolute and relative paths.

    **Parameters**

    :directory: string

        Directory you wish to ensure exists.

    **Example**

    >>> ensure_directory("~/test")
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def download_file(url, progress_bar=True, local_directory=None, 
                 local_filename=None):
    """
    Downloads a file from the provided URL.  
    
    **Parameters**

    :url: string

        The web address of the file to download.
        
    :progress_bar: boolean, optional

        Will generate a progress bar for the user as the file downloads. 
        iPython/Jupyter friendly.
        Default: True

    :local_directory: string, optional

        Absolute or relative path of a local directory where the file 
        will be downloaded.  If set to None, will default to current
        working directory.
        Default: None

    :local_filename: string, optional

        Local filename where the file will be downloaded.  If set to None, 
        will default to filename of downloaded file.
        Default: None
    
    **Example**
    
    >>> download_file("http://trident-project.org/data/ion_table/hm2012_lr.h5.gz")
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
    Uncompress a file using gunzip.
    
    **Parameters**

    :in_filename: string

        The filename of the gzipped file.

    :out_filename: string, optional

        The filename where the unzipped file will be saved.  If set to None,
        this will be the ``in_filename`` without the ``.gz`` suffix.
        Default: None

    :cleanup: boolean, optional
        
        Remove the zipped version of the file after unzipping it?
        Default: True.
    
    **Example**
    
    >>> gunzip_file("hm2012_lr.h5.gz")
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
    Compress a file using gzip.
    
    **Parameters**

    :in_filename: string

        The filename of the file to be gzipped.

    :out_filename: string, optional

        The filename where the zipped file will be saved.  If set to None,
        this will be the ``in_filename`` with a ``.gz`` suffix.
        Default: None

    :cleanup: boolean, optional
        
        Remove the unzipped version of the file after zipping it?
        Default: True.
    
    **Example**
    
    >>> gzip_file("hm2012_lr.h5")
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
    Parse the Trident local configuration file.  This function is called
    whenever Trident is imported, and it assures that Trident knows where
    the default ion table datafiles exist.  If a ``config.tri`` file doesn't 
    exist in ``$HOME/.trident`` or in the current working directory, then
    Trident will launch the :class:`~trident.create_config` function to
    try to automatically generate one for the user.  For more information
    on this process, see the installation documentation.
    """
    # Assure the ~/.trident directory exists, and read in the config file.
    home = expanduser("~")
    directory = os.path.join(home, '.trident')
    config_filename = os.path.join(directory, 'config.tri')

    # If config file exists in current directory, use it instead of file in
    # $HOME/.trident.  Stopgap for situations where user cannot access $HOME
    local_filename = os.path.join(os.getcwd(), 'config.tri')
    if os.path.exists(local_filename):
        config_filename = local_filename
    try:
        parser = SafeConfigParser()
        parser.read(config_filename)
        ion_table_dir = parser.get('Trident', 'ion_table_dir')
        ion_table_file = parser.get('Trident', 'ion_table_file')
    except:
        ion_table_dir, ion_table_file = create_config()

    if not os.path.exists(os.path.join(ion_table_dir, 
                                       ion_table_file)):
        print("")
        print("No ion table data file found in %s" % ion_table_dir)
        ion_table_file = get_datafiles(ion_table_dir)
        parser.set('Trident', 'ion_table_file', ion_table_file)
        with open(config_filename, 'w') as configfile:
            parser.write(configfile)
    return ion_table_dir, ion_table_file

def create_config():
    """
    Create a Trident configuration file using interaction with the user.
    This function is called by :class:`~trident.parse_config` if it appears 
    that the configuration has not yet been set up for the user.  It will 
    attempt to create a configuration file and download an ion table 
    datafile from the web.  It does this using user interaction from the 
    python prompt.
    """
    default_dir = expanduser('~/.trident')
    trident()
    print("It appears that this is your first time using Trident.  To finalize your")
    print("Trident installation, you must:")
    print(" * create a `~/.trident` directory")
    print(" * create a config.tri file in your `~/.trident` directory")
    print(" * download an ion table file for calculating ionization fractions")
    print("")
    print("You can do this manually by following the installation docs, or we can")
    print("do it automatically now if you have web access.")
    print("")
    print("Would you like to do this automatically? ([y]/n)")
    value = raw_input().rstrip()
    if not value == '' and not value == 'y':
        sys.exit('Instructions at http://trident.readthedocs.org/en/latest/installation.html')

    print("")
    print("Where would you like Trident to store the ion table file?")
    print("[%s]" % default_dir)

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
        print("Using %s" % datadir)
        print("")
    except:
        print("Cannot create directory %s" % datadir)
        raise

    # Get the datafile from the web
    datafile = get_datafiles(datadir=datadir)

    # Create the config file and make it to the datadir and datafiles chosen
    config = SafeConfigParser()
    config.add_section('Trident')
    config.set('Trident', 'ion_table_dir', datadir)
    config.set('Trident', 'ion_table_file', datafile)
    config_filename = expanduser('~/.trident/config.tri')
    with open(config_filename, 'w') as configfile:
        config.write(configfile)

    print("")
    print("Installation complete.  I recommend verifying your installation")
    print("to assure that everything is working.  Try: trident.verify()")

    # Return the ion_table_dir and ion_table_file so they
    # can be set as trident global variables for future use by ion_balance 
    # classes/functions
    return datadir, datafile

def get_datafiles(datadir=None, url=None):
    """
    Downloads an ion table datafile through interactive prompts.  

    **Parameters**

    :datadir: string, optional

        The directory to which to download the ion table datafile.
        If None, uses ``$HOME/.trident``.
        Default: None

    :url: string, optional

        The url to contact to get the data files.  If None, uses
        ``http://trident-project.org/data/ion_table/``     
        Default: None

    **Example**

    >>> get_datafiles()
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
        print("Cannot seem to access %s; Do you have internet access?" % url)
        raise

    line_list = page.split('\n')
    i = 1 # counter
    filenames = []
    print("The following data files are available:")
    for line in line_list:
       if not line.startswith('<!--X'): # identifies files by comments
            continue
       line = line[6:-3]  # strips off the HTML comment chars
       filenames.append(line.split()[0])
       print("%d) %s" % (i, line))
       i += 1

    # User chooses which file
    print("")
    print("Which number would you like to download and use? [1]")
    value = raw_input().rstrip()
    if value == '':
        value = '1'
    while 1:
        try:
            filename = filenames[int(value)-1] 
            break
        except IndexError:
            print("%d is not a valid option.  Please pick one of the listed numbers.")

    # Actually downloads the file to a temporary directory before gunzipping.
    fileurl = os.path.join(url, filename)
    tempdir = tempfile.mkdtemp()
    print("")
    download_file(fileurl, local_directory=tempdir)
    print("  Unzipping file: %s" % filename)
    gunzip_file(os.path.join(tempdir, filename), 
                out_filename=os.path.join(datadir, filename[:-3]))
    shutil.rmtree(tempdir) 
    return filename[:-3]

def trident():
    """
    Print a Trident ASCII logo to the screen.
    """
    print("""
MMMMMMMMMMMMMMMMMMM.............................................................
M...B....C....D...M.MMMMMMMMMM.......MM........MM.....................MM7.......
M...MM...M...MM...M.....MM.....................MM.....................MM7.......
M....MM..M..MM....M.....MM...MMMMMM..MM..MMMMMMMM..MMMMMMM..MMDMMMMM MMMMMMM....
M....MM..M..MM....M.....MM...MM ..MM.MM. MM....MM.,MM...MM+.MM....MM..MM7.......
M.....MMMMMMM.....M.....MM...MM......MM. MM....MM.8MMMMMMMD.MM....MM..MM7.......
M....... M........M.....MM...MM......MM..MM...8MM. MM...._..MM....MM..MMZ.MM ...
M........M........M.....MM...MM......MM..MMMMM.MM. MMMMMMM..MM....MM...MMMMM....
MMMMMMMMMMMMMMMMMMM.............................................................
""")

def trident_path():
    """
    Return the path where the trident source is installed.
    Useful for identifying where data files are (e.g. path/data).  Note that
    ion table datafiles are downloaded separate and placed in another
    location according to the ~/.trident/config.tri file.

    **Example**

    >>> print trident_path()
    """
    path_list = os.path.dirname(__file__).split('/')[:-1]
    path_list.append('trident')
    return '/'.join(path_list)

def ion_table_filepath():
    """
    Return the path of the default trident ion table datafile.
    """
    ion_table_dir, ion_table_file = parse_config()
    return os.path.join(ion_table_dir, ion_table_file)

def make_onezone_dataset(density=1e-26, temperature=1000, metallicity=0.3, 
                         domain_width=10.):
                           
    """
    Create a one-zone hydro dataset for use as test data.  The dataset
    consists of a single cubicle cell of gas with hydro quantities specified in
    the function kwargs.  It makes an excellent test dataset through 
    which to send a sightline and test Trident's capabilities for making 
    absorption spectra.

    Using the defaults and passing a ray through the full domain should
    result in a spectrum with a good number of absorption features.

    **Parameters**

    :density: float, optional

        The gas density value of the dataset in g/cm**3
        Default: 1e-26

    :temperature: float, optional

        The gas temperature value of the dataset in K
        Default: 10**3

    :metallicity: float, optional

        The gas metallicity value of the dataset in Zsun
        Default: 0.3

    :domain_width: float, optional

        The width of the dataset in kpc
        Default: 10.

    **Returns**

    **Example**

    Create a simple one-zone dataset, pass a ray through it, and generate
    a COS spectrum for that ray.

    >>> import trident
    >>> ds = trident.make_onezone_dataset()
    >>> ray = trident.make_simple_ray(ds,
    ...         start_position=ds.domain_left_edge,
    ...         end_position=ds.domain_right_edge,
    ...         fields=['density', 'temperature', 'metallicity'])
    >>> sg = trident.SpectrumGenerator('COS')
    >>> sg.make_spectrum(ray)
    >>> sg.plot_spectrum('spec_raw.png')
    """
    one = np.ones([1,1,1])
    zero = np.zeros([1,1,1])
    dens = np.array([[[density]]])*one
    temp = np.array([[[temperature]]])*one
    metal = np.array([[[metallicity]]])*one
    domain_width *= 1e3*pc/cm
    bbox = np.array([[0., domain_width], [0., domain_width], [0., domain_width]])

    data = {'density':dens, 'temperature':temp, 'metallicity':metal*Zsun,
            'velocity_x':zero, 'velocity_y':zero, 'velocity_z':zero}
    return load_uniform_grid(data, one.shape, length_unit='cm',
                              mass_unit='g', bbox=bbox)

def make_onezone_ray(density=1e-26, temperature=1000, metallicity=0.3, 
                     length=10, redshift=0, filename='ray.h5', 
                     column_densities=None):
    """
    Create a one-zone ray object for use as test data.  The ray
    consists of a single absorber of hydrodynamic characteristics 
    specified in the function kwargs.  It makes an excellent test dataset 
    to test Trident's capabilities for making absorption spectra.

    You can specify the column densities of different ions explicitly using
    the column_densities keyword, or you can let Trident calculate the 
    different ion columns internally from the density, temperature, and 
    metallicity fields.

    Using the defaults will produce a ray that should result in a spectrum 
    with a good number of absorption features.

    **Parameters**

    :density: float, optional

        The gas density value of the ray in g/cm**3
        Default: 1e-26

    :temperature: float, optional

        The gas temperature value of the ray in K
        Default: 10**3

    :metallicity: float, optional

        The gas metallicity value of the ray in Zsun
        Default: 0.3

    :length: float, optional

        The length of the ray in kpc
        Default: 10.

    :redshift: float, optional

        The redshift of the ray
        Default: 0

    :filename: string, optional

        The filename to which to save the ray to disk.  Due to the 
        mechanism for passing rays, the ray data must be saved to disk.
        Default: 'ray.h5'

    :column_densities: dict, optional

        The user can create a dictionary which adds more number density ion 
        fields to the ray.  Each key in the dictionary should be the desired 
        ion field name according to the field name format:
        i.e.  "<ELEMENT>_p<IONSTATE>_number_density" 
        e.g. neutral hydrogen = "H_p0_number_density".  
        The corresponding value for each key should be the desired column 
        density of that ion in cm**-2.  See example below.
        Default: None

    **Returns**

        A YT LightRay object

    **Example**

    Create a one-zone ray, and generate a COS spectrum from that ray.

    >>> import trident
    >>> ray = trident.make_onezone_ray()
    >>> sg = trident.SpectrumGenerator('COS')
    >>> sg.make_spectrum(ray)
    >>> sg.plot_spectrum('spec_raw.png')

    Create a one-zone ray with an HI column density of 1e21 (DLA) and generate
    a COS spectrum from that ray for just the Lyman alpha line.

    >>> import trident
    >>> ds = trident.make_onezone_ray(column_densities={'H_number_density': 1e21})
    >>> sg = trident.SpectrumGenerator('COS')
    >>> sg.make_spectrum(ray, lines=['Ly a'])
    >>> sg.plot_spectrum('spec_raw.png')
    """
    from yt import save_as_dataset
    length = YTArray([length], "kpc")
    data = {"density"            : YTArray([density], "g/cm**3"),
            "metallicity"        : YTArray([metallicity], "Zsun"),
            "dl"                 : length,
            "temperature"        : YTArray([temperature], "K"),
            "redshift"           : np.array([redshift]),
            "redshift_eff"       : np.array([redshift]),
            "velocity_los"       : YTArray([0.], "cm/s"),
            "x": length/2, "dx": length,
            "y": length/2, "dy": length,
            "z": length/2, "dz": length
            }

    extra_attrs = {"data_type": "yt_light_ray"}
    field_types = dict([(field, "grid") for field in data.keys()])

    # Add additional number_density fields to dataset
    if column_densities:
        for k,v in column_densities.iteritems():
            # Assure we add X_number_density for neutral ions
            # instead of X_p0_number_density
            key_string_list = k.split('_')
            if key_string_list[1] == 'p0':
                k = '_'.join([key_string_list[0], key_string_list[2], 
                              key_string_list[3]])
            v = YTArray([v], 'cm**-2')
            data[k] = v/length
            field_types[k] = 'grid'

    ds = {"current_time": 0.,
          "current_redshift": 0.,
          "cosmological_simulation": 0.,
          "domain_left_edge": np.zeros(3)*length,
          "domain_right_edge": np.ones(3)*length,
          "periodicity": [True]*3}

    save_as_dataset(ds, filename, data, field_types=field_types,
                       extra_attrs=extra_attrs)

    # load dataset and make spectrum
    ray = load(filename)
    return ray

def verify():
    """
    Verify that the bulk of Trident's functionality is working.  First, it
    ensures that the user has a configuration file and ion table datafile, 
    and creates/downloads these files if they do not exist.  Next, it
    creates a single-cell grid-based dataset in memory, generates a ray 
    by sending a sightline through that dataset, then makes a spectrum from 
    the ray object.  It saves all data to a tempdir before deleting it.
    """
    parse_config()
    from trident.spectrum_generator import SpectrumGenerator
    from trident.ray_generator import make_simple_ray
    print("")
    print("Creating single-cell dataset")
    print("----------------------------")
    print("")
    try:
        ds = make_onezone_dataset()
    except:
        print("Failed to create single-cell dataset")
        raise

    print("")
    print("Creating ray object through single-cell dataset")
    print("-----------------------------------------------")
    print("")
    tempdir = tempfile.mkdtemp()
    try:
        ray = make_simple_ray(ds,
                start_position=ds.domain_left_edge,
                end_position=ds.domain_right_edge,
                data_filename=os.path.join(tempdir, 'ray.h5'),
                fields=['density', 'temperature', 'metallicity'])
    except:
        print("Failed to create ray object")
        raise

    print("")
    print("Create spectrum with Lyman alpha, Mg II, and O VI lines")
    print("-------------------------------------------------------")
    print("")
    sg = SpectrumGenerator('COS')
    sg.make_spectrum(ray, lines=['Ly a', 'Mg II', 'O VI'])
    sg.save_spectrum(os.path.join(tempdir, 'spec_raw.h5'))
    sg.plot_spectrum(os.path.join(tempdir, 'spec_raw.png'))
    # Test other post processing of the spectrum
    sg.add_qso_spectrum()
    sg.add_milky_way_foreground()
    sg.apply_lsf()
    sg.add_gaussian_noise(30)
    sg.save_spectrum(os.path.join(tempdir, 'spec_final.h5'))
    sg.plot_spectrum(os.path.join(tempdir, 'spec_final.png'))
    print("Removing all temporary data files...")
    shutil.rmtree(tempdir) 
    print("")
    print("Congratulations, you have verified that Trident is installed correctly.")
    print("Now let's science!")
    print("")
