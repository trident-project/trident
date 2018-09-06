"""
Trident config

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
from six.moves.configparser import \
    SafeConfigParser
from six.moves import \
    input
import shutil
import tempfile
import sys

from trident.utilities import \
    ensure_directory, \
    get_datafiles, \
    make_onezone_dataset

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

    >>> print(trident_path())
    """
    # os.path.split(__file__)returns a tuple with [0] the path to the file
    # and [1] the filename.
    # Here, __file__ refers to this file (config.py)
    return os.path.split(__file__)[0]

def create_config():
    """
    Create a Trident configuration file using interaction with the user.
    This function is called by :class:`~trident.parse_config` if it appears
    that the configuration has not yet been set up for the user.  It will
    attempt to create a configuration file and download an ion table
    datafile from the web.  It does this using user interaction from the
    python prompt.
    """
    default_dir = os.path.expanduser('~/.trident')
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
    value = input().rstrip()
    if not value == '' and not value == 'y':
        sys.exit('Instructions at http://trident.readthedocs.org/en/latest/installation.html')

    print("")
    print("Where would you like Trident to store the ion table file?")
    print("[%s]" % default_dir)

    # First assure that the .trident directory is created for storing
    # the config file.
    ensure_directory(default_dir)

    datadir = input().rstrip()
    if datadir == '':
        datadir = default_dir
    datadir = os.path.expanduser(datadir)

    # Try to create data directory if it doesn't exist
    try:
        ensure_directory(datadir)
        print("Using %s" % datadir)
        print("")
    except BaseException:
        print("Cannot create directory %s" % datadir)
        raise

    # Get the datafile from the web
    datafile = get_datafiles(datadir=datadir)

    # Create the config file and make it to the datadir and datafiles chosen
    config = SafeConfigParser()
    config.add_section('Trident')
    config.set('Trident', 'ion_table_dir', datadir)
    config.set('Trident', 'ion_table_file', datafile)
    config_filename = os.path.expanduser('~/.trident/config.tri')
    with open(config_filename, 'w') as configfile:
        config.write(configfile)

    print("")
    print("Installation complete.  I recommend verifying your installation")
    print("to assure that everything is working.  Try: trident.verify()")

    # Return the config file path so we can load it and get parameters.
    return config_filename

def parse_config(variable=None):
    """
    Parse the Trident local configuration file.  This function is called
    whenever Trident is imported, and it assures that Trident knows where
    the default ion table datafiles exist.  If a ``config.tri`` file doesn't
    exist in ``$HOME/.trident`` or in the current working directory, then
    Trident will launch the :class:`~trident.create_config` function to
    try to automatically generate one for the user.  For more information
    on this process, see the installation documentation.

    **Parameters**

    :variable: string, optional

        If you wish to get the value a variable is set to in the config
        file, specify that variable name here.  Will return the result
        value of that variable. Default: None
    """
    # Assure the ~/.trident directory exists, and read in the config file.
    home = os.path.expanduser("~")
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
    except BaseException:
        config_filename = create_config()
        parser = SafeConfigParser()
        parser.read(config_filename)
        ion_table_dir = parser.get('Trident', 'ion_table_dir')
        ion_table_file = parser.get('Trident', 'ion_table_file')

    ion_table_dir = os.path.abspath(os.path.expanduser(ion_table_dir))
    if not os.path.exists(os.path.join(ion_table_dir,
                                       ion_table_file)):
        print("")
        print("No ion table data file found in %s" % ion_table_dir)
        ion_table_file = get_datafiles(ion_table_dir)
        parser.set('Trident', 'ion_table_file', ion_table_file)
        with open(config_filename, 'w') as configfile:
            parser.write(configfile)

    # value to return depends on what was set for "variable"
    if variable is None:
        return ion_table_dir, ion_table_file
    else:
        return parser.get('Trident', variable)

def verify(save=False):
    """
    Verify that the bulk of Trident's functionality is working.  First, it
    ensures that the user has a configuration file and ion table datafile,
    and creates/downloads these files if they do not exist.  Next, it
    creates a single-cell grid-based dataset in memory, generates a ray
    by sending a sightline through that dataset, then makes a spectrum from
    the ray object.  It saves all data to a tempdir before deleting it.

    **Parameters**

    :save: boolean, optional

        By default, verify saves all of its outputs to a temporary directory
        and then removes it upon completion.  If you would like to see the
        resulting data from verify(), set this to be True and it will save
        a light ray, and raw and processed spectra in the current working
        directory.
        Default: False

    **Example**

    Verify Trident works.

    >>> import trident
    >>> trident.verify()
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
    except BaseException:
        print("Failed to create single-cell dataset")
        raise

    print("")
    print("Creating ray object through single-cell dataset")
    print("-----------------------------------------------")
    print("")

    if save:
        tempdir = '.'
    else:
        tempdir = tempfile.mkdtemp()

    try:
        ray = make_simple_ray(ds,
                start_position=ds.domain_left_edge,
                end_position=ds.domain_right_edge,
                data_filename=os.path.join(tempdir, 'ray.h5'),
                fields=['density', 'temperature', 'metallicity'])
    except BaseException:
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

    if not save:
        print("Removing all temporary data files...")
        shutil.rmtree(tempdir)

    print("")
    print("Congratulations, you have verified that Trident is installed correctly.")
    print("Now let's science!")
    print("")

# Each time Trident is imported, we determine the settings from the config
# file or try to create a config file.  But don't do this on readthedocs, or
# it will fail in the build. In readthedocs environment, just set a dummy
# filepath so readthedocs can parse the docstrings OK.


on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    ion_table_dir = trident_path()
    ion_table_file = '__init__.py'
    ion_table_filepath = os.path.join(ion_table_dir, ion_table_file)
else:
    ion_table_dir, ion_table_file = parse_config()
    ion_table_filepath = os.path.join(ion_table_dir, ion_table_file)
