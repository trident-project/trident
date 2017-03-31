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
    except:
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

ion_table_dir, ion_table_file = parse_config()
ion_table_filepath = os.path.join(ion_table_dir, ion_table_file)

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
    config_filename = os.path.expanduser('~/.trident/config.tri')
    with open(config_filename, 'w') as configfile:
        config.write(configfile)

    print("")
    print("Installation complete.  I recommend verifying your installation")
    print("to assure that everything is working.  Try: trident.verify()")

    # Return the config file path so we can load it and get parameters.
    return config_filename
