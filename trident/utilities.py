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

import gzip
import os
from os.path import \
    expanduser
import six
import requests
import tempfile
import shutil
from yt.funcs import \
    get_pbar
import numpy as np
from six.moves import input
from yt.units import \
    cm, \
    pc, \
    Zsun
from yt import \
    load_uniform_grid, \
    YTArray, \
    load
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.funcs import \
    mylog

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
        local_filename = url.split(os.sep)[-1]
    if local_directory is None:
        local_directory = '.'
    ensure_directory(local_directory)

    # open local file handle
    filepath = os.path.join(local_directory, local_filename)
    filehandle = open(filepath, 'wb')

    # Get information about remote filesize
    r = requests.get(url, stream=True)
    filesize = int(r.headers["content-length"])/2**10  # in kB
    if progress_bar:
        pbar = get_pbar("Downloading file: %s" % local_filename, filesize)
    filesize_dl = 0
    chunk_size = 8192

    # Download file in chunks and update statusbar until done with transfer
    for content in r.iter_content(chunk_size):
        filesize_dl += len(content)/2**10
        filehandle.write(content)
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
    except BaseException:
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
    value = input().rstrip()
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

    extra_attrs = {"data_type": "yt_light_ray", "dimensionality": 3}
    field_types = dict([(field, "grid") for field in data.keys()])

    # Add additional number_density fields to dataset
    if column_densities:
        for k,v in six.iteritems(column_densities):
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

def import_check():
    """
    """
    # Avoid astropy error when importing from trident package directory.
    plist = os.path.dirname(os.path.abspath(__file__)).split(os.sep)
    package_path = os.sep.join(plist[:-1])
    if os.getcwd() == package_path:
        raise RuntimeError(
            """

The Trident package does not work correctly when imported from its
installation directory.  Please try moving to another directory.""")

def _determine_dataset_sampling_type(ds):
    """
    Determine whether the dataset is particle-based or grid-based.

    LightRays datasets are reloaded as particle type regardless of the
    underlying frontend, and they should always be treated as grid.
    """
    # Particle-based datasets like Gadget, Gizmo, Gasoline, Changa
    if isinstance(ds.index, ParticleIndex):
        part_type = True
    else:
        part_type = False

    # LightRays should always be treated as grid datasets.
    # Sometimes data_type is not defined.
    if getattr(ds, "data_type", None) == "yt_light_ray":
        part_type = False
    if part_type:
        return "particle"
    else:
        return "cell"

def _check_sampling_types_match(ds, ftype):
    """
    Checks if sampling_type of field and dataset matches.

    Users may not always put the correct ftype for their datasets, and this
    can lead to errors if a user specifies a grid-based field for a
    particle-dataset, which will lead to adding ion fields to deposited (and
    already interpolated) fields with unexpected behavior.  This alerts the
    user when it thinks this has happened.
    """
    # Determine the sampling type (e.g., "particle" or "cell") for the
    # dataset as a whole.
    sampling_type = _determine_dataset_sampling_type(ds)

    # Determine the sampling type (e.g., "particle" or "cell") for the
    # fields of "ftype" specified by the user.
    ds.index
    field_sampling_type = ds.field_info[ftype, 'density'].sampling_type

    if sampling_type != field_sampling_type:
        mylog.warning("===================================================")
        mylog.warning("MISMATCH BETWEEN SAMPLING_TYPE OF FTYPE AND DATASET")
        mylog.warning("sampling_type of (%s, 'density') = %s" % (ftype, field_sampling_type))
        mylog.warning("sampling_type of dataset = %s" % sampling_type)
        mylog.warning("THIS IS PROBABLY UNDESIRED BEHAVIOR.  PLEASE CHOOSE A DIFFERENT FTYPE.")
        mylog.warning("===================================================")

    return sampling_type
