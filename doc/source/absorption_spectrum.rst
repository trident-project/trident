.. _absorption_spectrum:

AbsorptionSpectrum
==================

For documentation on the main interface to spectrum creation in Trident,
see :ref:`spectrum-generation`.

The :class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
is the internal class for creating absorption spectra in Trident from
:class:`~trident.LightRay` objects. The
:class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
and its workhorse method
:meth:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum.make_spectrum`
return two arrays, one with wavelengths, the other with the normalized
flux values at each of the wavelength values.  It can also output a text file
listing all important lines.

Method for Creating Absorption Spectra
--------------------------------------

Once a :class:`~trident.LightRay`
has been created traversing a dataset using the :ref:`light-ray-generator`,
a series of arrays store the various fields of the gas parcels (represented
as cells) intersected along the ray.
:class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
steps through each element of the
:class:`~trident.LightRay`'s
arrays and calculates the column density for desired ion by multiplying its
number density with the path length through the cell.  Using these column
densities along with temperatures to calculate thermal broadening, voigt
profiles are deposited on to a featureless background spectrum.  By default,
the peculiar velocity of the gas is included as a doppler redshift in addition
to any cosmological redshift of the data dump itself.

Subgrid Deposition
^^^^^^^^^^^^^^^^^^

For features not resolved (i.e. possessing narrower width than the spectral
resolution),
:class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
performs subgrid deposition.  The subgrid deposition algorithm creates a number
of smaller virtual bins, by default the width of the virtual bins is 1/10th
the width of the spectral feature.  The Voigt profile is then deposited
into these virtual bins where it is resolved, and then these virtual bins
are numerically integrated back to the resolution of the original spectral bin
size, yielding accurate equivalent widths values.
:class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
informs the user how many spectral features are deposited in this fashion.

Creating an Absorption Spectrum
-------------------------------

Initialization
^^^^^^^^^^^^^^

To instantiate an
:class:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum`
object, the arguments required are the
minimum and maximum wavelengths (assumed to be in Angstroms), and the number
of wavelength bins to span this range (including the endpoints)

.. code-block:: python

   from trident.absorption_spectrum.absorption_spectrum import AbsorptionSpectrum

   sp = AbsorptionSpectrum(900.0, 1800.0, 10001)

Adding Features to the Spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Absorption lines and continuum features can then be added to the spectrum.
To add a line, you must know some properties of the line: the rest wavelength,
f-value, gamma value, and the atomic mass in amu of the atom.  That line must
be tied in some way to a field in the dataset you are loading, and this field
must be added to the LightRay object when it is created.  Below, we will
add the H Lyman-alpha line, which is tied to the neutral hydrogen field
('H_number_density').

.. code-block:: python

  my_label = 'HI Lya'
  field = 'H_number_density'
  wavelength = 1215.6700 # Angstroms
  f_value = 4.164E-01
  gamma = 6.265e+08
  mass = 1.00794

  sp.add_line(my_label, field, wavelength, f_value, gamma, mass, label_threshold=1.e10)

In the the call to
:meth:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum.add_line`
the ``field`` argument tells the spectrum generator which
field from the ray data to use to calculate the column density.  The
``label_threshold`` keyword tells the spectrum generator to add all lines
above a column density of 10 :superscript:`10` cm :superscript:`-2` to the
text line list output at the end.  If None is provided, as is the default,
no lines of this type will be added to the text list.

Continuum features with optical depths that follow a power law can be added
with the
:meth:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum.add_continuum`
function.  Like adding lines, you must specify details like the wavelength
and the field in the dataset and LightRay that is tied to this feature.
The wavelength refers to the location at which the continuum begins to be 
applied to the dataset, and as it moves to lower wavelength values, the 
optical depth value decreases according to the defined power law.  The 
normalization value is the column density of the linked field which results
in an optical depth of 1 at the defined wavelength.  Below, we add the hydrogen 
Lyman continuum.

.. code-block:: python

  my_label = 'HI Lya'
  field = 'H_number_density'
  wavelength = 912.323660 # Angstroms
  normalization = 1.6e17
  index = 3.0

  sp.add_continuum(my_label, field, wavelength, normalization, index)

Making the Spectrum
^^^^^^^^^^^^^^^^^^^

Once all the lines and continuua are added, the spectrum is made with the
:meth:`~trident.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum.make_spectrum`
function.

.. code-block:: python

  wavelength, flux = sp.make_spectrum('lightray.h5',
                                      output_file='spectrum.fits',
                                      line_list_file='lines.txt')

A spectrum will be made using the specified ray data and the wavelength and
flux arrays will also be returned.  If you set the optional
``use_peculiar_velocity`` keyword to False, the lines will not incorporate
doppler redshifts to shift the deposition of the line features.

Three output file formats are supported for writing out the spectrum: fits,
hdf5, and ascii.  The file format used is based on the extension provided
in the ``output_file`` keyword: ``.fits`` for a fits file,
``.h5`` for an hdf5 file, and anything else for an ascii file.

.. note:: To write out a fits file, you must install the `astropy <http://www.astropy.org>`_ python library in order to access the astropy.io.fits module.  You can usually do this by simply running `pip install astropy` at the command line.

Generating Spectra in Parallel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spectrum generation is parallelized using a multi-level
strategy where each absorption line is deposited by a different processor.
If the number of available processors is greater than the number of lines,
then the deposition of individual lines will be divided over multiple
processors.

Absorption spectrum creation can be run in parallel simply by adding the following
to the top of the script and running with ``mpirun``.

.. code-block:: python

   import yt
   yt.enable_parallelism()

For more information on parallelism in yt, see
`Parallel Computation With yt
<http://yt-project.org/docs/dev/analyzing/parallel_computation.html>`__.
