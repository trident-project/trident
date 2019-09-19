"""
Trident is a python-based package extending yt for generating synthetic
observations of hydrodynamical datasets.

Please visit the Trident and yt websites for more information:
http://trident-project.org
http://yt-project.org
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = "1.2"

# Must run import_check() before anything else is imported to avoid
# astropy error when importing trident in trident package directory
from trident.utilities import import_check
import_check()

from trident.config import \
    parse_config, \
    create_config, \
    trident, \
    trident_path, \
    verify

from trident.ion_balance import \
    add_ion_fields, \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field, \
    solar_abundance, \
    atomic_mass

from trident.instrument import \
    Instrument

from trident.line_database import \
    LineDatabase, \
    Line

from trident.lsf import \
    LSF

from trident.plotting import \
    plot_spectrum

from trident.spectrum_generator import \
    SpectrumGenerator, \
    valid_instruments, \
    load_spectrum

from trident.utilities import \
    make_onezone_dataset, \
    make_onezone_ray

from trident.ray_generator import \
    make_simple_ray, \
    make_compound_ray

from trident.roman import \
    to_roman, \
    from_roman

from trident.light_ray import \
    LightRay

# Making installation path global
path = trident_path()
