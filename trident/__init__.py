"""
Trident is a python-based package extending yt for generating synthetic
observations of hydrodynamical datasets.

Please visit the Trident and yt websites for more information:
http://trident-project.org
http://yt-project.org
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2106, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = "0.5.1"

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
    valid_instruments

from trident.utilities import \
    parse_config, \
    create_config, \
    trident, \
    trident_path, \
    make_onezone_dataset, \
    make_onezone_ray, \
    verify

from trident.ray_generator import \
    make_simple_ray, \
    make_compound_ray

from trident.roman import \
    to_roman, \
    from_roman

# Import yt's LightRay class
from yt.analysis_modules.cosmological_observation.api import \
    LightRay

# Making installation path global
path = trident_path()
