"""
Trident is a python-based package extending yt for generating synthetic
observations of hydrodynamical datasets.

Please visit the Trident and yt websites for more information:
http://trident-project.org
http://yt-project.org
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from ion_balance import \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field, \
    solar_abundance, \
    atomic_mass

from instrument import \
    Instrument

from line_database import \
    LineDatabase, \
    Line

from lsf import \
    LSF

from plotting import \
    plot_spectrum

from spectrum_generator import \
    SpectrumGenerator, \
    valid_instruments

from utilities import \
    print_logo

# Import yt's LightRay class
from yt.analysis_modules.cosmological_observation.api import \
    LightRay

from ray_generator import \
    make_simple_ray, \
    make_compound_ray

print_logo()
