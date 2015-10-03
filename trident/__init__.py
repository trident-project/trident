"""
Trident is a python-based package extending yt for generating synthetic
observations of hydrodynamical datasets.

Please visit the Trident and yt websites for more information:
http://trident-project.org
http://yt-project.org
"""

# Import all structures from trident.py
from code import *

from ion_balance import \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field, \
    solar_abundance, \
    atomic_mass

from line_database import LineDatabase, Line

# Import yt's LightRay class
from yt.analysis_modules.cosmological_observation.api import LightRay
