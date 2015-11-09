"""
LSF class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import os
import roman
import sys
from yt.convenience import load
from yt.analysis_modules.absorption_spectrum.api import AbsorptionSpectrum
from yt.funcs import mylog, YTArray
from matplotlib import pyplot
from ion_balance import \
    add_ion_fraction_field, \
    add_ion_number_density_field, \
    add_ion_density_field, \
    add_ion_mass_field, \
    atomic_mass
from line_database import LineDatabase

class LSF():
    """
    Line Spread Function class

    The user must define either a filename or a function and a width

    Parameters

    function : string, optional
        the function defining the LSF kernel.
        valid functions are "boxcar" or "gaussian"

    width : int, optional
        the width of the LSF kernel

    filename : string, optional
        the filename of a textfile for a user-specified kernel. each line
        in the textfile is the non-normalized flux value of the kernel
    """
    def __init__(self, function=None, width=None, filename=None):
        self.kernel = []
        # if filename is defined, use it
        if filename is not None:
            # Check to see if the file is in the local dir
            if os.path.isfile(filename):
                lsf_file = open(filename, 'r')
            # otherwise use the file in the lsf_kernels dir
            else:
                filename2 = os.path.join(os.path.dirname(__file__), "..",
                                         "data", "lsf_kernels", filename)
                if os.path.isfile(filename2):
                    lsf_file = open(filename2, 'r')
                else:
                    sys.exit("filename must be in local directory or in",
                             "trident/data/lsf_kernels directory")
            for line in lsf_file:
                self.kernel.append(float(line.split()[1]))
            lsf_file.close()
        elif function is not None and width is not None:
            if function == 'boxcar':
               self.kernel = np.ones(width)
            #XXX Define more functional forms
            #elif function == 'gaussian':
            #   self.kernel = np.gaussian(width)
        else:
            sys.exit("Either filename OR function+width must be specified.")
