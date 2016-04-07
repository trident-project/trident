"""
LSF class and member functions for Line Spread Functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
import sys
from yt.funcs import \
    mylog
from trident.utilities import \
    trident_path

class LSF(object):
    """
    Line Spread Function class

    The user must define either a filename or a function and a width

    Parameters

    function : string, optional
        the function defining the LSF kernel.
        valid functions are "boxcar" or "gaussian"

    width : int, optional
        the width of the LSF kernel.  

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
                filename2 = os.path.join(trident_path(), "data", \
                                         "lsf_kernels", filename)
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
                if width % 2 == 0:
                    mylog.warn("LSF kernel must have an odd length. Reducing kernel size by 1.")
                    width -= 1
                self.kernel = np.ones(width)/width
            elif function == 'gaussian':
                from astropy.convolution import Gaussian1DKernel
                self.kernel = Gaussian1DKernel(width)
        else:
            sys.exit("Either filename OR function+width must be specified.")

    def __repr__(self):
        if self.filename is not None:
            return "from file %s" % self.filename
        else:
            return "%s of width %d" % (function, width)
