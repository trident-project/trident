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
    A class representing a spectrograph's line spread function.

    A line spread function can be defined either by a function and width
    *or* by a filename containing a custom kernel.

    **Parameters**

    :function: string, optional

        The function defining the LSF kernel.
        valid functions are "boxcar" or "gaussian"

    :width: int, optional

        The width of the LSF kernel in bins.  

    :filename: string, optional

        The filename of a textfile for a user-specified kernel. Each line
        in the textfile contains a normalized flux value of the kernel.
        For examples, see contents of ``trident.__path__/data/lsf_kernels``
        Trident searches for these files either in the aforementioned 
        directory or in the execution directory.

    **Examples**

    Generate an LSF based on a text file:

    >>> LSF(filename='avg_COS.txt')

    Generate a boxcar-based LSF:

    >>> LSF(function='boxcar', width=30)

    Generate a gaussian-based LSF:

    >>> LSF(function='guassian', width=7)
    """
    def __init__(self, function=None, width=None, filename=None):
        self.kernel = []
        self.filename = filename
        self.function = function
        self.width = width
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
                    raise RuntimeError("LSF filename not found in current " +
                        "directory or in %s/data/lsf_kernels directory" % 
                        trident_path())
            for line in lsf_file:
                self.kernel.append(float(line.split()[1]))
            lsf_file.close()
            self.kernel = np.array(self.kernel)
            self.width = self.kernel.size
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
            raise RuntimeError("Either LSF filename OR function+width must be specified.")

    def __repr__(self):
        if self.filename is not None:
            return "LSF from file %s" % self.filename
        else:
            return "LSF %s of width %d" % (self.function, self.width)
