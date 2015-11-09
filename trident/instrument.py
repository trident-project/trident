"""
Instrument class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

class Instrument(object):
    """
    An instrument template for specifying a spectrograph/telescope pair

    **Parameters**

    lambda_min : int
        Minimum desired wavelength for generated spectrum (in angstroms)

    lambda_max : int
        Maximum desired wavelength for generated spectrum (in angstroms)

    n_lambda : int
        Number of desired wavelength bins for the spectrum
        Default: None

    dlambda : float
        Desired bin width for the spectrum
        <note from Devin: which one supercedes the others?>
        Default: None

    lsf_kernel : string
        The filename for the LSF kernel
        Default: None

    name : string
        Name assigned to the Instrument object
        Default: None

    """
    def __init__(self, lambda_min, lambda_max, n_lambda=None,
                 dlambda=None, lsf_kernel=None, name=None):
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.lsf_kernel = lsf_kernel
        if n_lambda is None and dlambda is None:
            raise RuntimeError("Either n_lambda or dlambda must be set to "
                               "specify the binsize")
        elif dlambda is not None:
            # Add 1 to n_lambda to make sure dlambda preserved
            n_lambda = ((lambda_max - lambda_min) / dlambda) + 1
        self.n_lambda = n_lambda
        if name is not None:
            self.name = name
