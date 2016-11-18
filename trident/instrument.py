"""
Instrument class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

class Instrument(object):
    """
    An instrument class for specifying a spectrograph/telescope pair

    **Parameters**

    :lambda_min: int

        Minimum desired wavelength for generated spectrum in angstroms

    :lambda_max: int

        Maximum desired wavelength for generated spectrum in angstroms

    :n_lambda: int

        Number of desired wavelength bins for the spectrum
        Setting dlambda overrides n_lambda value
        Default: None

    :dlambda: float

        Desired bin width for the spectrum in angstroms
        Setting dlambda overrides n_lambda value
        Default: None

    :lsf_kernel: string

        The filename for the :class:`~trident.LSF` kernel
        Default: None

    :name: string

        Name assigned to the :class:`~trident.Instrument` object
        Default: None

    """
    def __init__(self, lambda_min, lambda_max, n_lambda=None,
                 dlambda=None, lsf_kernel=None, name=None):
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.lsf_kernel = lsf_kernel
        self.name = name
        if n_lambda is None and dlambda is None:
            raise RuntimeError("Either n_lambda or dlambda must be set to "
                               "specify the binsize")
        elif dlambda is not None:
            # adding 1 here to assure we cover full lambda range
            n_lambda = (lambda_max - lambda_min) / float(dlambda) + 1
        self.n_lambda = n_lambda
        if dlambda is None:
            # adding 1 here to assure we cover full lambda range
            dlambda = (lambda_max - lambda_min) / float(n_lambda - 1)
        self.dlambda = dlambda

    def __repr__(self):
        disp = "<Instrument>:\n"
        disp += "    name: %s\n" % self.name
        disp += "    lambda_min: %f\n" % self.lambda_min
        disp += "    lambda_max: %f\n" % self.lambda_max
        disp += "    n_lambda: %d\n" % self.n_lambda
        disp += "    dlambda: %f\n" % self.dlambda
        disp += "    lsf_kernel: %s" % self.lsf_kernel
        return disp
