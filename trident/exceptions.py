"""
Trident exceptions

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

class TridentError(Exception):
    """
    Trident Base Class Error
    """
    pass

class NotConfiguredError(TridentError):
    """
    When Trident hasn't been configured properly with a config file.
    """
    pass

class NoIonBalanceTableError(TridentError):
    """
    When Trident is missing its ion balance file and needs it.
    """
    pass

