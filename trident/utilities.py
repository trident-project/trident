"""
Miscellaneous Utilities for Trident

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os

def trident():
    """
    Prints a nice ASCII logo!
    """
    print """
MMMMMMMMMMMMMMMMMMM.............................................................
M...B....C....D...M.MMMMMMMMMM.......MM........MM.....................MM7.......
M...MM...M...MM...M.....MM.....................MM.....................MM7.......
M....MM..M..MM....M.....MM...MMMMMM..MM..MMMMMMMM..MMMMMMM..MMDMMMMM MMMMMMM....
M....MM..M..MM....M.....MM...MM ..MM.MM. MM....MM.,MM...MM+.MM....MM..MM7.......
M.....MMMMMMM.....M.....MM...MM......MM. MM....MM.8MMMMMMMD.MM....MM..MM7.......
M....... M........M.....MM...MM......MM..MM...8MM. MM...._..MM....MM..MMZ.MM ...
M........M........M.....MM...MM......MM..MMMMM.MM. MMMMMMM..MM....MM...MMMMM....
MMMMMMMMMMMMMMMMMMM.............................................................
Please send questions to trident-project-users@googlegroups.com
"""

def trident_path():
    """
    A function returning the path of the trident installation directory.
    """
    return '/'.join(os.path.dirname(__file__).split('/')[:-1])
