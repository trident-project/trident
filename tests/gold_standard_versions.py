"""
Testing utilities for Trident

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import print_function
from trident import path as trident_path
import os

def get_gold_standard_version():
    f = open(os.path.join(trident_path, '../.travis.yml'), 'r')
    for line in f:
        line = line.lstrip()
        if line.startswith('YT_GOLD'):
            line_list = line.split('=')
            yt_gold = line_list[1]
        if line.startswith('TRIDENT_GOLD'):
            line_list = line.split('=')
            trident_gold = line_list[1]
    f.close()
    print()
    print('Latest Gold Standard Commit Tags\n')
    print('yt = %s' % yt_gold, end='')
    print('Trident = %s' % trident_gold)
    print('To update to them, `git checkout <tag>` in appropriate repository')

if __name__ == "__main__":
    get_gold_standard_version()
