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

from trident import path as trident_path
import os
import re

def get_gold_standard_version():
    f = open(os.path.join(trident_path, '../.circleci/config.yml'), 'r')
    for line in f:
        line = line.lstrip()
        pattern_YT = "echo 'YT_GOLD=(\S+)'"
        pattern_trident = "echo 'TRIDENT_GOLD=(\S+)'"
        match_YT = re.search(pattern_YT, line)
        match_trident = re.search(pattern_trident, line)
        if match_YT:
            yt_gold = match_YT.group(1)
        if match_trident:
            trident_gold = match_trident.group(1)
    f.close()
    print()
    print('Latest Gold Standard Commit Tags\n')
    print('yt = %s' % yt_gold)
    print('Trident = %s\n' % trident_gold)
    print('To update to them, `git checkout <tag>` in appropriate repository')


if __name__ == "__main__":
    get_gold_standard_version()
