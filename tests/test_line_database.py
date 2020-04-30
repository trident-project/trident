"""
Tests for Line Database classes and functions

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from trident.line_database import \
    uniquify, \
    Line, \
    LineDatabase
from trident.config import \
    trident_path
import numpy as np
from shutil import copyfile
import os.path

def test_uniquify():
    """
    Test uniquify function
    """
    my_list = list(range(3))
    my_list[1] = 0
    l_unique = uniquify(my_list)
    assert len(l_unique) == len(my_list)-1

    arr = np.arange(3)
    arr[1] = 0
    arr_unique = uniquify(arr)
    assert len(arr_unique) == len(arr)-1

    arr = np.ones(5)
    arr_unique = uniquify(arr)
    assert len(arr_unique) == 1

def test_line():
    HI = Line('H', 'I', 1216, 1.5, 2.3, identifier='Ly a')
    assert HI.field == "H_p0_number_density"
    print(HI)

def test_line_database_from_stored_file():
    ld = LineDatabase('lines.txt')
    print(ld)

def test_line_database_from_local_file():
    line_file = os.path.join(trident_path(), 'data/line_lists/lines.txt')
    copyfile(line_file, 'test_lines.txt')
    ld = LineDatabase('test_lines.txt')
    print(ld)
    os.remove('test_lines.txt')

def test_line_database_from_input():
    ld = LineDatabase()
    HI = Line('H', 'I', 1216, 1.5, 2.3, identifier='Ly a')
    ld.add_line('H', 'I', 1216, 1.5, 2.3, identifier='Ly a')
    assert ld.lines_all[0].identifier == HI.identifier
    print(ld)

def test_select_lines_from_line_database():
    ld = LineDatabase('lines.txt')
    assert len(ld.select_lines('Ne')) == 8 # 8 listed Ne lines
    assert len(ld.select_lines('O', 'VI')) == 2 # 2 listed OVI lines
    assert len(ld.select_lines('H', 'I')) == 39 # 39 listed HI lines
    assert len(ld.select_lines('H', 'I', 1216)) == 1 # 1 listed HI lines
    assert len(ld.select_lines(identifier='Ly b')) == 1 # 1 listed Ly beta

def test_parse_subset_from_line_database():
    ld = LineDatabase('lines.txt')
    assert len(ld.parse_subset(['Al'])) == 3
    ld = LineDatabase('lines.txt')
    assert len(ld.parse_subset(['Mg II'])) == 2
    ld = LineDatabase('lines.txt')
    assert len(ld.parse_subset(['O V 630'])) == 1

