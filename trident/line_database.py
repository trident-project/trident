"""
Line, LineDatabase class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

import os
from sets import \
    Set
from yt.funcs import \
    mylog
from trident.utilities import \
    trident_path
from trident.roman import \
    fromRoman

def uniquify(list):
   # order preserving method for reducing duplicates in a list
   checked = []
   for val in list:
       if val not in checked:
           checked.append(val)
   return checked

class Line:
    """An individual atomic transition identified uniquely by element,
    ionic state, wavelength, gamma, oscillator strength, and identifier.

    **Parameters**

    element : string
        The element of the transition using element's symbol on periodic table

    ion_state : string
        The roman numeral representing the ionic state of the transition

    wavelength : float
        The wavelength of the transition in angstroms

    gamma : float
        The gamma of the transition in Hertz

    f_value: float
        The oscillator strength of the transition

    field : string, optional
        The default yt field name associated with the ion responsible for
        this line

    identifier : string, optional
        An optional identifier for the transition


    Example:

    Create a Line object for the neutral hydrogen 1215 Angstroms transition.

    >>> HI = Line('H', 'I', 1215.67, 469860000, 0.41641, 'Lya')

    """
    def __init__(self, element, ion_state, wavelength, gamma, f_value,
                 field=None, identifier=None):
        self.element = element
        self.ion_state = ion_state
        self.wavelength = wavelength
        self.gamma = gamma
        self.f_value = f_value
        self.name = '%s%s %d' % (element, ion_state, round(float(wavelength), 0))
        if identifier is None:
            identifier = self.name
        self.identifier = identifier

        # Automatically populate the field if not defined
        if field is None:
            ion_number = fromRoman(ion_state)
            if ion_number == 1:
                keyword = element
            else:
                keyword = "%s_p%d" %  (element, (ion_number-1))
            field = "%s_number_density" % keyword
        self.field = field

    def __repr__(self):
        return self.identifier

class LineDatabase:
    """
    LineDatabase is the class that holds all of the spectral lines to be used
    by the SpectrumGenerator.
    """
    def __init__(self, input_file=None):
        self.lines_all = []
        self.lines_subset = []
        self.input_file = input_file
        if input_file is not None:
            self.load_line_list_from_file(input_file)
        else:
            self.input_file = 'Manually Entered'

    def add_line(self, element, ion_state, wavelength, gamma,
                 f_value, field=None, identifier=None):
        """
        """
        self.lines_all.append(Line(element, ion_state, wavelength, gamma,
                                   f_value, field, identifier))

    def load_line_list_from_file(self, filename):
        # check to see if file exists in trident/data/line_lists
        # if not, look in cwd
        filename = os.path.join(trident_path(), "data", "line_lists", filename)
        if not os.path.isfile(filename):
            filename = filename.split('/')[-1]
        if not os.path.isfile(filename):
            raise RuntimeError("line_list %s is not found in local "
                               "directory or in trident/data/line_lists "
                               % (filename.split('/')[-1]))

        # Step through each line of text in file and add to database
        for line in file(filename).readlines():
            online = line.rstrip().split()
            if line.startswith("#") or len(online) < 5:
                continue

            element, ion_state, wavelength, gamma, f_value = online[:5]

            # optional identifier should be added if existent
            if len(online) > 5:
                identifier = " ".join(online[5:])
            else:
                identifier = None
            self.add_line(element, ion_state, wavelength, gamma, f_value,
                          identifier=identifier)

    def select_lines(self, element=None, ion_state=None, wavelength=None,
                    identifier=None):
        """
        """
        counter = 0
        for line in self.lines_all:
            # identifier set; use it to find line
            if identifier is not None:
                if line.identifier == identifier:
                    self.lines_subset.append(line)
                    counter += 1
            # element, ion, and wavelength set; use them to find line
            elif ion_state is not None and wavelength is not None:
                if line.element == element and line.ion_state == ion_state \
                   and round(float(line.wavelength), 0) == \
                       round(float(wavelength), 0):
                    self.lines_subset.append(line)
                    counter += 1
            # element and ion set; use them to find line
            elif ion_state is not None:
                if line.element == element and line.ion_state == ion_state:
                    self.lines_subset.append(line)
                    counter += 1
            # only element set; use it to find line
            else:
                if line.element == element:
                    self.lines_subset.append(line)
                    counter += 1
        return counter

    def parse_subset(self, subsets):
        """
        ["C", "C II", "C II 1402", "H I"]
        """
        # if no subsets specified, then use all lines available
        if subsets is None:
            self.lines_subset = self.lines_all
            mylog.info("Using all %d available lines in '%s'." % \
                       (len(self.lines_all), self.input_file))
            return self.lines_subset
        if isinstance(subsets, basestring):
            subsets = [subsets]
        for val in subsets:
            # try to add line based on identifier
            if self.select_lines(identifier=val) > 0:
                continue
            val = val.split()
            if len(val) == 1:
                # add all lines associated with an element
                if self.select_lines(val[0]) == 0:
                    mylog.info("No lines found in subset '%s'." % val[0])
            elif len(val) == 2:
                # add all lines associated with an ion
                if self.select_lines(val[0], val[1]) == 0:
                    mylog.info("No lines found in subset '%s %s'." % \
                               (val[0], val[1]))
            elif len(val) == 3:
                # add only one line
                if self.select_lines(val[0], val[1], val[2]) == 0:
                    mylog.info("No lines found in subset '%s %s %s'." %
                               (val[0], val[1], val[2]))

        # Get rid of duplicates in subset and re-sort
        self.lines_subset = uniquify(self.lines_subset)
        return self.lines_subset
