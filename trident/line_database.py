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
    from_roman

def uniquify(list):
   # order preserving method for reducing duplicates in a list
   checked = []
   for val in list:
       if val not in checked:
           checked.append(val)
   return checked

class Line:
    """A class representing an individual atomic transition.  Each Line object
    is uniquely identified by element, ionic state, wavelength, gamma, 
    oscillator strength, and identifier.

    **Parameters**

    :element: string

        The element of the transition using element's symbol on periodic table
        Example: 'H', 'C', 'Mg'

    :ion_state: string

        The roman numeral representing the ionic state of the transition
        Example: 'I' for neutral species, 'II' for singly ionized, etc.

    :wavelength: float
    
        The wavelength of the transition in angstroms
        Example: 1216 for Lyman alpha

    :gamma: float

        The gamma of the transition in Hertz

    :f_value: float

        The oscillator strength of the transition

    :field: string, optional

        The default yt field name associated with the ion responsible for
        this line
        Example: 'H_p1_number_density' for HII

    :identifier: string, optional

        An optional identifier for the transition
        Example: 'Ly a' for Lyman alpha

    **Example**

    Create a Line object for the neutral hydrogen 1215 Angstroms transition.

    >>> HI = Line('H', 'I', 1215.67, 469860000, 0.41641, 'Ly a')

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
            ion_number = from_roman(ion_state)
            if ion_number == 1:
                keyword = "%s" %  (element)
            else:
                keyword = "%s_p%d" %  (element, (ion_number-1))
            field = "%s_number_density" % keyword
        self.field = field

    def __repr__(self):
        return self.identifier

class LineDatabase:
    """
    Class for storing and selecting collections of spectral lines.  These lines
    will be used in the :class:`~trident.SpectrumGenerator` and 
    :class:`~trident.add_ion_fields()` functionality.

    Without arguments, the LineDatabase will be empty, and you must
    manually add individual lines to it using the 
    :class:`~trident.LineDatabase.add_line` function.
    If LineDatabase is provided with an optional :input_file:, it will 
    automatically add spectral lines for each corresponding line in the list.

    Once created, you can select a subset of the total lines present in 
    the database for further use.  Use the 
    :class:`~trident.LineDatabase.parse_subset`
    function to accomplish this.

    **Parameters**

    :input_file: string, optional

        An optional input_file can be provided to pre-store a list of Line
        objects.  input_file should be a tab delimited text file of the 
        format:

        element, ion_state, wavelength, gamma, f_value, (name)

        H,         I,        1215.67,  4.69e8, 4.16e-1, Ly a

    **Example**

    >>> # Create a LineDatabase using the lines present in lines.txt
    >>> ldb = LineDatabase('lines.txt')

    >>> # Parse ldb and only select Lyman alpha, Mg II and Fe lines
    >>> lines = ldb.parse_subset(lines=['H I 1216', 'Mg II', 'Fe'])
    >>> print lines

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
        Manually add a line to the :class:`~trident.LineDatabase`.

        **Parameters**

        :element: string

            The element of the transition using element's symbol on periodic table
            Example: 'H', 'C', 'Mg'
    
        :ion_state: string

            The roman numeral representing the ionic state of the transition
            Example: 'I' for neutral species, 'II' for singly ionized, etc.
    
        :wavelength: float

            The wavelength of the transition in angstroms
            Example: 1216 for Lyman alpha

        :gamma: float

            The gamma of the transition in Hertz

        :f_value: float

            The oscillator strength of the transition

        :field: string, optional

            The default yt field name associated with the ion responsible for
            this line
            Example: 'H_p1_number_density' for HII

        :identifier: string, optional

            An optional identifier for the transition
            Example: 'Ly a' for Lyman alpha


        **Example**


        >>> # Create a LineDatabase using the lines present in lines.txt
        >>> ldb = LineDatabase('lines.txt')

        >>> # Manually add the neutral hydrogen line to ldb
        >>> ldb.add_line('H', 'I', 1215.67, 469860000, 0.41641, 'Ly a')
        >>> print ldb.lines_all
        """
        self.lines_all.append(Line(element, ion_state, wavelength, gamma,
                                   f_value, field, identifier))

    def load_line_list_from_file(self, filename):
        """
        Load a line list from a file into the LineDatabase.
        Line list file is a tab-delimited text file in the format:

        element, ion_state, wavelength, gamma, f_value, (name)

        H,         I,        1215.67,  4.69e8, 4.16e-1, Ly a

        **Parameters**

        filename : string

            The filename of the list to add.  First looks in 
            trident.__path__/data/line_lists directory, then in cwd.
        """
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
        Select lines based on atom, ion state, identifier, and/or wavelength.
        Once you've created a LineDatabase, you can subselect certain lines 
        from it based on line characteristics.  Recommended to use 
        :class:`~trident.LineDatabase.parse_subset` instead which allows 
        selecting of multiple sets of lines simultaneously.

        **Parameters**
        
        :element: string, optional

            The element of the transition using element's symbol on periodic table
            Example: 'H', 'C', 'Mg'

        :ion_state: string, optional
        
            The roman numeral representing the ionic state of the transition
            Example: 'I' for neutral species, 'II' for singly ionized, etc.

        :wavelength: float, optional

            The wavelength of the transition in angstroms
            Example: 1216 for Lyman alpha

        :identifier: string, optional

            An optional identifier for the transition
            Example: 'Ly a' for Lyman alpha

        **Returns**
            
        :counter: int 

            A counter indicating how many lines were selected and added to 
            the line_subset list attribute.

        **Example**
        
        >>> ldb = LineDatabase('lines.txt')
        >>> ldb.select_lines(element='Mg', ion_state='II')
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

    def parse_subset(self, subsets=None):
        """
        Select multiple lines based on atom, ion state, identifier, and/or 
        wavelength.  Once you've created a LineDatabase, you can subselect 
        certain lines from it based on line characteristics.  Preferred to 
        use this method over :class:`~trident.LineDatabase.select_lines`.

        Will return the unique union of all lines matching the specified 
        subsets from the :class:`~trident.LineDatabase`.

        **Parameters**
        
        :subsets: list of strings, optional

            List strings matching possible lines.  Strings can be of the
            form:
            * Atom - Examples: "H", "C", "Mg"
            * Ion - Examples: "H I", "H II", "C IV", "Mg II"
            * Line - Examples: "H I 1216", "C II 1336", "Mg II 1240"
            * Identifier - Examples: "Ly a", "Ly b"

            If set to None, selects **all** lines in 
            :class:`~trident.LineDatabase`.
            Default: None

        **Returns**
            
        :line subset: list of :class:`trident.Line` objects

            A list of the Lines that were selected

        **Example**
        
        >>> # Get a list of all lines of Carbon, Mg II and Lyman alpha
        >>> ldb = LineDatabase('lines.txt')
        >>> lines = ldb.parse_subset(['C', 'Mg II', 'H I 1216'])
        >>> print lines
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

    def parse_subset_to_ions(self, subsets=None):
        """
        Select ions based on those needed to create specific lines.  
        Once you've created a LineDatabase, you can subselect 
        certain ions from it based on the line characteristics of atom,
        ion state, identifier, and/or wavelength.  Similar to 
        :class:`~trident.LineDatabase.parse_subset` but outputs a list of
        ion tuples (e.g. ('H', 1), ('Fe', 2)), instead of a list of 
        :class:`~trident.Line` objects.

        Will return the unique union of all ions matching the specified 
        subsets from the :class:`~trident.LineDatabase`.

        **Parameters**
        
        :subsets: list of strings, optional

            List strings matching possible lines.  Strings can be of the
            form:
            * Atom - Examples: "H", "C", "Mg"
            * Ion - Examples: "H I", "H II", "C IV", "Mg II"
            * Line - Examples: "H I 1216", "C II 1336", "Mg II 1240"
            * Identifier - Examples: "Ly a", "Ly b"
            
            If set to None, selects ions necessary to produce **all** lines
            in :class:`~trident.LineDatabase`.
            Default: None

        **Returns**
            
        :ion subset: list of ion tuples

            A list of the ions necessary to produce the desired lines
            Each ion tuple is of the form ('H', 1) = neutral hydrogen

        **Example**
        
        Get a list of all ions necessary to generate lines for Carbon, 
        Mg II and Lyman alpha

        >>> ldb = LineDatabase('lines.txt')
        >>> ions = ldb.parse_subset_to_ions(['C', 'Mg II', 'H I 1216'])
        >>> print ions
        """
        self.parse_subset(subsets)
        ions = []
        for line in self.lines_subset:
            ions.append((line.element, from_roman(line.ion_state)))
        ions = uniquify(ions)
        return ions
