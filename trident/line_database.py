import roman

class Line:
    """An individual atomic transition identified uniquely by element, 
    ionic state, wavelength, gamma, oscillator strength, and identifier.

    Parameters
    ----------

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
        self.ion = ion
        self.wavelength = wavelength
        self.gamma = gamma
        self.f_value = f_value
        self.name = element+ion_state+' %d' % round(wavelength, 0)
        if identifier is None:
            identifier = self.name
        self.identifier = identifier

        # Automatically populate the field if not defined
        if field is None:
            ion_number = roman.fromRoman(ion)
                if ion_number == 1:
                    keyword = element
                else:
                    keyword = "%s_p%d" %  (element, (ion_number-1))
            field = "%s_number_density" % keyword
        self.field = field

class LineDatabase:
    """
    LineDatabase is the class that holds all of the spectral lines to be used
    by the SpectrumGenerator.
    """
    def __init__(input_file=None):
        self.lines = []
        self.input_file = input_file
        if input_file is not None:
            self.load_line_list_from_file(input_file)

    def add_line(self, element, ion_state, wavelength, gamma, 
                 f_value, field=None, identifier=None):
        """
        """
        self.lines.append(Line(element, ion_state, wavelength, gamma, 
                               f_value, field, identifier))

    def load_line_list_from_file(self, filename):
        # check to see if file exists in cwd, if not, check in
        # trident/data/line_lists
        if not os.path.isfile(filename):
            filename = os.path.join(os.path.dirname(__file__), "..",
                                    "data", "line_lists", filename)
        if not os.path.isfile(filename):
            raise RuntimeError("line_list %s is not found in local "
                               "directory or in trident/data/line_lists "
                               % (filename.split('/')[-1]))

        # Step through each line of text in file and add to database
        for line in file(filename).readlines():
            online = line.split()
            if line.startswith("#") or len(online) < 5 or len(online) > 6: 
                continue

            element, ion_state, wavelength, gamma, f_value = online[:5]

            # optional identifier should be added if existent
            try:
                identifier = online[5]
            except:
                identifier = None
            self.add_line(element, ion_state, wavelength, gamma, f_value, 
                          identifier=identifier)
            
    def parse_subset(subsets):
        """
        ["C", "CII", "CII 1402", "HI"]
        """
        subset = []
        for val in subsets:
            element = None
            ion = None
            line = None
            if " " in val:
                val = val.split()
            else:
                val = [val]
            if len(val[0]) == 1:
                element = val
                

        mylog.info("Loading %d lines from %s." % (len(self.line_list), self.input_file))
        return subset
