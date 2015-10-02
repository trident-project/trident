class Line:
    """An individual atomic transition identified uniquely by element, 
    ionic state, wavelength, gamma, oscillator strength, and name.

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

    oscillator strength : float
        The oscillator strength of the transition

    name : string, optional
        An optional name for the transition


    Example: 

    Create a Line object for the neutral hydrogen 1215 Angstroms transition.

    >>> HI = Line('H', 'I', 1215.67, 469860000, 0.41641, 'Lya')
    
    """
    def __init__(self, element, ion_state, wavelength, gamma, oscillator, 
                 alt_name=None):
        self.element = element
        self.ion = ion
        self.wavelength = wavelength
        self.gamma = gamma
        self.oscillator = oscillator
        self.identifier = element+ion_state+'_%d' % round(wavelength, 0)
        self.name = name
