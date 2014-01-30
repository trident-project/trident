import h5py
import numpy as np

from yt.analysis_modules.absorption_spectrum.api import \
      AbsorptionSpectrum
      
atomic_mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941,
               'Be': 9.012182, 'B': 10.811, 'C': 12.0107,
               'N': 14.0067, 'O': 15.9994, 'F': 18.9984032,
               'Ne': 20.1797, 'Na': 22.989770, 'Mg': 24.3050,
               'Al': 26.981538, 'Si': 28.0855, 'P': 30.973761,
               'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
               'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955910,
               'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961,
               'Mn': 54.938049, 'Fe': 55.845, 'Co': 58.933200,
               'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409}

class SpectrumGenerator(AbsorptionSpectrum):
    def __init__(self, lambda_min, lambda_max, n_lambda):
        AbsorptionSpectrum.__init__(lambda_min, lambda_max, n_lambda)

    def _get_qso_spectrum(self, redshift=0.0, filename=None):
        """
        Read in the composite QSO spectrum and return an interpolated version 
        to fit the desired wavelength interval and binning.
        """

        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "data",
                                    "hstrq_sdss.asc")

        data = np.loadtxt(filename)
        qso_lambda = data[:, 0]
        qso_lambda += qso_lambda * redshift
        qso_flux = data[:, 1]

        index = np.digitize(self.lambda_bins, qso_lambda)
        np.clip(index, 1, qso_lambda.size - 1, out=index)
        slope = (qso_flux[index] - qso_flux[index - 1]) / \
          (qso_lambda[index] - qso_lambda[index - 1])
        my_flux = slope * (self.lambda_bins - qso_lambda[index]) + qso_flux[index]
        return my_flux

    def add_qso_spectrum(self, flux_field=None,
                         redshift=0.0, filename=None):
        if flux_field=None:
            flux_field = self.flux_field
        qso_spectrum = get_qso_spectrum(redshift=redshift,
                                        filename=filename)
        flux_field *= qso_spectrum

    def add_gaussian_noise(self, snr, n_bins=None, out=None):
        if n_bins is None:
            n_bins = self.lambda_bins.size
        if out is None:
            out = self.flux_field
        np.add(out, np.random.normal(loc=0.0, scale=1/float(snr), size=n_bins),
               out=out)
        return out

    def load_line_list(filename=None):
        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "data",
                                    "Nist_elem_list.txt")
        for line in file(filename).readlines():
            online = line.split()
            if line.startswith("#") or "--" in line or len(online) != 4: continue
            list_ion, wavelength, gamma, f_value = online
            label = list_ion
            ion = list_ion
            if "Ly" in list_ion:
                ion = "HI"
                label = "HI %s" % list_ion
            if "*" in ion:
                ion = ion[:ion.find("*")]
            element = ion[:2]
            if element[1].isupper():
                element = element[:1]        
            field = "%s_Cloudy_eq_NumberDensity_post" % ion
            sp.add_line(label, field, wavelength,
                        f_value, gamma, atomic_mass[element])
        print "Load %d lines from %s." % (len(self.line_list), filename)

if __name__ == "__main__":
    my_flux = get_qso_spectrum(900, 1800, 900)
    #line_list = load_line_list("short_list.txt")
    line_list = load_line_list()

    # Create spectrum object: lambda_min, lambda_max, number of bins
    lambda_min = 800 # angstroms
    lambda_max = 2000
    n_lambda = 10000
    sp = AbsorptionSpectrum(lambda_min, lambda_max, n_lambda)
    for line in line_list:

    wavelength, flux = sp.make_spectrum("mod_lightray_123456789.h5",
                                        output_file='spectrum.h5', 
                                        line_list_file='lines.txt')

    SNR = 30
    qso_flux = np.ones(n_lambda)
    add_qso_spectrum(wavelength, qso_flux, redshift=0.5)
    add_qso_spectrum(wavelength, flux, redshift=0.5)
    add_gaussian_noise(SNR, n_lambda, out=flux)
    np.clip(flux, 0, flux.max(), out=flux)

    from matplotlib import pyplot
    pyplot.plot(wavelength, flux)
    pyplot.plot(wavelength, qso_flux)
    pyplot.savefig("noise_qso.png")
    #pyplot.show()
