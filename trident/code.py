import h5py
import numpy as np
import os
from yt.analysis_modules.absorption_spectrum.api import \
      AbsorptionSpectrum
from yt.funcs import mylog
from matplotlib import pyplot
      
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
        AbsorptionSpectrum.__init__(self, lambda_min, lambda_max, n_lambda)

    def _get_qso_spectrum(self, redshift=0.0, filename=None):
        """
        Read in the composite QSO spectrum and return an interpolated version 
        to fit the desired wavelength interval and binning.
        """

        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "..", "data",
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

    def _get_milky_way_foreground(self, filename=None):
        """
        Read in the composite QSO spectrum and return an interpolated version 
        to fit the desired wavelength interval and binning.
        """

        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "..", "data",
                                    "superstack.dat")

        data = np.loadtxt(filename)
        MW_lambda = data[:, 0]
        MW_flux = data[:, 1]

        index = np.digitize(self.lambda_bins, MW_lambda)
        np.clip(index, 1, MW_lambda.size - 1, out=index)
        slope = (MW_flux[index] - MW_flux[index - 1]) / \
          (MW_lambda[index] - MW_lambda[index - 1])
        my_flux = slope * (self.lambda_bins - MW_lambda[index]) + MW_flux[index]
        # just set values that go beyond the data to 1
        my_flux[self.lambda_bins > 1799.9444] = 1.0
        return my_flux

    def add_milky_way_foreground(self, flux_field=None,
                                 filename=None):
        if flux_field is None:
            flux_field = self.flux_field
        MW_spectrum = self._get_milky_way_foreground(filename=filename)
        flux_field *= MW_spectrum

    def add_qso_spectrum(self, flux_field=None,
                         redshift=0.0, filename=None):
        if flux_field is None:
            flux_field = self.flux_field
        qso_spectrum = self._get_qso_spectrum(redshift=redshift,
                                              filename=filename)
        flux_field *= qso_spectrum

    def add_gaussian_noise(self, snr, n_bins=None, out=None, seed=None):
        np.random.seed(seed)
        if n_bins is None:
            n_bins = self.lambda_bins.size
        if out is None:
            out = self.flux_field
        np.add(out, np.random.normal(loc=0.0, scale=1/float(snr), size=n_bins),
               out=out)
        return out

    def apply_lsf(self, instrument, flux_field=None):
        instruments = {"COS":"avg_COS.txt"}
        if flux_field is None:
            flux_field = self.flux_field
        if instrument in instruments:
            lsf_filename = os.path.join(os.path.dirname(__file__), "..", "data",
                                        "lsf_kernels", instruments[instrument])
            lsf_file = open(lsf_filename, 'r')
            lsf_kernel = []
            for line in lsf_file:
                lsf_kernel.append(float(line.split()[1]))
            self.flux_field = np.convolve(lsf_kernel,flux_field,'same')
            lsf_file.close()
        else:
            raise RuntimeError("You have not chosen a valid instrument for your LSF.")

    def load_spectrum(self, filename=None):
        if not filename.endswith(".h5"):
            raise RuntimeError("Only hdf5 format supported for loading spectra.")
        in_file = h5py.File(filename, "r")
        self.lambda_bins = in_file['wavelength'].value
        self.flux_field = in_file['flux'].value
        in_file.close()
    
    def load_line_list(self, filename=None):
        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), "..", "data",
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

            if "Ly" in list_ion:
                field = "H_number_density" 
            else:
                field = "%s_Cloudy_eq_NumberDensity_post" % ion
                
            #field = "%s_Cloudy_eq_NumberDensity_post" % ion
            self.add_line(label, field, float(wavelength),
                          float(f_value), float(gamma),
                          atomic_mass[element], label_threshold=1e3)
        mylog.info("Load %d lines from %s." % (len(self.line_list), filename))

    def _write_spectrum_ascii(self, filename):
        print "Writing spectrum to ascii file: %s." % filename
        f = open(filename, 'w')
        f.write("# wavelength[A] flux rootflux\n")
        for i in xrange(self.lambda_bins.size):
            f.write("%e %e %e\n" % (self.lambda_bins[i],
                                    self.flux_field[i], 
                                    self.flux_field[i]**0.5))
        f.close()

    def _write_spectrum_fits(self, filename):
        """
        Write spectrum to a fits file.
        """
        try:
            import pyfits
        except:
            print "Could not import the pyfits module.  Please install pyfits."
            return

        print "Writing spectrum to fits file: %s." % filename
        col1 = pyfits.Column(name='wavelength', format='E', array=self.lambda_bins)
        col2 = pyfits.Column(name='flux', format='E', array=self.flux_field)
        col3 = pyfits.Column(name='rootflux', format='E', array=self.flux_field**0.5)
        cols = pyfits.ColDefs([col1, col2, col3])
        tbhdu = pyfits.new_table(cols)
        tbhdu.writeto(filename, clobber=True)

    def _write_spectrum_hdf5(self, filename):
        """
        Write spectrum to an hdf5 file.

        """
        print "Writing spectrum to hdf5 file: %s." % filename
        output = h5py.File(filename, 'w')
        output.create_dataset('wavelength', data=self.lambda_bins)
        output.create_dataset('flux', data=self.flux_field)
        output.create_dataset('rootflux', data=self.flux_field**0.5)
        output.close()

def plot_spectrum(wavelength, flux, filename):
    # number of rows and columns
    n_rows = 1
    n_columns = 1

    # blank space between edge of figure and active plot area
    top_buffer = 0.05
    bottom_buffer = 0.15
    left_buffer = 0.05
    right_buffer = 0.03
    
    # blank space between plots
    hor_buffer = 0.05
    vert_buffer = 0.05

    # calculate the height and width of each panel
    panel_width = ((1.0 - left_buffer - right_buffer - 
                    ((n_columns-1)*hor_buffer)) / n_columns)
    panel_height = ((1.0 - top_buffer - bottom_buffer - 
                     ((n_rows-1)*vert_buffer)) / n_rows)

    # create a figure (figsize is in inches)
    pyplot.figure(figsize=(12, 4))

    # loop over all panels
    for i in range(n_rows * n_columns):

        # get the row and column number
        my_row = i / n_columns
        my_column = i % n_columns

        # calculate the position of the bottom, left corner of this plot
        left_side = left_buffer + (my_column * panel_width) + \
          my_column * hor_buffer
        top_side = 1.0 - (top_buffer + (my_row * panel_height) + \
                          my_row * vert_buffer)
        bottom_side = top_side - panel_height

        # create an axes object on which we will make the plot
        my_axes = pyplot.axes((left_side, bottom_side, panel_width, panel_height))

        my_axes.plot(wavelength, flux)
        my_axes.xaxis.set_label_text("$\\lambda$ [$\\AA$]")
        my_axes.yaxis.set_label_text("relative flux")
        my_axes.set_xlim(wavelength.min(), wavelength.max())
        my_axes.set_ylim(0, flux.max()*1.1)
        pyplot.savefig(filename)
