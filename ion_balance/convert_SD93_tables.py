import SD93IonBalance as IB
import numpy as na
import string
import h5py
import sys

ascii_files = sys.argv[1:]
output_file = 'tables/ion_balance.h5'
output = h5py.File(output_file,'w')

for file in ascii_files:
    prefix = file.find("_ion_balance.txt")
    if prefix > 0:
        i_table = IB.SD93IonBalanceTable(file)
        element = file[prefix-2:prefix]
        if not element[0].isalpha(): element = element[1]
        element = string.capitalize(element)
        print "Add %s to %s as %s." % (file, output_file, element)
        output.create_dataset(element, data=i_table.ion_fraction, dtype='>f4')
        output[element].attrs['Temperature'] = i_table.temperature.astype('>f4')
        del i_table
    else:
        print "Skipping %s." % file

output.close()
