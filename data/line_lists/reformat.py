from collections import defaultdict
import numpy as np
import roman

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

lines = file("all_lines.txt", "r").readlines()

all_lines = defaultdict(list)
for line in lines:
    line = line.strip()
    online = line.split()
    ion = online[0]
    if ion[1].islower():
        element = ion[:2]
        ion_state = ion[2:]
    else:
        element = ion[:1]
        ion_state = ion[1:]
    if ion_state.endswith("*"):
        ion_state = ion_state[:-1]
        altname = "%s %s* %d" % (element, ion_state, np.round(float(online[1])))
    else:
        altname = ""

    all_lines[element].append({"element": element,
                               "ion": roman.fromRoman(ion_state),
                               "wavelength": online[1],
                               "freq": (1 / float(online[1])),
                               "gamma": online[2],
                               "f_value": online[3],
                               "altname": altname})

elements = sorted(all_lines.keys(),
                  key=lambda x: atomic_mass[x])

f = file("final_lines.txt", "w")
f.write("#%-9s%16s%16s%16s%16s\n" % ("Ion", "Wavelength [A]", "gamma", "f_value", "alt. name"))
ly_no = 1
for element in elements:
    ion_lines = sorted(all_lines[element], key=lambda x: (x["ion"], x["freq"]))
    for ion_line in ion_lines:
        if element == "H":
            ion_line["altname"] = "Ly %d" % ly_no
            ly_no += 1
        output_ion = "%s %s" % (ion_line["element"], roman.toRoman(ion_line["ion"]))
        f.write("%-10s%16f%16e%16e%16s\n" % (output_ion, float(ion_line["wavelength"]),
                                             float(ion_line["gamma"]), float(ion_line["f_value"]),
                                             ion_line["altname"]))
f.close()
