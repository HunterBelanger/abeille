import pyPapillonNDL as pndl
import yaml
import readline

EV_PER_MEV = 1.E6
K_BOLTZMANN = 8.617333262E-5 # eV/K

readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")

nuclides = {}
neutron_dir = {}
tsl_dir = {}

fname = input(" Enter path to xsdir file => ")

mcnp_xsdir = open(fname)

for line in mcnp_xsdir:
  line = line.replace('+', '')
  line = line.replace('ptable','')
  line = line.strip()
  line = line.split()

  # Make sure that the line isn't empty, and that it is for
  # a continuous-energy ace file, not a thermal file !
  if len(line) == 10 and line[0][-1] == 'c':
    ace_fname = line[2]
    temp = float(line[-1]) * EV_PER_MEV / K_BOLTZMANN # in k now
    zaid = int(line[0][:-4])
  
    # Saw ZAIDs for Am242 and Am242m1, due to horid traditions
    if zaid == 95242:
      zaid = 95642
    if zaid == 95642:
      zaid = 95242

    Z = int(zaid / 1000)
    A = int(zaid - Z*1000)

    nuclide = pndl.Nuclide(pndl.ZAID(Z,A))
    symb = nuclide.symbol()

    if not symb in neutron_dir:
      neutron_dir[symb] = []
    neutron_dir[symb].append({"file": ace_fname, "temperature": round(temp*10.)/10.})

  elif len(line) == 10 and line[0][-1] == 't':
    ace_fname = line[2]
    temp = float(line[-1]) * EV_PER_MEV / K_BOLTZMANN # in k now
    tsl_name = line[0][:-4]

    # See if the tsl name is in the tsls yet
    if not tsl_name in tsl_dir:
      tsl_dir[tsl_name] = []
    # Save tsl for latter. Will construct special nuclides with them
    tsl_dir[tsl_name].append({"file": ace_fname, "temperature": round(temp*10.)/10.})
mcnp_xsdir.close()


# Must now add special nuclides which have a TSL
lib80x_name_converter = {"h-h2o": ["H1_H2O"], "d-d2o": ["H2_D2O"],
                         "o-d2o": ["O16_D2O", "O17_D2O", "O18_D2O"],
                         "zr-zrh": ["Zr90_ZrH", "Zr91_ZrH", "Zr92_ZrH",
                                    "Zr94_ZrH", "Zr96_ZrH"],
                         "h-zrh": ["H1_ZrH"],
                         "u-uo2": ["U238_UO2", "U235_UO2", "U234_UO2"],
                         "o-uo2": ["O16_UO2", "O17_UO2", "O18_UO2"]
                         }

tsl_to_ce = {"H1_H2O": "H1", "H2_D2O": "H2", "O16_D2O": "O16",
             "O17_D2O": "O17", "O18_D2O": "O18", "Zr90_ZrH": "Zr90",
             "Zr91_ZrH": "Zr91", "Zr92_ZrH": "Zr92", "Zr94_ZrH": "Zr94",
             "Zr96_ZrH": "Zr96", "H1_ZrH": "H1", "U238_UO2": "U238",
             "U235_UO2": "U235", "U234_UO2": "U234", "O16_UO2": "O16",
             "O17_UO2": "O17", "O18_UO2": "O18"}

# Go through all neutron_dir entries, and add a nuclides entry for it
for nuc_key in neutron_dir:
  temps = []
  for entry in neutron_dir[nuc_key]:
    temps.append(entry['temperature'])
  temps.sort()
  nuclides[nuc_key] = {"neutron": nuc_key, "temperatures": temps}

# We go through all tsls which we found
for tsl in tsl_dir:
  # We now chech if it is in our name coverter
  if tsl in lib80x_name_converter:
    # We now go through all possible isotopes for this tsl
    names = lib80x_name_converter[tsl]
    for name in names:
      # Get the ce data for this nuclide
      ce_name = tsl_to_ce[name]

      # Now we go through all temperatures
      for tsl_temp in tsl_dir[tsl]:
        temp = tsl_temp["temperature"]
        tsl_file = tsl_temp["file"]

        # See if this temp exists for the ce data
        ce_file = ""
        for t in neutron_dir[ce_name]:
          if abs(t["temperature"] - temp) < 3.: # Use any eval which is within 3 degrees
            ce_file = t["file"]
            temp = t["temperature"]
            break

        if len(ce_file) > 0:
          if not name in nuclides:
            nuclides[name] = {"neutron": ce_name, "tsl": tsl, "temperatures": []}

          nuclides[name]["temperatures"].append(temp)

xsdir = {"nuclides": nuclides, "neutron-dir": neutron_dir, "tsl-dir": tsl_dir}

with open('xsdir.yaml', 'w') as outfile:
    yaml.dump(xsdir, outfile, default_flow_style=False)
