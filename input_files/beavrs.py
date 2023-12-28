import yaml
import numpy as np

# Control Rod Step
# In range [0, 574], 0 is all rods out, 574 is all rods full in.
# Use S = 15 for HZP
S = 15

dS = 1.58193

# Steps for each bank
SD = max(0, 228-S)

if SD < 113:
  SC = max(0, 228-S+113+3)
else:
  SC = 228

if SC < 113:
  SB = max(0, 228-S+113*2+5)
else:
  SB = 228

if SB < 113:
  SA = max(0, 228-S+113*3+7)
else:
  SA = 228

# Shutdown banks are not coupled. Here, 0 is full out, 228 is full in
Sa = 0
Sb = 0
Sc = 0
Sd = 0
Se = 0

AZC = float(SA)*dS # Control Rod Bank A Z Correction
BZC = float(SB)*dS # Control Rod Bank B Z Correction
CZC = float(SC)*dS # Control Rod Bank C Z Correction
DZC = float(SD)*dS # Control Rod Bank D Z Correction

aZC = float(228 - Sa)*dS # Control Rod Shutdown Bank A Z Correction
bZC = float(228 - Sb)*dS # Control Rod Shutdown Bank B Z Correction
cZC = float(228 - Sc)*dS # Control Rod Shutdown Bank C Z Correction
dZC = float(228 - Sd)*dS # Control Rod Shutdown Bank D Z Correction
eZC = float(228 - Se)*dS # Control Rod Shutdown Bank E Z Correction


# Lists to hold all info
materials = []
surfaces = []
cells = []
universes = []

#===============================================================================
# Materials

# Fuel 1.6% Enriched
F16 = 1
materials.append({"name": "Fuel-1.6%",
                  "id": 1,
                  "temperature": 800.,
                  "density-units": "sum",
                  "fractions": "atoms", 
                  "composition": [{"nuclide": "O16",  "fraction": 4.5897e-02},
                                  {"nuclide": "O17",  "fraction": 1.7436e-05},
                                  {"nuclide": "O18",  "fraction": 9.2032e-05},
                                  {"nuclide": "U234", "fraction": 3.0131e-06},
                                  {"nuclide": "U235", "fraction": 3.7503e-04},
                                  {"nuclide": "U238", "fraction": 2.2625e-02}
                                 ]})

# Fuel 2.4% Enriched
F24 = 2
materials.append({"name": "Fuel-2.4%",
                  "id": 2,
                  "temperature": 800.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "O16",  "fraction": 4.5830e-02},
                                  {"nuclide": "O17",  "fraction": 1.7411e-05},
                                  {"nuclide": "O18",  "fraction": 9.1898e-05},
                                  {"nuclide": "U234", "fraction": 4.4842e-06},
                                  {"nuclide": "U235", "fraction": 5.5814e-04},
                                  {"nuclide": "U238", "fraction": 2.2407e-02}
                                 ]})

# Fuel 3.1% Enriched
F31 = 3
materials.append({"name": "Fuel-3.1%",
                  "id": 3,
                  "temperature": 800.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "O16",  "fraction": 4.5853e-02},
                                  {"nuclide": "O17",  "fraction": 1.7420e-05},
                                  {"nuclide": "O18",  "fraction": 9.1942e-05},
                                  {"nuclide": "U234", "fraction": 5.7987e-06},
                                  {"nuclide": "U235", "fraction": 7.2175e-04},
                                  {"nuclide": "U238", "fraction": 2.2253e-02}
                                 ]})

# Fuel 3.2% Enriched
F32 = 4
materials.append({"name": "Fuel-3.2%",
                  "id": 4,
                  "temperature": 800.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "O16",  "fraction": 4.6029e-02},
                                  {"nuclide": "O17",  "fraction": 1.7487e-05},
                                  {"nuclide": "O18",  "fraction": 9.2296e-05},
                                  {"nuclide": "U234", "fraction": 5.9959e-06},
                                  {"nuclide": "U235", "fraction": 7.4630e-04},
                                  {"nuclide": "U238", "fraction": 2.2317e-02}
                                 ]})

# Fuel 3.4% Enriched
F34 = 5
materials.append({"name": "Fuel-3.4%",
                  "id": 5,
                  "temperature": 800.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "O16",  "fraction": 4.6110e-02},
                                  {"nuclide": "O17",  "fraction": 1.7517e-05},
                                  {"nuclide": "O18",  "fraction": 9.2459e-05},
                                  {"nuclide": "U234", "fraction": 6.4018e-06},
                                  {"nuclide": "U235", "fraction": 7.9681e-04},
                                  {"nuclide": "U238", "fraction": 2.2307e-02}
                                 ]})

# Air
Air = 6
materials.append({"name": "Air",
                  "id": 6,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Ar36", "fraction": 7.8730e-09},
                                  {"nuclide": "Ar38", "fraction": 1.4844e-09},
                                  {"nuclide": "Ar40", "fraction": 2.3506e-06},
                                  {"nuclide": "C12",  "fraction": 6.7539e-08},
                                  {"nuclide": "C13",  "fraction": 7.5658e-10},
                                  {"nuclide": "N14",  "fraction": 1.9680e-04},
                                  {"nuclide": "N15",  "fraction": 7.2354e-07},
                                  {"nuclide": "O16",  "fraction": 5.2866e-05},
                                  {"nuclide": "O17",  "fraction": 2.0084e-08},
                                  {"nuclide": "O18",  "fraction": 1.0601e-07}
                                 ]})

# Borosilicate Glass
BSiGlass = 7
materials.append({"name": "Borosilicate Glass",
                  "id": 7,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Al27", "fraction": 1.7352e-03},
                                  {"nuclide": "B10",  "fraction": 9.6506e-04},
                                  {"nuclide": "B11",  "fraction": 3.9189e-03},
                                  {"nuclide": "O16",  "fraction": 4.6514e-02},
                                  {"nuclide": "O17",  "fraction": 1.7671e-05},
                                  {"nuclide": "O18",  "fraction": 9.3268e-05},
                                  {"nuclide": "Si28", "fraction": 1.6926e-02},
                                  {"nuclide": "Si29", "fraction": 8.5944e-04},
                                  {"nuclide": "Si30", "fraction": 5.6654e-04}
                                 ]})

# Ag-In-Cd Control Rods
AgInCd = 8
materials.append({"name": "Ag-In-Cd Control Rods",
                  "id": 8,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Ag107", "fraction": 2.3523e-02},
                                  {"nuclide": "Ag109", "fraction": 2.1854e-02},
                                  {"nuclide": "Cd106", "fraction": 3.3882e-05},
                                  {"nuclide": "Cd108", "fraction": 2.4166e-05},
                                  {"nuclide": "Cd110", "fraction": 3.3936e-04},
                                  {"nuclide": "Cd111", "fraction": 3.4821e-04},
                                  {"nuclide": "Cd112", "fraction": 6.5611e-04},
                                  {"nuclide": "Cd113", "fraction": 3.3275e-04},
                                  {"nuclide": "Cd114", "fraction": 7.8252e-04},
                                  {"nuclide": "Cd116", "fraction": 2.0443e-04},
                                  {"nuclide": "In113", "fraction": 3.4219e-04},
                                  {"nuclide": "In115", "fraction": 7.6511e-03}
                                 ]})

# B4C Control Rods
B4C = 9
materials.append({"name": "B4C Control Rods",
                  "id": 9,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "B10", "fraction": 1.5206e-02},
                                  {"nuclide": "B11", "fraction": 6.1514e-02},
                                  {"nuclide": "C12", "fraction": 1.8972e-02},
                                  {"nuclide": "C13", "fraction": 2.1252e-04}
                                 ]})

# He
He = 10
materials.append({"name": "He",
                  "id": 10,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "He3", "fraction": 4.8089e-10},
                                  {"nuclide": "He4", "fraction": 2.4044e-04}
                                 ]})

# Inconcel 718 
Inconel = 11
materials.append({"name": "Inconel 718",
                  "id": 11,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Cr50", "fraction": 7.8239e-04},
                                  {"nuclide": "Cr52", "fraction": 1.5088e-02},
                                  {"nuclide": "Cr53", "fraction": 1.7108e-03},
                                  {"nuclide": "Cr54", "fraction": 4.2586e-04},
                                  {"nuclide": "Fe54", "fraction": 1.4797e-03},
                                  {"nuclide": "Fe56", "fraction": 2.3229e-02},
                                  {"nuclide": "Fe57", "fraction": 5.3645e-04},
                                  {"nuclide": "Fe58", "fraction": 7.1392e-05},
                                  {"nuclide": "Mn55", "fraction": 7.8201e-04},
                                  {"nuclide": "Ni58", "fraction": 2.9320e-02},
                                  {"nuclide": "Ni60", "fraction": 1.1294e-02},
                                  {"nuclide": "Ni61", "fraction": 4.9094e-04},
                                  {"nuclide": "Ni62", "fraction": 1.5653e-03},
                                  {"nuclide": "Ni64", "fraction": 3.9864e-04},
                                  {"nuclide": "Si28", "fraction": 5.6757e-04},
                                  {"nuclide": "Si29", "fraction": 2.8820e-05},
                                  {"nuclide": "Si30", "fraction": 1.8998e-05}
                                 ]})

# SS304
SS304 = 12
materials.append({"name": "Stainless Steel 304",
                  "id": 12,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Cr50", "fraction": 7.6778e-04},
                                  {"nuclide": "Cr52", "fraction": 1.4806e-02},
                                  {"nuclide": "Cr53", "fraction": 1.6789e-03},
                                  {"nuclide": "Cr54", "fraction": 4.1791e-04},
                                  {"nuclide": "Fe54", "fraction": 3.4620e-03},
                                  {"nuclide": "Fe56", "fraction": 5.4345e-02},
                                  {"nuclide": "Fe57", "fraction": 1.2551e-03},
                                  {"nuclide": "Fe58", "fraction": 1.6703e-04},
                                  {"nuclide": "Mn55", "fraction": 1.7604e-03},
                                  {"nuclide": "Ni58", "fraction": 5.6089e-03},
                                  {"nuclide": "Ni60", "fraction": 2.1605e-03},
                                  {"nuclide": "Ni61", "fraction": 9.3917e-05},
                                  {"nuclide": "Ni62", "fraction": 2.9945e-04},
                                  {"nuclide": "Ni64", "fraction": 7.6261e-05},
                                  {"nuclide": "Si28", "fraction": 9.5281e-04},
                                  {"nuclide": "Si29", "fraction": 4.8381e-05},
                                  {"nuclide": "Si30", "fraction": 3.1893e-05}
                                 ]})

# Zircaloy 4
Zirc4 = 13
materials.append({"name": "Zircaloy 4",
                  "id": 13,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Cr50",  "fraction": 3.2962e-06},
                                  {"nuclide": "Cr52",  "fraction": 6.3564e-05},
                                  {"nuclide": "Cr53",  "fraction": 7.2076e-06},
                                  {"nuclide": "Cr54",  "fraction": 1.7941e-06},
                                  {"nuclide": "Fe54",  "fraction": 8.6698e-06},
                                  {"nuclide": "Fe56",  "fraction": 1.3610e-04},
                                  {"nuclide": "Fe57",  "fraction": 3.1431e-06},
                                  {"nuclide": "Fe58",  "fraction": 4.1829e-07},
                                  {"nuclide": "O16",   "fraction": 3.0744e-04},
                                  {"nuclide": "O17",   "fraction": 1.1680e-07},
                                  {"nuclide": "O18",   "fraction": 6.1648e-07},
                                  {"nuclide": "Sn112", "fraction": 4.6735e-06},
                                  {"nuclide": "Sn114", "fraction": 3.1799e-06},
                                  {"nuclide": "Sn115", "fraction": 1.6381e-06},
                                  {"nuclide": "Sn116", "fraction": 7.0055e-05},
                                  {"nuclide": "Sn117", "fraction": 3.7003e-05},
                                  {"nuclide": "Sn118", "fraction": 1.1669e-04},
                                  {"nuclide": "Sn119", "fraction": 4.1387e-05},
                                  {"nuclide": "Sn120", "fraction": 1.5697e-04},
                                  {"nuclide": "Sn122", "fraction": 2.2308e-05},
                                  {"nuclide": "Sn124", "fraction": 2.7897e-05},
                                  {"nuclide": "Zr90",  "fraction": 2.1828e-02},
                                  {"nuclide": "Zr91",  "fraction": 4.7601e-03},
                                  {"nuclide": "Zr92",  "fraction": 7.2759e-03},
                                  {"nuclide": "Zr94",  "fraction": 7.3734e-03},
                                  {"nuclide": "Zr96",  "fraction": 1.1879e-03}
                                 ]})

# Borated Water
BH2O = 14
RATIO_LW = (4.9456e-02)/((4.9456e-02)+(7.7035e-06))
RATIO_HW = 1. - RATIO_LW
materials.append({"name": "Borated Water",
                  "id": 14,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "B10",     "fraction": 7.9714e-06},
                                  {"nuclide": "B11",     "fraction": 3.2247e-05},
                                  {"nuclide": "H1_H2O",  "fraction": 4.9456e-02},
                                  {"nuclide": "O16",     "fraction": RATIO_LW*2.4673e-02},
                                  {"nuclide": "O17",     "fraction": RATIO_LW*9.3734e-06},
                                  {"nuclide": "O18",     "fraction": RATIO_LW*4.9474e-05},
                                  {"nuclide": "H2_D2O",  "fraction": 7.7035e-06},
                                  {"nuclide": "O16_D2O", "fraction": RATIO_HW*2.4673e-02},
                                  {"nuclide": "O17_D2O", "fraction": RATIO_HW*9.3734e-06},
                                  {"nuclide": "O18_D2O", "fraction": RATIO_HW*4.9474e-05}
                                 ]})

# Nozzle / Support Plate Borated Water
NS_BH2O = 15
RATIO_LW = (6.5512e-02)/((6.5512e-02)+(1.0204e-05))
RATIO_HW = 1. - RATIO_LW
materials.append({"name": "Nozzle/Support Plate Borated Water",
                  "id": 15,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "B10",     "fraction": 1.0559e-05},
                                  {"nuclide": "B11",     "fraction": 4.2716e-05},
                                  {"nuclide": "H1_H2O",  "fraction": 6.5512e-02},
                                  {"nuclide": "O16",     "fraction": RATIO_LW*3.2683e-02},
                                  {"nuclide": "O17",     "fraction": RATIO_LW*1.2416e-05},
                                  {"nuclide": "O18",     "fraction": RATIO_LW*6.5535e-05},
                                  {"nuclide": "H2_D2O",  "fraction": 1.0204e-05},
                                  {"nuclide": "O16_D2O", "fraction": RATIO_HW*3.2683e-02},
                                  {"nuclide": "O17_D2O", "fraction": RATIO_HW*1.2416e-05},
                                  {"nuclide": "O18_D2O", "fraction": RATIO_HW*6.5535e-05}
                                 ]})

# Nozzle / Support Plate Stainless Steel
NS_SS = 16
materials.append({"name": "Nozzle/Support Splate Stainless Steel",
                  "id": 16,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Cr50", "fraction": 3.5223e-04},
                                  {"nuclide": "Cr52", "fraction": 6.7924e-03},
                                  {"nuclide": "Cr53", "fraction": 7.7020e-04},
                                  {"nuclide": "Cr54", "fraction": 1.9172e-04},
                                  {"nuclide": "Fe54", "fraction": 1.5882e-03},
                                  {"nuclide": "Fe56", "fraction": 2.4931e-02},
                                  {"nuclide": "Fe57", "fraction": 5.7578e-04},
                                  {"nuclide": "Fe58", "fraction": 7.6625e-05},
                                  {"nuclide": "Mn55", "fraction": 8.0762e-04},
                                  {"nuclide": "Ni58", "fraction": 2.5731e-03},
                                  {"nuclide": "Ni60", "fraction": 9.9117e-04},
                                  {"nuclide": "Ni61", "fraction": 4.3085e-05},
                                  {"nuclide": "Ni62", "fraction": 1.3738e-04},
                                  {"nuclide": "Ni64", "fraction": 3.4985e-05},
                                  {"nuclide": "Si28", "fraction": 4.3711e-04},
                                  {"nuclide": "Si29", "fraction": 2.2195e-05},
                                  {"nuclide": "Si30", "fraction": 1.4631e-05}
                                 ]})

# Carbon Steel
CS = 17
materials.append({"name": "Nozzle/Support Splate Stainless Steel",
                  "id": 17,
                  "temperature": 575.,
                  "density-units": "sum",
                  "fractions": "atoms",
                  "composition": [{"nuclide": "Al27",  "fraction": 4.3523e-05},
                                  {"nuclide": "B10",   "fraction": 2.5833e-06},
                                  {"nuclide": "B11",   "fraction": 1.0450e-05},
                                  {"nuclide": "C12",   "fraction": 1.0442e-03},
                                  {"nuclide": "C13",   "fraction": 1.1697e-05},
                                  {"nuclide": "Ca40",  "fraction": 1.7043e-05},
                                  {"nuclide": "Ca42",  "fraction": 1.1375e-07},
                                  {"nuclide": "Ca43",  "fraction": 2.3734e-08},
                                  {"nuclide": "Ca44",  "fraction": 3.6673e-07},
                                  {"nuclide": "Ca46",  "fraction": 7.0322e-10},
                                  {"nuclide": "Ca48",  "fraction": 3.2875e-08},
                                  {"nuclide": "Cr50",  "fraction": 1.3738e-05},
                                  {"nuclide": "Cr52",  "fraction": 2.6493e-04},
                                  {"nuclide": "Cr53",  "fraction": 3.0041e-05},
                                  {"nuclide": "Cr54",  "fraction": 7.4778e-06},
                                  {"nuclide": "Cu63",  "fraction": 1.0223e-04},
                                  {"nuclide": "Cu65",  "fraction": 4.5608e-05},
                                  {"nuclide": "Fe54",  "fraction": 4.7437e-03},
                                  {"nuclide": "Fe56",  "fraction": 7.4465e-02},
                                  {"nuclide": "Fe57",  "fraction": 1.7197e-03},
                                  {"nuclide": "Fe58",  "fraction": 2.2886e-04},
                                  {"nuclide": "Mn55",  "fraction": 6.4126e-04},
                                  {"nuclide": "Mo100", "fraction": 2.9814e-05},
                                  {"nuclide": "Mo92",  "fraction": 4.4822e-05},
                                  {"nuclide": "Mo94",  "fraction": 2.8110e-05},
                                  {"nuclide": "Mo95",  "fraction": 4.8567e-05},
                                  {"nuclide": "Mo96",  "fraction": 5.1015e-05},
                                  {"nuclide": "Mo97",  "fraction": 2.9319e-05},
                                  {"nuclide": "Mo98",  "fraction": 7.4327e-05},
                                  {"nuclide": "Nb93",  "fraction": 5.0559e-06},
                                  {"nuclide": "Ni58",  "fraction": 4.0862e-04},
                                  {"nuclide": "Ni60",  "fraction": 1.5740e-04},
                                  {"nuclide": "Ni61",  "fraction": 6.8420e-06},
                                  {"nuclide": "Ni62",  "fraction": 2.1815e-05},
                                  {"nuclide": "Ni64",  "fraction": 5.5557e-06},
                                  {"nuclide": "P31",   "fraction": 3.7913e-05},
                                  {"nuclide": "S32",   "fraction": 3.4808e-05},
                                  {"nuclide": "S33",   "fraction": 2.7420e-07},
                                  {"nuclide": "S34",   "fraction": 1.5368e-06},
                                  {"nuclide": "S36",   "fraction": 5.3398e-09},
                                  {"nuclide": "Si28",  "fraction": 6.1702e-04},
                                  {"nuclide": "Si29",  "fraction": 3.1330e-05},
                                  {"nuclide": "Si30",  "fraction": 2.0653e-05},
                                  {"nuclide": "Ti46",  "fraction": 1.2144e-06},
                                  {"nuclide": "Ti47",  "fraction": 1.0952e-06},
                                  {"nuclide": "Ti48",  "fraction": 1.0851e-05},
                                  {"nuclide": "Ti49",  "fraction": 7.9634e-07},
                                  {"nuclide": "Ti50",  "fraction": 7.6249e-07},
                                  {"nuclide": "V50",   "fraction": 1.1526e-07},
                                  {"nuclide": "V51",   "fraction": 4.5989e-05}
                                 ]})

#===============================================================================
# Geometry

# ALL AXIAL SURFACES
LZC = -230. # Lattice Z Correction
surfaces.append({"id": 96, "type": "zplane", "z0": 460.0000+LZC}) # Highest Extent
surfaces.append({"id": 95, "type": "zplane", "z0": 431.8760+LZC}) # Top of Upper Nozzle
surfaces.append({"id": 94, "type": "zplane", "z0": 423.0490+LZC}) # Bottom of Upper Nozzle
surfaces.append({"id": 93, "type": "zplane", "z0": 421.5320+LZC}) # Top of BPRA Rod Plenum
surfaces.append({"id": 92, "type": "zplane", "z0": 419.7040+LZC}) # Top of Fuel Rod
surfaces.append({"id": 91, "type": "zplane", "z0": 417.1640+LZC}) # Top of Fuel Rod Plenum
surfaces.append({"id": 90, "type": "zplane", "z0": 415.5580+LZC}) # Top of Control Rod Plenum
surfaces.append({"id": 81, "type": "zplane", "z0": 415.1640+LZC}) # Grid 8 Top
surfaces.append({"id": 80, "type": "zplane", "z0": 411.8060+LZC}) # Grid 8 Bottom
surfaces.append({"id": 76, "type": "zplane", "z0": 403.7780+LZC}) # Bottom of Control Rod Plenum 
surfaces.append({"id": 75, "type": "zplane", "z0": 402.5080+LZC}) # Top of Active Fuel
surfaces.append({"id": 74, "type": "zplane", "z0": 401.2380+LZC}) # Top of Active Absorber
surfaces.append({"id": 73, "type": "zplane", "z0": 400.6380+LZC}) # Control Rod Step 228 
surfaces.append({"id": 71, "type": "zplane", "z0": 364.7250+LZC}) # Grid 7 Top
surfaces.append({"id": 70, "type": "zplane", "z0": 359.0100+LZC}) # Grid 7 Bottom
surfaces.append({"id": 61, "type": "zplane", "z0": 312.5280+LZC}) # Grid 6 Top
surfaces.append({"id": 60, "type": "zplane", "z0": 306.8130+LZC}) # Grid 6 Bottom
surfaces.append({"id": 51, "type": "zplane", "z0": 260.3310+LZC}) # Grid 5 Top
surfaces.append({"id": 50, "type": "zplane", "z0": 254.6160+LZC}) # Grid 5 Bottom
surfaces.append({"id": 41, "type": "zplane", "z0": 208.1340+LZC}) # Grid 4 Top
surfaces.append({"id": 40, "type": "zplane", "z0": 202.4190+LZC}) # Grid 4 Bottom
surfaces.append({"id": 31, "type": "zplane", "z0": 155.9370+LZC}) # Grid 3 Top
surfaces.append({"id": 30, "type": "zplane", "z0": 150.2220+LZC}) # Grid 3 Bottom
surfaces.append({"id": 21, "type": "zplane", "z0": 103.7400+LZC}) # Grid 2 Top
surfaces.append({"id": 20, "type": "zplane", "z0": 098.0250+LZC}) # Grid 2 Bottom
surfaces.append({"id": 15, "type": "zplane", "z0": 041.8280+LZC}) # Bottom of Lower Absorber (AIC)
surfaces.append({"id": 14, "type": "zplane", "z0": 040.5580+LZC}) # Bottom of Active Absorber
surfaces.append({"id": 13, "type": "zplane", "z0": 040.5200+LZC}) # Grid 1 Top
surfaces.append({"id": 12, "type": "zplane", "z0": 039.9580+LZC}) # Control Rod Step 0
surfaces.append({"id": 11, "type": "zplane", "z0": 038.6600+LZC}) # Bot. of BPRA Rod
surfaces.append({"id": 10, "type": "zplane", "z0": 037.1621+LZC}) # Grid 1 Bottom
surfaces.append({"id":  9, "type": "zplane", "z0": 036.7480+LZC}) # Bottom of Active Fuel
surfaces.append({"id":  8, "type": "zplane", "z0": 035.0000+LZC}) # Bottom of Fuel Rod
surfaces.append({"id":  7, "type": "zplane", "z0": 020.0000+LZC}) # Bottom of Support Plate
surfaces.append({"id":  6, "type": "zplane", "z0": 000.0000+LZC}) # Lowest Extent

#-------------------------------------------------------------------------
# Pin Spacer Surfaces for Inconel Grid
surfaces.append({"id": 101, "type": "xplane", "x0": -0.61015})
surfaces.append({"id": 102, "type": "xplane", "x0":  0.61015})
surfaces.append({"id": 103, "type": "yplane", "y0": -0.61015})
surfaces.append({"id": 104, "type": "yplane", "y0":  0.61015})

# Pin Spacer Surfaces for Zircaloy Grid
surfaces.append({"id": 105, "type": "xplane", "x0": -0.61049})
surfaces.append({"id": 106, "type": "xplane", "x0":  0.61049})
surfaces.append({"id": 107, "type": "yplane", "y0": -0.61049})
surfaces.append({"id": 108, "type": "yplane", "y0":  0.61049})

#-------------------------------------------------------------------------
# Fuel Pin
surfaces.append({"id": 109, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.39218}) # Fuel Rad
surfaces.append({"id": 110, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.40005}) # He Rad
surfaces.append({"id": 111, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.45720}) # Zircaloy Rad

cells.append({"id": 101, "region": "-109", "material": F16}) # Fuel 1.6%
cells.append({"id": 102, "region": "-109", "material": F24}) # Fuel 2.4%
cells.append({"id": 103, "region": "-109", "material": F31}) # Fuel 3.1%
cells.append({"id": 104, "region": "-109", "material": F32}) # Fuel 3.2%
cells.append({"id": 105, "region": "-109", "material": F34}) # Fuel 3.4%
cells.append({"id": 106, "region": "+109 & -110", "material": He}) # He Gap
cells.append({"id": 107, "region": "+110 & -111", "material": Zirc4}) # Cladding
cells.append({"id": 108, "region": "+111", "material": BH2O}) # Water, no grid
cells.append({"id": 109, "region": "+111 & ~(-101 U +102 U -103 U +104)", "material": BH2O}) # Water, with Inconcel Grid
cells.append({"id": 110, "region": "-101 U +102 U -103 U +104", "material": Inconel}) # Inconel Grid
cells.append({"id": 111, "region": "+111 & ~(-105 U +106 U -107 U +108)", "material": BH2O}) # Water, with Zircaloy Grid
cells.append({"id": 112, "region": "-105 U +106 U -107 U +108", "material": Zirc4}) # Zircaloy Grid
cells.append({"id": 113, "region": "", "material": BH2O}) # Water
cells.append({"id": 114, "region": "-111", "material": Zirc4}) # Cladding pin bottom
cells.append({"id": 115, "region": "-111", "material": NS_SS}) # Nozel
cells.append({"id": 116, "region": "+111", "material": NS_BH2O}) # Nozel water

# Upper Fuel Pin Plenum with "spring"
surfaces.append({"id": 117, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.06459}) # Inconel Rad
surfaces.append({"id": 118, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.40005}) # He Rad
surfaces.append({"id": 119, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.45720}) # Zircaloy Rad
cells.append({"id": 117, "region": "-117", "material": Inconel}) # Spring
cells.append({"id": 118, "region": "+117 & -118", "material": He})
cells.append({"id": 119, "region": "+118 & -119", "material": Zirc4})

universes.append({"id": 1601, "cells": [101, 106, 107, 108]}) # 1.6%, no grid
universes.append({"id": 2401, "cells": [102, 106, 107, 108]}) # 2.4%, no grid
universes.append({"id": 3101, "cells": [103, 106, 107, 108]}) # 3.1%, no grid
universes.append({"id": 3201, "cells": [104, 106, 107, 108]}) # 3.2%, no grid
universes.append({"id": 3401, "cells": [105, 106, 107, 108]}) # 3.4%, no grid

universes.append({"id": 1602, "cells": [101, 106, 107, 109, 110]}) # 1.6%, inconel grid
universes.append({"id": 2402, "cells": [102, 106, 107, 109, 110]}) # 2.4%, inconel grid
universes.append({"id": 3102, "cells": [103, 106, 107, 109, 110]}) # 3.1%, inconel grid
universes.append({"id": 3202, "cells": [104, 106, 107, 109, 110]}) # 3.2%, inconel grid
universes.append({"id": 3402, "cells": [105, 106, 107, 109, 110]}) # 3.4%, inconel grid

universes.append({"id": 1603, "cells": [101, 106, 107, 111, 112]}) # 1.6%, zircaloy grid
universes.append({"id": 2403, "cells": [102, 106, 107, 111, 112]}) # 2.4%, zircaloy grid
universes.append({"id": 3103, "cells": [103, 106, 107, 111, 112]}) # 3.1%, zircaloy grid
universes.append({"id": 3203, "cells": [104, 106, 107, 111, 112]}) # 3.2%, zircaloy grid
universes.append({"id": 3403, "cells": [105, 106, 107, 111, 112]}) # 3.4%, zircaloy grid

universes.append({"id": 1, "cells": [113]}) # Infinite water universe
universes.append({"id": 2, "cells": [115, 116]}) # Nozzel/support universe
universes.append({"id": 3, "cells": [114, 108]}) # Zirc pin top/bottom
universes.append({"id": 4, "cells": [117, 118, 119, 108]}) # Plenum pin cell
universes.append({"id": 5, "cells": [117, 118, 119, 109, 110]}) # Plenum pin cell with inconel grid

# Axial cells for fuel pins
cells.append({"id": 150, "region": "+95", "universe": 1}) # Water on top of everything
cells.append({"id": 151, "region": "+94 & -95", "universe": 2}) # Nozzel
cells.append({"id": 152, "region": "+92 & -94", "universe": 1}) # Water between nozel and pin
cells.append({"id": 153, "region": "+91 & -92", "universe": 3}) # Zirc pin top
cells.append({"id": 154, "region": "+81 & -91", "universe": 4}) # Plenum pin above grid 8
cells.append({"id": 155, "region": "+80 & -81", "universe": 5}) # Plenum pin with grid 8
cells.append({"id": 156, "region": "+75 & -80", "universe": 4}) # Plenum pin bellow grid 8

cells.append({"id": 160, "region": "+71 & -75", "universe": 1601}) # 1.6% AG7
cells.append({"id": 161, "region": "+70 & -71", "universe": 1603}) # 1.6% G7
cells.append({"id": 162, "region": "+61 & -70", "universe": 1601}) # 1.6% AG6
cells.append({"id": 163, "region": "+60 & -61", "universe": 1603}) # 1.6% G6
cells.append({"id": 164, "region": "+51 & -60", "universe": 1601}) # 1.6% AG5
cells.append({"id": 165, "region": "+50 & -51", "universe": 1603}) # 1.6% G5
cells.append({"id": 166, "region": "+41 & -50", "universe": 1601}) # 1.6% AG4
cells.append({"id": 167, "region": "+40 & -41", "universe": 1603}) # 1.6% G4
cells.append({"id": 168, "region": "+31 & -40", "universe": 1601}) # 1.6% AG3
cells.append({"id": 169, "region": "+30 & -31", "universe": 1603}) # 1.6% G3
cells.append({"id": 170, "region": "+21 & -30", "universe": 1601}) # 1.6% AG2
cells.append({"id": 171, "region": "+20 & -21", "universe": 1603}) # 1.6% G2
cells.append({"id": 172, "region": "+13 & -20", "universe": 1601}) # 1.6% AG1
cells.append({"id": 173, "region": "+10 & -13", "universe": 1602}) # 1.6% G1
cells.append({"id": 174, "region": "+9  & -10", "universe": 1601}) # 1.6% BG1

cells.append({"id": 180, "region": "+71 & -75", "universe": 2401}) # 2.4% AG7
cells.append({"id": 181, "region": "+70 & -71", "universe": 2403}) # 2.4% G7
cells.append({"id": 182, "region": "+61 & -70", "universe": 2401}) # 2.4% AG6
cells.append({"id": 183, "region": "+60 & -61", "universe": 2403}) # 2.4% G6
cells.append({"id": 184, "region": "+51 & -60", "universe": 2401}) # 2.4% AG5
cells.append({"id": 185, "region": "+50 & -51", "universe": 2403}) # 2.4% G5
cells.append({"id": 186, "region": "+41 & -50", "universe": 2401}) # 2.4% AG4
cells.append({"id": 187, "region": "+40 & -41", "universe": 2403}) # 2.4% G4
cells.append({"id": 188, "region": "+31 & -40", "universe": 2401}) # 2.4% AG3
cells.append({"id": 189, "region": "+30 & -31", "universe": 2403}) # 2.4% G3
cells.append({"id": 190, "region": "+21 & -30", "universe": 2401}) # 2.4% AG2
cells.append({"id": 191, "region": "+20 & -21", "universe": 2403}) # 2.4% G2
cells.append({"id": 192, "region": "+13 & -20", "universe": 2401}) # 2.4% AG1
cells.append({"id": 193, "region": "+10 & -13", "universe": 2402}) # 2.4% G1
cells.append({"id": 194, "region": "+9  & -10", "universe": 2401}) # 2.4% BG1

cells.append({"id": 200, "region": "+71 & -75", "universe": 3101}) # 3.1% AG7
cells.append({"id": 201, "region": "+70 & -71", "universe": 3103}) # 3.1% G7
cells.append({"id": 202, "region": "+61 & -70", "universe": 3101}) # 3.1% AG6
cells.append({"id": 203, "region": "+60 & -61", "universe": 3103}) # 3.1% G6
cells.append({"id": 204, "region": "+51 & -60", "universe": 3101}) # 3.1% AG5
cells.append({"id": 205, "region": "+50 & -51", "universe": 3103}) # 3.1% G5
cells.append({"id": 206, "region": "+41 & -50", "universe": 3101}) # 3.1% AG4
cells.append({"id": 207, "region": "+40 & -41", "universe": 3103}) # 3.1% G4
cells.append({"id": 208, "region": "+31 & -40", "universe": 3101}) # 3.1% AG3
cells.append({"id": 209, "region": "+30 & -31", "universe": 3103}) # 3.1% G3
cells.append({"id": 210, "region": "+21 & -30", "universe": 3101}) # 3.1% AG2
cells.append({"id": 211, "region": "+20 & -21", "universe": 3103}) # 3.1% G2
cells.append({"id": 212, "region": "+13 & -20", "universe": 3101}) # 3.1% AG1
cells.append({"id": 213, "region": "+10 & -13", "universe": 3102}) # 3.1% G1
cells.append({"id": 214, "region": "+9  & -10", "universe": 3101}) # 3.1% BG1

cells.append({"id": 220, "region": "+71 & -75", "universe": 3201}) # 3.2% AG7
cells.append({"id": 221, "region": "+70 & -71", "universe": 3203}) # 3.2% G7
cells.append({"id": 222, "region": "+61 & -70", "universe": 3201}) # 3.2% AG6
cells.append({"id": 223, "region": "+60 & -61", "universe": 3203}) # 3.2% G6
cells.append({"id": 224, "region": "+51 & -60", "universe": 3201}) # 3.2% AG5
cells.append({"id": 225, "region": "+50 & -51", "universe": 3203}) # 3.2% G5
cells.append({"id": 226, "region": "+41 & -50", "universe": 3201}) # 3.2% AG4
cells.append({"id": 227, "region": "+40 & -41", "universe": 3203}) # 3.2% G4
cells.append({"id": 228, "region": "+31 & -40", "universe": 3201}) # 3.2% AG3
cells.append({"id": 229, "region": "+30 & -31", "universe": 3203}) # 3.2% G3
cells.append({"id": 230, "region": "+21 & -30", "universe": 3201}) # 3.2% AG2
cells.append({"id": 231, "region": "+20 & -21", "universe": 3203}) # 3.2% G2
cells.append({"id": 232, "region": "+13 & -20", "universe": 3201}) # 3.2% AG1
cells.append({"id": 233, "region": "+10 & -13", "universe": 3202}) # 3.2% G1
cells.append({"id": 234, "region": "+9  & -10", "universe": 3201}) # 3.2% BG1

cells.append({"id": 240, "region": "+71 & -75", "universe": 3201}) # 3.4% AG7
cells.append({"id": 241, "region": "+70 & -71", "universe": 3203}) # 3.4% G7
cells.append({"id": 242, "region": "+61 & -70", "universe": 3201}) # 3.4% AG6
cells.append({"id": 243, "region": "+60 & -61", "universe": 3203}) # 3.4% G6
cells.append({"id": 244, "region": "+51 & -60", "universe": 3201}) # 3.4% AG5
cells.append({"id": 245, "region": "+50 & -51", "universe": 3203}) # 3.4% G5
cells.append({"id": 246, "region": "+41 & -50", "universe": 3201}) # 3.4% AG4
cells.append({"id": 247, "region": "+40 & -41", "universe": 3203}) # 3.4% G4
cells.append({"id": 248, "region": "+31 & -40", "universe": 3201}) # 3.4% AG3
cells.append({"id": 249, "region": "+30 & -31", "universe": 3203}) # 3.4% G3
cells.append({"id": 250, "region": "+21 & -30", "universe": 3201}) # 3.4% AG2
cells.append({"id": 251, "region": "+20 & -21", "universe": 3203}) # 3.4% G2
cells.append({"id": 252, "region": "+13 & -20", "universe": 3201}) # 3.4% AG1
cells.append({"id": 253, "region": "+10 & -13", "universe": 3202}) # 3.4% G1
cells.append({"id": 254, "region": "+9 & -10", "universe": 3201}) # 3.4% BG1

cells.append({"id": 157, "region": "+8 & -9", "universe": 3}) # Zirc pin bottom
cells.append({"id": 158, "region": "+7 & -8", "universe": 2}) # Support plate
cells.append({"id": 159, "region": "-7", "universe": 1}) # Water bellow support plate

universes.append({"id": 1600, "cells": [150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                        160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
                                        170, 171, 172, 173, 174]}) # 1.6% Fuel Pin Universe

universes.append({"id": 2400, "cells": [150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                        180, 181, 182, 183, 184, 185, 186, 187, 188, 189,
                                        190, 191, 192, 193, 194]}) # 2.4% Fuel Pin Universe

universes.append({"id": 3100, "cells": [150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                        200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
                                        210, 211, 212, 213, 214]}) # 3.1% Fuel Pin Universe

universes.append({"id": 3200, "cells": [150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                        220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
                                        230, 231, 232, 233, 234]}) # 3.2% Fuel Pin Universe

universes.append({"id": 3400, "cells": [150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                        240, 241, 242, 243, 244, 245, 246, 247, 248, 249,
                                        250, 251, 252, 253, 254]}) # 3.4% Fuel Pin Universe

# Used Cells 100-259
# Used Universes 1-5, 1600-1603, 2400-2404, 3100-3103, 3200-3203, 3400-3403
# Used Surfaces 1-119

#---------------------------------------------------------------------------
# Empty Guide Tube at Dashpot
surfaces.append({"id": 301, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.50419}) # Water rad
surfaces.append({"id": 302, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.54610}) # Zircaloy rad

cells.append({"id": 300, "region": "-301", "material": BH2O}) # Water in GT at Dashpot
cells.append({"id": 301, "region": "-302 & +301", "material": Zirc4})
cells.append({"id": 302, "region": "+302", "material": BH2O}) # Water outside GT at Dashpot, no grid
cells.append({"id": 303, "region": "+302 & ~(-101 U +102 U -103 U +104)", "material": BH2O}) # Water outside GT at Dashpot, inconel grid
cells.append({"id": 304, "region": "-101 U +102 U -103 U +104", "material": Inconel}) # Inconel grid

universes.append({"id": 31, "cells": [300, 301, 302]}) # GT at Dashpot, no grid
universes.append({"id": 32, "cells": [300, 301, 303, 304]}) # GT at Dashpot, inconel grid

# Empty Guide Tube above Dashpot
surfaces.append({"id": 303, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.56134}) # Water Rad
surfaces.append({"id": 304, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.60198}) # Zircaloy Rad

cells.append({"id": 305, "region": "-303", "material": BH2O}) # Water in GT above dashpot
cells.append({"id": 306, "region": "-304 & + 303", "material": Zirc4}) # Zirc GT above dashpot
cells.append({"id": 307, "region": "+304", "material": BH2O}) # Water around GT, no grid
cells.append({"id": 308, "region": "+304 & ~(-101 U +102 U -103 U +104)", "material": BH2O}) # Water around GT, inconel grid
cells.append({"id": 309, "region": "+304 & ~(-105 U +106 U -107 U +108)", "material": BH2O}) # Water around GT, zircaloy grid
cells.append({"id": 310, "region": "-105 U +106 U -107 U +108", "material": Zirc4}) # zircaloy grid
cells.append({"id": 311, "region": "", "material": NS_BH2O}) # Infinite water for the nozzel of GT

universes.append({"id": 33, "cells": [305, 306, 307]}) # GT above dashpot, no grid
universes.append({"id": 34, "cells": [305, 306, 308, 304]}) # GT above dashpot, inconel grid
universes.append({"id": 35, "cells": [305, 306, 309, 310]}) # GT above dashpot, zircaloy grid
universes.append({"id": 36, "cells": [311]}) # Nozzel/support for GT

cells.append({"id": 320, "region": "+95", "universe": 1}) # Water above nozzel
cells.append({"id": 321, "region": "+94 & -95", "universe": 36}) # Nozzel
cells.append({"id": 322, "region": "+81 & -94", "universe": 33}) # GT AG8
cells.append({"id": 323, "region": "+80 & -81", "universe": 34}) # GT G8
cells.append({"id": 324, "region": "+71 & -80", "universe": 33}) # GT AG7
cells.append({"id": 325, "region": "+70 & -71", "universe": 35}) # GT G7
cells.append({"id": 326, "region": "+61 & -70", "universe": 33}) # GT AG6
cells.append({"id": 327, "region": "+60 & -61", "universe": 35}) # GT G6
cells.append({"id": 328, "region": "+51 & -60", "universe": 33}) # GT AG5
cells.append({"id": 329, "region": "+50 & -51", "universe": 35}) # GT G5
cells.append({"id": 330, "region": "+41 & -50", "universe": 33}) # GT AG4
cells.append({"id": 331, "region": "+40 & -41", "universe": 35}) # GT G4
cells.append({"id": 332, "region": "+31 & -40", "universe": 33}) # GT AG3
cells.append({"id": 333, "region": "+30 & -31", "universe": 35}) # GT G3
cells.append({"id": 334, "region": "+21 & -30", "universe": 33}) # GT AG2
cells.append({"id": 335, "region": "+20 & -21", "universe": 35}) # GT G2
cells.append({"id": 336, "region": "+13 & -20", "universe": 33}) # GT AG1
cells.append({"id": 337, "region": "+12 & -13", "universe": 34}) # GT G1, ADP
cells.append({"id": 338, "region": "+10 & -12", "universe": 32}) # GT G1, BDP
cells.append({"id": 339, "region": "+8  & -10", "universe": 31}) # GT BG1, BDP
cells.append({"id": 340, "region": "+7  & -8",  "universe": 36}) # Support plate water
cells.append({"id": 341, "region": "-7",        "universe": 1}) # Water bellow support

universes.append({"id": 300, "cells": [320, 321, 322, 323, 324, 325, 326, 327, 328, 329,
                                       330, 331, 332, 333, 334, 335, 336, 337, 338, 339,
                                       340, 341]}) # GT Universe

#---------------------------------------------------------------------------
# Instrument Tube Pin
surfaces.append({"id": 401, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.43688}) # Air Rad
surfaces.append({"id": 402, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.48387}) # Zircaloy 1 Rad
surfaces.append({"id": 403, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.56134}) # Water Rad
surfaces.append({"id": 404, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.60198}) # Zircaloy 2 Rad

cells.append({"id": 40, "region": "-401", "material": Air}) # Air in IT
cells.append({"id": 41, "region": "+401 & -402", "material": Zirc4})
cells.append({"id": 42, "region": "+402 & -403", "material": BH2O})
cells.append({"id": 43, "region": "+403 & -404", "material": Zirc4})
cells.append({"id": 44, "region": "+404", "material": BH2O}) # Water, no grid
cells.append({"id": 45, "region": "+404 & ~(-101 U +102 U -103 U +104)", "material": BH2O}) # Water, inconel grid
cells.append({"id": 46, "region": "-101 U +102 U -103 U +104", "material": Inconel})
cells.append({"id": 47, "region": "+404 & ~(-105 U +106 U -107 U +108)", "material": BH2O}) # Water, zircaloy grid
cells.append({"id": 48, "region": "-105 U +106 U -107 U +108", "material": Zirc4})
cells.append({"id": 49, "region": "+404", "material": NS_BH2O}) # Nozzel water around IT

universes.append({"id": 40, "cells": [40, 41, 42, 43, 44]}) # IT, no grid
universes.append({"id": 41, "cells": [40, 41, 42, 43, 45, 46]}) # IT, inconel grid
universes.append({"id": 42, "cells": [40, 41, 42, 43, 47, 48]}) # IT, zircaloy grid
universes.append({"id": 43, "cells": [40, 41, 42, 43, 49]}) # IT, no grid, nozzel

# Bare Instrument Thimble Pin
surfaces.append({"id": 405, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.43688}) # Air Rad
surfaces.append({"id": 406, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.48387}) # Zircaloy Rad

cells.append({"id": 50, "region": "+402", "material": BH2O}) # Water around Bare IT

universes.append({"id": 44, "cells": [40, 41, 50]}) # Bare IT with water

cells.append({"id": 420, "region": "+95", "universe": 1}) # Water above nozzel
cells.append({"id": 421, "region": "+94 & -95", "universe": 36}) # Nozzel (from GT)
cells.append({"id": 422, "region": "+81 & -94", "universe": 40}) # IT AG8
cells.append({"id": 423, "region": "+80 & -81", "universe": 41}) # IT G8
cells.append({"id": 424, "region": "+71 & -80", "universe": 40}) # IT AG7
cells.append({"id": 425, "region": "+70 & -71", "universe": 42}) # IT G7
cells.append({"id": 426, "region": "+61 & -70", "universe": 40}) # IT AG6
cells.append({"id": 427, "region": "+60 & -61", "universe": 42}) # IT G6
cells.append({"id": 428, "region": "+51 & -60", "universe": 40}) # IT AG5
cells.append({"id": 429, "region": "+50 & -51", "universe": 42}) # IT G5
cells.append({"id": 430, "region": "+41 & -50", "universe": 40}) # IT AG4
cells.append({"id": 431, "region": "+40 & -41", "universe": 42}) # IT G4
cells.append({"id": 432, "region": "+31 & -40", "universe": 40}) # IT AG3
cells.append({"id": 433, "region": "+30 & -31", "universe": 42}) # IT G3
cells.append({"id": 434, "region": "+21 & -30", "universe": 40}) # IT AG2
cells.append({"id": 435, "region": "+20 & -21", "universe": 42}) # IT G2
cells.append({"id": 436, "region": "+13 & -20", "universe": 40}) # IT AG1
cells.append({"id": 437, "region": "+10 & -13", "universe": 41}) # IT G1
cells.append({"id": 438, "region": "+8  & -10", "universe": 40}) # IT BG1
cells.append({"id": 439, "region": "+7  & -8",  "universe": 43}) # Support plate bare IT
cells.append({"id": 440, "region": "-7",        "universe": 44}) # Water bare IT

universes.append({"id": 400, "cells": [420, 421, 422, 423, 424, 425, 426, 427, 428, 429,
                                       430, 431, 432, 433, 434, 435, 436, 437, 438, 439,
                                       440]}) # Instrument Tube Universe

#---------------------------------------------------------------------------
# BP above Dashpot
surfaces.append({"id": 501, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.21400}) # Air Rad
surfaces.append({"id": 502, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.23051}) # SS304 Rad
surfaces.append({"id": 503, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.24130}) # He Rad
surfaces.append({"id": 504, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.42672}) # Borosilicate Glass Rad
surfaces.append({"id": 505, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.43688}) # He Rad
surfaces.append({"id": 506, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.48387}) # SS304 Rad
surfaces.append({"id": 507, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.56134}) # Water Rad
surfaces.append({"id": 508, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.60198}) # Zircaloy Rad

cells.append({"id": 51, "region": "-501", "material": Air}) # Core of BP above DP
cells.append({"id": 52, "region": "+501 & -502", "material": SS304})
cells.append({"id": 53, "region": "+502 & -503", "material": He})
cells.append({"id": 54, "region": "+503 & -504", "material": BSiGlass})
cells.append({"id": 55, "region": "+504 & -505", "material": He})
cells.append({"id": 56, "region": "+505 & -506", "material": SS304})
cells.append({"id": 57, "region": "+506 & -507", "material": BH2O})
cells.append({"id": 58, "region": "+507 & -508", "material": Zirc4})
cells.append({"id": 59, "region": "+508", "material": BH2O}) # Water no grid
cells.append({"id": 60, "region": "+508 & ~(-105 U +106 U -107 U +108)", "material": BH2O}) # Wate with zirc grid
cells.append({"id": 61, "region": "-105 U +106 U -107 U +108", "material": Zirc4}) # Zirc Grid

universes.append({"id": 51, "cells": [51, 52, 53, 54, 55, 56, 57, 58, 59]}) # BP no grid
universes.append({"id": 52, "cells": [51, 52, 53, 54, 55, 56, 57, 58, 60, 61]}) # BP zirc grid

cells.append({"id": 62, "region": "+502 & -505", "material": He}) # Plenum He
cells.append({"id": 63, "region": "+508 & ~(-101 U +102 U -103 U +104)", "material": BH2O}) # BP plenum water with inconel grid
cells.append({"id": 64, "region": "-101 U +102 U -103 U +104", "material": Inconel}) # Inconel Grid

universes.append({"id": 53, "cells": [51, 52, 62, 56, 57, 58, 59]}) # BP plenum no grid
universes.append({"id": 54, "cells": [51, 52, 62, 56, 57, 58, 63, 64]}) # BP plenum inconel grid

cells.append({"id": 65, "region": "-506", "material": SS304}) # SS pin
cells.append({"id": 66, "region": "+506", "material": NS_BH2O}) # SS pin surounded by nozzel water

universes.append({"id": 55, "cells": [65, 57, 58, 59]}) # SS pin in GT
universes.append({"id": 56, "cells": [65, 66]}) # SS pin nozzel
universes.append({"id": 57, "cells": [65, 57, 58, 63, 64]}) # SS pin in GT with inconel grid

cells.append({"id": 67, "region": "+506 & -301", "material": BH2O}) # Water between SS pin and Dashpot GT

universes.append({"id": 58, "cells": [65, 67, 301, 303, 304]}) # SS pin in DP GT with inconel grid

cells.append({"id": 70, "region": "+95", "universe": 1}) # Water
cells.append({"id": 71, "region": "+94 & -95", "universe": 56}) # SS pin in nozzel
cells.append({"id": 72, "region": "+93 & -94", "universe": 55}) # SS pin in GT
cells.append({"id": 73, "region": "+81 & -93", "universe": 53}) # BP plenum no grid
cells.append({"id": 74, "region": "+80 & -81", "universe": 54}) # BP plenum with inconel grid
cells.append({"id": 75, "region": "+74 & -80", "universe": 53}) # BP plenum no grid
cells.append({"id": 76, "region": "+71 & -74", "universe": 51}) # BP AG7
cells.append({"id": 77, "region": "+70 & -71", "universe": 52}) # BP G7
cells.append({"id": 78, "region": "+61 & -70", "universe": 51}) # BP AG6
cells.append({"id": 79, "region": "+60 & -61", "universe": 52}) # BP G6
cells.append({"id": 80, "region": "+51 & -60", "universe": 51}) # BP AG5
cells.append({"id": 81, "region": "+50 & -51", "universe": 52}) # BP G5
cells.append({"id": 82, "region": "+41 & -50", "universe": 51}) # BP AG4
cells.append({"id": 83, "region": "+40 & -41", "universe": 52}) # BP G4
cells.append({"id": 84, "region": "+31 & -40", "universe": 51}) # BP AG3
cells.append({"id": 85, "region": "+30 & -31", "universe": 52}) # BP G3
cells.append({"id": 86, "region": "+21 & -30", "universe": 51}) # BP AG2
cells.append({"id": 87, "region": "+20 & -21", "universe": 52}) # BP G2
cells.append({"id": 88, "region": "+14 & -20", "universe": 51}) # BP AG1
cells.append({"id": 89, "region": "+13 & -14", "universe": 55}) # SS pin in GT
cells.append({"id": 90, "region": "+12 & -13", "universe": 57}) # SS pin in GT with inconel grid
cells.append({"id": 91, "region": "+11 & -12", "universe": 58}) # SS pin in DP GT with inconel grid
cells.append({"id": 92, "region": "+10 & -11", "universe": 32}) # DP GT with inconel grid
cells.append({"id": 93, "region": "+8  & -10", "universe": 31}) # DP GT
cells.append({"id": 94, "region": "+7 & -8", "universe": 36}) # Support water bellow DP GT
cells.append({"id": 95, "region": "-7", "universe": 1}) # Water bellow

universes.append({"id": 500, "cells": [70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                       80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
                                       90, 91, 92, 93, 94, 95]}) # BP Universe

#-----------------------------------------------------------------------------
# Control Rod Thimble ONLY

# Control Rod Pin Upper
surfaces.append({"id": 901, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.37338}) # B4C Rad
surfaces.append({"id": 902, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.38608}) # He Rad
surfaces.append({"id": 903, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.48387}) # SS304 Rad

# Control Rod Pin Lower
surfaces.append({"id": 904, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.38227}) # Ag-In-CD Rad
surfaces.append({"id": 905, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.38608}) # He Rad

# Control Rod Pin Spacer
surfaces.append({"id": 906, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.37845}) # SS304 Rad
surfaces.append({"id": 907, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.37845}) # He Rad

# Control Rod Pin Plenum
surfaces.append({"id": 908, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.06459}) # Inconcel Rad
surfaces.append({"id": 909, "type": "zcylinder", "x0": 0., "y0": 0., "r": 0.38608}) # He Rad


cells.append({"id": 900, "region": "-903", "material": SS304}) # SS pin
universes.append({"id": 900, "cells": [900]}) # SS pin uni

cells.append({"id": 901, "region": "-901", "material": B4C})
cells.append({"id": 902, "region": "+901 & -902", "material": He})
cells.append({"id": 903, "region": "+902 & -903", "material": SS304})
universes.append({"id": 901, "cells": [901, 902, 903]}) # CR Pin Upper Absorber

cells.append({"id": 904, "region": "-904", "material": AgInCd})
cells.append({"id": 905, "region": "+904 & -905", "material": He})
cells.append({"id": 906, "region": "+905 & -903", "material": SS304})
universes.append({"id": 902, "cells": [904, 905, 906]}) # CR Pin Lower Absorber

cells.append({"id": 907, "region": "-906", "material": SS304})
cells.append({"id": 908, "region": "+906 & -907", "material": He})
cells.append({"id": 909, "region": "+907 & -903", "material": SS304})
universes.append({"id": 903, "cells": [907, 908, 909]}) # CR Pin Spacer

cells.append({"id": 910, "region": "-908", "material": Inconel})
cells.append({"id": 911, "region": "+908 & -909", "material": He})
cells.append({"id": 912, "region": "+909 & -903", "material": SS304})
universes.append({"id": 904, "cells": [910, 911, 912]}) # CR Pin Plenum

cells.append({"id": 913, "region": "-903", "material": BH2O})
universes.append({"id": 905, "cells": [913]}) # CR Water universe (cylinder only)

cells.append({"id": 914, "region": "-903", "material": NS_BH2O})
universes.append({"id": 906, "cells": [914]}) # CR Nozzel water universe (cylinder only)

def make_cr_thimble_uni(ZC, id):
  surfaces.append({"id": id+1, "type": "zplane", "z0": 415.558+LZC+ZC}) # Top of CR Plenum
  surfaces.append({"id": id+2, "type": "zplane", "z0": 403.778+LZC+ZC}) # Bottom of CR Plenum
  surfaces.append({"id": id+3, "type": "zplane", "z0": 402.508+LZC+ZC}) # Bottom of Spacer
  surfaces.append({"id": id+4, "type": "zplane", "z0": 143.428+LZC+ZC}) # Bottom of Upper Absorber
  surfaces.append({"id": id+5, "type": "zplane", "z0": 041.828+LZC+ZC}) # Bottom of Lower Absorber
  surfaces.append({"id": id+6, "type": "zplane", "z0": 039.958+LZC+ZC}) # Bottom of CR

  cells.append({"id": id+1, "region": "+{}".format(id+1), "universe": 900}) # SS pin top
  cells.append({"id": id+2, "region": "+{} & -{}".format(id+2, id+1), "universe": 904}) # CR upper plenum
  cells.append({"id": id+3, "region": "+{} & -{}".format(id+3, id+2), "universe": 903}) # CR Spacer
  cells.append({"id": id+4, "region": "+{} & -{}".format(id+4, id+3), "universe": 901}) # CR Upper absorber
  cells.append({"id": id+5, "region": "+{} & -{}".format(id+5, id+4), "universe": 902}) # CR Upper absorber
  cells.append({"id": id+6, "region": "+{} & -{}".format(id+6, id+5), "universe": 900}) # SS pin bottom
  cells.append({"id": id+7, "region": "+8 & -{}".format(id+6), "universe": 905}) # Water bellow CR
  cells.append({"id": id+8, "region": "+7 & -8", "universe": 906}) # Water bellow CR
  cells.append({"id": id+9, "region": "-7", "universe": 905}) # Water to bottom

  universes.append({"id": id, "cells": list(range(id+1, id+10))})
  return id

A_thmbl_uni = make_cr_thimble_uni(AZC, 99900)
B_thmbl_uni = make_cr_thimble_uni(BZC, 99800)
C_thmbl_uni = make_cr_thimble_uni(CZC, 99700)
D_thmbl_uni = make_cr_thimble_uni(DZC, 99600)
a_thmbl_uni = make_cr_thimble_uni(aZC, 99500)
b_thmbl_uni = make_cr_thimble_uni(bZC, 99400)
c_thmbl_uni = make_cr_thimble_uni(cZC, 99300)
d_thmbl_uni = make_cr_thimble_uni(dZC, 99200)
e_thmbl_uni = make_cr_thimble_uni(eZC, 99100)


cells.append({"id": 915, "region": "+903", "universe": 300}) # This is the GT uni, with middle cut out
cells.append({"id": 921, "region": "-903", "universe": A_thmbl_uni}) # CR A core
cells.append({"id": 922, "region": "-903", "universe": B_thmbl_uni}) # CR B core
cells.append({"id": 923, "region": "-903", "universe": C_thmbl_uni}) # CR C core
cells.append({"id": 924, "region": "-903", "universe": D_thmbl_uni}) # CR D core
cells.append({"id": 925, "region": "-903", "universe": a_thmbl_uni}) # CR Shutdown A core
cells.append({"id": 926, "region": "-903", "universe": b_thmbl_uni}) # CR Shutdown B core
cells.append({"id": 927, "region": "-903", "universe": c_thmbl_uni}) # CR Shutdown C core
cells.append({"id": 928, "region": "-903", "universe": d_thmbl_uni}) # CR Shutdown D core
cells.append({"id": 929, "region": "-903", "universe": e_thmbl_uni}) # CR Shutdown E core

CRA = 911
CRB = 912
CRC = 913
CRD = 914
CRa = 915
CRb = 916
CRc = 917
CRd = 918
CRe = 919
universes.append({"id": CRA, "cells": [915, 921]}) # CR A
universes.append({"id": CRB, "cells": [915, 922]}) # CR B
universes.append({"id": CRC, "cells": [915, 923]}) # CR C
universes.append({"id": CRD, "cells": [915, 924]}) # CR D
universes.append({"id": CRa, "cells": [915, 925]}) # CR Shutdown A
universes.append({"id": CRb, "cells": [915, 926]}) # CR Shutdown B
universes.append({"id": CRc, "cells": [915, 927]}) # CR Shutdown C
universes.append({"id": CRd, "cells": [915, 928]}) # CR Shutdown D
universes.append({"id": CRe, "cells": [915, 929]}) # CR Shutdown E


#------------------------------------------------------------------------------
# Assembly Grid Sleves

# Inner surface of sleve
surfaces.append({"id": 1201, "type": "xplane", "x0": -10.70864})
surfaces.append({"id": 1202, "type": "xplane", "x0":  10.70864})
surfaces.append({"id": 1203, "type": "yplane", "y0": -10.70864})
surfaces.append({"id": 1204, "type": "yplane", "y0":  10.70864})

# Outer surface of sleve
surfaces.append({"id": 1205, "type": "xplane", "x0": -10.74798})
surfaces.append({"id": 1206, "type": "xplane", "x0":  10.74798})
surfaces.append({"id": 1207, "type": "yplane", "y0": -10.74798})
surfaces.append({"id": 1208, "type": "yplane", "y0":  10.74798})

cells.append({"id": 1200, "region": "+1201 & -1202 & +1203 & -1204", "material": BH2O}) # Water inside sleve (filled with assembly)
cells.append({"id": 1201, "region": "+1205 & -1206 & +1207 & -1208 & ~(+1201 & -1202 & +1203 & -1204)", "material": SS304}) # SS sleve
cells.append({"id": 1202, "region": "+1205 & -1206 & +1207 & -1208 & ~(+1201 & -1202 & +1203 & -1204)", "material": Zirc4}) # Zircaloy sleve
cells.append({"id": 1203, "region": "-1205 U +1206 U -1207 U +1208", "material": BH2O}) # Water outside sleve

universes.append({"id": 1201, "cells": [1200, 1201, 1203]}) # SS sleve uni
universes.append({"id": 1202, "cells": [1200, 1202, 1203]}) # Zircaloy sleve uni

cells.append({"id": 1210, "region": "+95", "material": BH2O}) # Water on top of nozzel
cells.append({"id": 1211, "region": "+94 & -95", "material": NS_BH2O}) # Nozzel water
cells.append({"id": 1212, "region": "+81 & -94", "material": BH2O}) # Water above G8
cells.append({"id": 1213, "region": "+80 & -81", "universe": 1201}) # SS sleve G8
cells.append({"id": 1214, "region": "+71 & -80", "material": BH2O}) # Water
cells.append({"id": 1215, "region": "+70 & -71", "universe": 1202}) # Zirc sleve G7
cells.append({"id": 1216, "region": "+61 & -70", "material": BH2O}) # Water
cells.append({"id": 1217, "region": "+60 & -61", "universe": 1202}) # Zirc sleve G6
cells.append({"id": 1218, "region": "+51 & -60", "material": BH2O}) # Water
cells.append({"id": 1219, "region": "+50 & -51", "universe": 1202}) # Zirc sleve G5
cells.append({"id": 1220, "region": "+41 & -50", "material": BH2O}) # Water
cells.append({"id": 1221, "region": "+40 & -41", "universe": 1202}) # Zirc sleve G4
cells.append({"id": 1222, "region": "+31 & -40", "material": BH2O}) # Water
cells.append({"id": 1223, "region": "+30 & -31", "universe": 1202}) # Zirc sleve G3
cells.append({"id": 1224, "region": "+21 & -30", "material": BH2O}) # Water
cells.append({"id": 1225, "region": "+20 & -21", "universe": 1202}) # Zirc sleve G2
cells.append({"id": 1226, "region": "+13 & -20", "material": BH2O}) # Water
cells.append({"id": 1227, "region": "+10 & -13", "universe": 1201}) # SS sleve G1
cells.append({"id": 1228, "region": "+8 & -10", "material": BH2O}) # Water
cells.append({"id": 1229, "region": "+7 & -8", "material": NS_BH2O}) # Nozzel water
cells.append({"id": 1230, "region": "-7", "material": BH2O}) # Water down to bottom

universes.append({"id": 1200, "cells": [1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219,
                                        1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229,
                                        1230]}) # Universe for assembly sleves

def L(i,F,G,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "Normal Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
                  F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
                  F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

def L6(i,F,G,B,Dir,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  if Dir not in ['N', 'S', 'E', 'W']:
    raise RuntimeError("Invalid direction")
  
  if Dir == 'N':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, B, F, F, G, F, F, B, F, F, F, F, F,
            F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]

  elif Dir == 'S':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
            F, F, F, F, F, B, F, F, G, F, F, B, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]

  elif Dir == 'E':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, B, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, B, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, B, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, B, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]

  elif Dir == 'W':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, B, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, B, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, B, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, B, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]


  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "Normal Assembly",
    "universes": unis})
  return 10000+i

def L4(i,F,G,B,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "4BA Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, G, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

def L8(i,F,G,B,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "8BA Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, B, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, B, F, F, I, F, F, B, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, B, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

def L12_C1(i,F,G,B,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "12BA Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, B, F, F, G, F, F, B, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, G, F, F, I, F, F, G, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, B, F, F, G, F, F, B, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

def L12_C2(i,F,G,B,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "12BA Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, B, F, F, G, F, F, B, F, F, F, F, F,
                  F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, B, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, G, F, F, B, F, F, I, F, F, B, F, F, G, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, B, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
                  F, F, F, F, F, B, F, F, G, F, F, B, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

def L15(i,F,G,B,Dir,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  if Dir not in ['NW', 'NE', 'SE', 'SW']:
    raise RuntimeError("Invalid direction")
  
  if Dir == 'NW':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
            F, F, F, B, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, B, F, F, B, F, F, B, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, B, F, F, I, F, F, B, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, B, F, F, B, F, F, B, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]

  elif Dir == 'NE':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, B, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, B, F, F, B, F, F, B, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, B, F, F, I, F, F, B, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, B, F, F, B, F, F, B, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]

  elif Dir == 'SE':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, B, F, F, B, F, F, B, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, B, F, F, I, F, F, B, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, G, F, F, B, F, F, B, F, F, B, F, F, B, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, B, F, F, F,
            F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]

  elif Dir == 'SW':
    unis = [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, G, F, F, G, F, F, G, F, F, F, F, F,
            F, F, F, G, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, B, F, F, B, F, F, B, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, B, F, F, I, F, F, B, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, B, F, F, B, F, F, B, F, F, B, F, F, G, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, B, F, F, F, F, F, F, F, F, F, G, F, F, F,
            F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
            F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F]


  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "Normal Assembly",
    "universes": unis})
  return 10000+i

def L16(i,F,G,B,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "16BA Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, I, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, G, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

def L20(i,F,G,B,I=None):
  # If I isn't given, this assembly doesn't have an instrument tube,
  # so we set it be to a guide tube.
  if I is None:
    I = G

  universes.append({"id": 10000+i, "type": "rectlinear", "shape": [17,17,1],
    "pitch": [1.25984, 1.25984, 460.], "outer": 1200,
    "origin": [0., 0., 0.], "name": "20BA Assembly",
    "universes": [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, B, F, F, G, F, F, B, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, G, F, F, I, F, F, G, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, B, F, F, B, F, F, G, F, F, B, F, F, B, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, B, F, F, F, F, F, F, F, F, F, B, F, F, F,
                  F, F, F, F, F, B, F, F, B, F, F, B, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
                 ]})
  return 10000+i

#===============================================================================
# Geometry

#-------------------------------------------------------------------------------
# Outside Universe
# Surfaces for the pressure vessel
surfaces.append({"id": 17001, "type": "zcylinder", "x0": 0., "y0": 0., "r": 219.15}) # PV Liner inner radius
surfaces.append({"id": 17002, "type": "zcylinder", "x0": 0., "y0": 0., "r": 219.71}) # PV Liner outer radius
surfaces.append({"id": 17003, "type": "zcylinder", "x0": 0., "y0": 0., "r": 241.3, "boundary": "vacuum"}) # PV outer radius
surfaces.append({"id": 17004, "type": "zplane", "z0": 0., "boundary": "vacuum"}) # Bottom of reactor
surfaces.append({"id": 17005, "type": "zplane", "z0": 460., "boundary": "vacuum"}) # Top of reactor
# Neutron reflectors
surfaces.append({"id": 17006, "type": "zcylinder", "x0": 0., "y0": 0., "r": 194.84}) # Inner radius of neutron sheild
surfaces.append({"id": 17007, "type": "zcylinder", "x0": 0., "y0": 0., "r": 201.63}) # Outer radius of neutron sheild
surfaces.append({"id": 17008, "type": "plane", "A": float(np.tan(29.*np.pi/180.)), "B": -1, "C": 0., "D": 0.})
surfaces.append({"id": 17009, "type": "plane", "A": float(np.tan(61.*np.pi/180.)), "B": -1, "C": 0., "D": 0.})
surfaces.append({"id": 17010, "type": "plane", "A": -float(np.tan(29.*np.pi/180.)), "B": -1, "C": 0., "D": 0.})
surfaces.append({"id": 17011, "type": "plane", "A": -float(np.tan(61.*np.pi/180.)), "B": -1, "C": 0., "D": 0.})
# Core barrel
surfaces.append({"id": 17012, "type": "zcylinder", "x0": 0., "y0": 0., "r": 187.96}) # Inner radius of Core barrel
surfaces.append({"id": 17013, "type": "zcylinder", "x0": 0., "y0": 0., "r": 193.675}) # Outer radius of Core barrel

cells.append({"id": 17001, "region": "-17002 & +17001 & +17004 & -17005", "material": SS304}) # PV Liner
cells.append({"id": 17002, "region": "-17003 & +17002 & +17004 & -17005", "material": CS}) # PV

cells.append({"id": 17003, "region": "-17012 & +17004 & -17005", "material": BH2O}) # Water inside core barel
cells.append({"id": 17004, "region": "-17013 & +17012 & +17004 & -17005", "material": SS304}) # Core barel

cells.append({"id": 17005, "region": "-17001 & +17013 & +17004 & -17005 & -17008 & +17010", "material": BH2O}) # N Water reflector
cells.append({"id": 17006, "region": "-17001 & +17013 & +17004 & -17005 & +17008 & -17010", "material": BH2O}) # S Water reflector
cells.append({"id": 17007, "region": "-17001 & +17013 & +17004 & -17005 & -17009 & -17011", "material": BH2O}) # E Water reflector
cells.append({"id": 17008, "region": "-17001 & +17013 & +17004 & -17005 & +17009 & +17011", "material": BH2O}) # W Water reflector

cells.append({"id": 17009, "region": "+17013 & -17006 & +17004 & -17005 & -17008 & +17009", "material": BH2O}) # SE inner water
cells.append({"id": 17010, "region": "+17006 & -17007 & +17004 & -17005 & -17008 & +17009", "material": SS304}) # SE sheild
cells.append({"id": 17011, "region": "+17007 & -17001 & +17004 & -17005 & -17008 & +17009", "material": BH2O}) # SE outer water

cells.append({"id": 17012, "region": "+17013 & -17006 & +17004 & -17005 & +17008 & -17009", "material": BH2O}) # NW inner water
cells.append({"id": 17013, "region": "+17006 & -17007 & +17004 & -17005 & +17008 & -17009", "material": SS304}) # NW sheild
cells.append({"id": 17014, "region": "+17007 & -17001 & +17004 & -17005 & +17008 & -17009", "material": BH2O}) # NW outer water

cells.append({"id": 17015, "region": "+17013 & -17006 & +17004 & -17005 & +17010 & -17011", "material": BH2O}) # SW inner water
cells.append({"id": 17016, "region": "+17006 & -17007 & +17004 & -17005 & +17010 & -17011", "material": SS304}) # SW sheild
cells.append({"id": 17017, "region": "+17007 & -17001 & +17004 & -17005 & +17010 & -17011", "material": BH2O}) # SW outer water

cells.append({"id": 17018, "region": "+17013 & -17006 & +17004 & -17005 & -17010 & +17011", "material": BH2O}) # NE inner water
cells.append({"id": 17019, "region": "+17006 & -17007 & +17004 & -17005 & -17010 & +17011", "material": SS304}) # NE sheild
cells.append({"id": 17020, "region": "+17007 & -17001 & +17004 & -17005 & -17010 & +17011", "material": BH2O}) # NE outer water

universes.append({"id": 17001, "cells": [17001, 17002, 17003, 17004, 17005, 17006, 17007, 17008, 17009, 17010,
  17011, 17012, 17013, 17014, 17015, 17016, 17017, 17018, 17019, 17020]})

#=============================================================================
# Lattices

EMPTY_____ = -1
F16_______ = L(1,1600,300)
F16______A = L(2,1600,CRA)
F16______B = L(3,1600,CRB)
F16______C = L(4,1600,CRC)
F16______D = L(5,1600,CRD)
F16______b = L(6,1600,CRb)
F16______c = L(7,1600,CRc)
F16______d = L(8,1600,CRd)
F16______e = L(9,1600,CRe)
F24_______ = L(10,2400,300)
F24______D = L(11,2400,CRD)
F31_______ = L(12,3100,300)
F31______a = L(13,3100,CRa)
F31_20____ = L20(14,3100,300,500)
F31_16____ = L16(15,3100,300,500)
F31_06N___ = L6(16,3100,300,500,'N')
F31_06S___ = L6(17,3100,300,500,'S')
F31_06E___ = L6(18,3100,300,500,'E')
F31_06W___ = L6(19,3100,300,500,'W')
F31_15NW__ = L15(20,3100,300,500,'NW')
F31_15NE__ = L15(21,3100,300,500,'NE')
F31_15SW__ = L15(22,3100,300,500,'SW')
F31_15SE__ = L15(23,3100,300,500,'SE')
F24_16____ = L16(24,2400,300,500)
F24_12____ = L12_C1(25,2400,300,500)

#                 R           P           N           M            L           K           J           H           G           F           E          D           C           B           A
core_unis = [EMPTY_____, EMPTY_____, EMPTY_____, EMPTY_____, F31_______, F31_06S___, F31_______, F31_06S___, F31_______, F31_06S___, F31_______, EMPTY_____, EMPTY_____, EMPTY_____, EMPTY_____, #  1

             EMPTY_____, EMPTY_____, F31_______, F31______a, F31_16____, F16______B, F31_20____, F16______C, F31_20____, F16______B, F31_16____, F31______a, F31_______, EMPTY_____, EMPTY_____, #  2

             EMPTY_____, F31_______, F31_15SE__, F24_16____, F16______d, F24_16____, F16______b, F24_16____, F16______b, F24_16____, F16______c, F24_16____, F31_15SW__, F31_______, EMPTY_____, #  3

             EMPTY_____, F31______a, F24_16____, F24______D, F24_16____, F16_______, F24_12____, F16______e, F24_12____, F16_______, F24_16____, F24______D, F24_16____, F31______a, EMPTY_____, #  4

             F31_______, F31_16____, F16______c, F24_16____, F16_______, F24_12____, F16_______, F24_12____, F16_______, F24_12____, F16_______, F24_16____, F16______d, F31_16____, F31_______, #  5

             F31_06E___, F16______B, F24_16____, F16_______, F24_12____, F16______C, F24_12____, F16______A, F24_12____, F16______C, F24_12____, F16_______, F24_16____, F16______B, F31_06W___, #  6

             F31_______, F31_20____, F16______b, F24_12____, F16_______, F24_12____, F16_______, F24_16____, F16_______, F24_12____, F16_______, F24_12____, F16______b, F31_20____, F31_______, #  7

             F31_06E___, F16______C, F24_16____, F16______e, F24_12____, F16______A, F24_16____, F16______D, F24_16____, F16______A, F24_12____, F16______e, F24_16____, F16______C, F31_06W___, #  8

             F31_______, F31_20____, F16______b, F24_12____, F16_______, F24_12____, F16_______, F24_16____, F16_______, F24_12____, F16_______, F24_12____, F16______b, F31_20____, F31_______, #  9

             F31_06E___, F16______B, F24_16____, F16_______, F24_12____, F16______C, F24_12____, F16______A, F24_12____, F16______C, F24_12____, F16_______, F24_16____, F16______B, F31_06W___, # 10

             F31_______, F31_16____, F16______d, F24_16____, F16_______, F24_12____, F16_______, F24_12____, F16_______, F24_12____, F16_______, F24_16____, F16______c, F31_16____, F31_______, # 11

             EMPTY_____, F31______a, F24_16____, F24______D, F24_16____, F16_______, F24_12____, F16______e, F24_12____, F16_______, F24_16____, F24______D, F24_16____, F31______a, EMPTY_____, # 12

             EMPTY_____, F31_______, F31_15NE__, F24_16____, F16______c, F24_16____, F16______b, F24_16____, F16______b, F24_16____, F16______d, F24_16____, F31_15NW__, F31_______, EMPTY_____, # 13

             EMPTY_____, EMPTY_____, F31_______, F31______a, F31_16____, F16______B, F31_20____, F16______C, F31_20____, F16______B, F31_16____, F31______a, F31_______, EMPTY_____, EMPTY_____, # 14

             EMPTY_____, EMPTY_____, EMPTY_____, EMPTY_____, F31_______, F31_06N___, F31_______, F31_06N___, F31_______, F31_06N___, F31_______, EMPTY_____, EMPTY_____, EMPTY_____, EMPTY_____] # 15

universes.append({"id": 20000, "type": "rectlinear", "shape": [15,15,1],
  "pitch": [21.50364, 21.50364, 460.], "outer": 17001,
  "origin": [0., 0., 230.], "name": "Normal Assembly",
  "universes": core_unis})
CORE = 20000


#=============================================================================
# Tallies
tallies = []

tallies.append({})
tallies[-1]["quantity"] = "flux"
tallies[-1]["estimator"] = "track-length"
tallies[-1]["name"] = "flux_beavrs"
tallies[-1]["low"] = [-204.28458, -204.28458,   0.]
tallies[-1]["hi"] =  [ 204.28458,  204.28458, 460.]
tallies[-1]["shape"] = [323,323,50]
tallies[-1]['energy-bounds'] = [1.E-11,0.625E-6, 20.]

tallies.append({})
tallies[-1]["quantity"] = "flux"
tallies[-1]["estimator"] = "track-length"
tallies[-1]["name"] = "flux_spectrum_beavrs"
tallies[-1]["low"] = [-190., -190., 0.]
tallies[-1]["hi"] = [190., 190., 460.]
tallies[-1]["shape"] = [1,1,1]
Ebounds = []
Ebounds_array = list(np.logspace(np.log10(1.E-11), np.log10(20.), 2000))
for Ebound in Ebounds_array:
  Ebounds.append(float(Ebound))
tallies[-1]['energy-bounds'] = Ebounds

tallies.append({})
tallies[-1]["quantity"] = "fission"
tallies[-1]["estimator"] = "track-length"
tallies[-1]["name"] = "fission_beavrs"
tallies[-1]["low"] = [-161.2773, -161.2773, 36.748]
tallies[-1]["hi"] =  [ 161.2773,  161.2773, 402.508]
tallies[-1]["shape"] = [255,255,1]
tallies[-1]['energy-bounds'] = [1.E-11, 20.]

#=============================================================================
# Sources
sources = []

sources.append({})
sources[0]["spatial"] = {"type": "box", "low": [-160., -160.,  50.],
                                         "hi": [ 160.,  160., 400.]}
sources[0]["direction"] = {"type": "isotropic"}
sources[0]["energy"] = {"type": "watt", "a": 0.977, "b": 2.546} # U235 Watt spectrum from MCNP manual
sources[0]["fissile-only"] = True
sources[0]["weight"] = 1.

#=============================================================================
# Entropy
entropy = {"low": [-161.2773, -161.2773,  35.],
           "hi" : [ 161.2773,  161.2773, 410.],
           "shape": [8, 8, 8]}

#=============================================================================
# Settings
settings = {}
settings['nparticles'] = 100000
settings['ngenerations'] = 5300
settings['nignored'] = 300
settings['simulation'] = 'k-eigenvalue'
settings['transport'] = 'surface-tracking'
settings['dbrc-nuclides'] = ['U234', 'U235', 'U238']
#settings['max-run-time'] = 2875 # 48hr 55min max run time

#=============================================================================
# Plots
plots = []

plots.append({"type": "slice", "basis": "xy", "resolution": [2500, 2500],
  "origin": [0., 0., 330.],
  "dimensions": [2.*241.4, 2.*241.4],
  "color": "material", "name": "core"})

plots.append({"type": "slice", "basis": "xz", "resolution": [2500, 2500],
  "origin": [0., 0., 230.],
  "dimensions": [2.*241.4, 461.],
  "color": "material", "name": "core_long"})

#=============================================================================
# Write Chenille input file
beavrs = {}
beavrs['materials'] = materials
beavrs['surfaces'] = surfaces
beavrs['cells'] = cells
beavrs['universes'] = universes
beavrs['root-universe'] = CORE
beavrs['sources'] = sources
beavrs['tallies'] = tallies
beavrs['entropy'] = entropy
beavrs['settings'] = settings
beavrs['plots'] = plots

with open('beavrs.yaml', 'w') as input_file:
  yaml.dump(beavrs, input_file)
