import abeille as ab

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


Fuel16 = ab.CEMaterial(name="Fuel 1.6%")
Fuel16.temperature = 560.
Fuel16.add_nuclide("O16",  4.5897e-02)
Fuel16.add_nuclide("O17",  1.7436e-05)
Fuel16.add_nuclide("O18",  9.2032e-05)
Fuel16.add_nuclide("U234", 3.0131e-06)
Fuel16.add_nuclide("U235", 3.7503e-04)
Fuel16.add_nuclide("U238", 2.2625e-02)

Fuel24 = ab.CEMaterial(name="Fuel 2.4%")
Fuel24.temperature = 560.
Fuel24.add_nuclide("O16",  4.5830e-02)
Fuel24.add_nuclide("O17",  1.7411e-05)
Fuel24.add_nuclide("O18",  9.1898e-05)
Fuel24.add_nuclide("U234", 4.4842e-06)
Fuel24.add_nuclide("U235", 5.5814e-04)
Fuel24.add_nuclide("U238", 2.2407e-02)

Fuel31 = ab.CEMaterial(name="Fuel 3.1%")
Fuel31.temperature = 560.
Fuel31.add_nuclide("O16",  4.5853e-02)
Fuel31.add_nuclide("O17",  1.7420e-05)
Fuel31.add_nuclide("O18",  9.1942e-05)
Fuel31.add_nuclide("U234", 5.7987e-06)
Fuel31.add_nuclide("U235", 7.2175e-04)
Fuel31.add_nuclide("U238", 2.2253e-02)

Air = ab.CEMaterial(name="Air")
Air.temperature = 560.
Air.add_nuclide("Al27", 1.7352e-03)
Air.add_nuclide("B10",  9.6506e-04)
Air.add_nuclide("B11",  3.9189e-03)
Air.add_nuclide("O16",  4.6514e-02)
Air.add_nuclide("O17",  1.7671e-05)
Air.add_nuclide("O18",  9.3268e-05)
Air.add_nuclide("Si28", 1.6926e-02)
Air.add_nuclide("Si29", 8.5944e-04)
Air.add_nuclide("Si30", 5.6654e-04)

AgInCd = ab.CEMaterial(name="Ag-In-Cd")
AgInCd.temperature = 560.
AgInCd.add_nuclide("Ag107", 2.3523e-02)
AgInCd.add_nuclide("Ag109", 2.1854e-02)
AgInCd.add_nuclide("Cd106", 3.3882e-05)
AgInCd.add_nuclide("Cd108", 2.4166e-05)
AgInCd.add_nuclide("Cd110", 3.3936e-04)
AgInCd.add_nuclide("Cd111", 3.4821e-04)
AgInCd.add_nuclide("Cd112", 6.5611e-04)
AgInCd.add_nuclide("Cd113", 3.3275e-04)
AgInCd.add_nuclide("Cd114", 7.8252e-04)
AgInCd.add_nuclide("Cd116", 2.0443e-04)
AgInCd.add_nuclide("In113", 3.4219e-04)
AgInCd.add_nuclide("In115", 7.6511e-03)

B4C = ab.CEMaterial(name="B4C")
B4C.temperature = 560.
B4C.add_nuclide("B10", 1.5206e-02)
B4C.add_nuclide("B11", 6.1514e-02)
B4C.add_nuclide("C12", 1.8972e-02)
B4C.add_nuclide("C13", 2.1252e-04)

He = ab.CEMaterial(name="He")
He.temperature = 560.
He.add_nuclide("He3", 4.8089e-10)
He.add_nuclide("He4", 2.4044e-04)

Inconel = ab.CEMaterial(name="Inconel 718")
Inconel.temperature = 560.
Inconel.add_nuclide("Cr50", 7.8239e-04)
Inconel.add_nuclide("Cr52", 1.5088e-02)
Inconel.add_nuclide("Cr53", 1.7108e-03)
Inconel.add_nuclide("Cr54", 4.2586e-04)
Inconel.add_nuclide("Fe54", 1.4797e-03)
Inconel.add_nuclide("Fe56", 2.3229e-02)
Inconel.add_nuclide("Fe57", 5.3645e-04)
Inconel.add_nuclide("Fe58", 7.1392e-05)
Inconel.add_nuclide("Mn55", 7.8201e-04)
Inconel.add_nuclide("Ni58", 2.9320e-02)
Inconel.add_nuclide("Ni60", 1.1294e-02)
Inconel.add_nuclide("Ni61", 4.9094e-04)
Inconel.add_nuclide("Ni62", 1.5653e-03)
Inconel.add_nuclide("Ni64", 3.9864e-04)
Inconel.add_nuclide("Si28", 5.6757e-04)
Inconel.add_nuclide("Si29", 2.8820e-05)
Inconel.add_nuclide("Si30", 1.8998e-05)

SS304 = ab.CEMaterial(name="Stainless Steel 304")
SS304.temperature = 560.
SS304.add_nuclide("Cr50", 7.6778e-04)
SS304.add_nuclide("Cr52", 1.4806e-02)
SS304.add_nuclide("Cr53", 1.6789e-03)
SS304.add_nuclide("Cr54", 4.1791e-04)
SS304.add_nuclide("Fe54", 3.4620e-03)
SS304.add_nuclide("Fe56", 5.4345e-02)
SS304.add_nuclide("Fe57", 1.2551e-03)
SS304.add_nuclide("Fe58", 1.6703e-04)
SS304.add_nuclide("Mn55", 1.7604e-03)
SS304.add_nuclide("Ni58", 5.6089e-03)
SS304.add_nuclide("Ni60", 2.1605e-03)
SS304.add_nuclide("Ni61", 9.3917e-05)
SS304.add_nuclide("Ni62", 2.9945e-04)
SS304.add_nuclide("Ni64", 7.6261e-05)
SS304.add_nuclide("Si28", 9.5281e-04)
SS304.add_nuclide("Si29", 4.8381e-05)
SS304.add_nuclide("Si30", 3.1893e-05)

Zirc4 = ab.CEMaterial(name="Zircaloy 4")
Zirc4.temperature = 560.
Zirc4.add_nuclide("Cr50",  3.2962e-06)
Zirc4.add_nuclide("Cr52",  6.3564e-05)
Zirc4.add_nuclide("Cr53",  7.2076e-06)
Zirc4.add_nuclide("Cr54",  1.7941e-06)
Zirc4.add_nuclide("Fe54",  8.6698e-06)
Zirc4.add_nuclide("Fe56",  1.3610e-04)
Zirc4.add_nuclide("Fe57",  3.1431e-06)
Zirc4.add_nuclide("Fe58",  4.1829e-07)
Zirc4.add_nuclide("O16",   3.0744e-04)
Zirc4.add_nuclide("O17",   1.1680e-07)
Zirc4.add_nuclide("O18",   6.1648e-07)
Zirc4.add_nuclide("Sn112", 4.6735e-06)
Zirc4.add_nuclide("Sn114", 3.1799e-06)
Zirc4.add_nuclide("Sn115", 1.6381e-06)
Zirc4.add_nuclide("Sn116", 7.0055e-05)
Zirc4.add_nuclide("Sn117", 3.7003e-05)
Zirc4.add_nuclide("Sn118", 1.1669e-04)
Zirc4.add_nuclide("Sn119", 4.1387e-05)
Zirc4.add_nuclide("Sn120", 1.5697e-04)
Zirc4.add_nuclide("Sn122", 2.2308e-05)
Zirc4.add_nuclide("Sn124", 2.7897e-05)
Zirc4.add_nuclide("Zr90",  2.1828e-02)
Zirc4.add_nuclide("Zr91",  4.7601e-03)
Zirc4.add_nuclide("Zr92",  7.2759e-03)
Zirc4.add_nuclide("Zr94",  7.3734e-03)
Zirc4.add_nuclide("Zr96",  1.1879e-03)

# Borated Water
RATIO_LW = (4.9456e-02)/((4.9456e-02)+(7.7035e-06))
RATIO_HW = 1. - RATIO_LW
BH2O = ab.CEMaterial(name="Borated Water")
BH2O.temperature = 560.
BH2O.add_nuclide("B10",     7.9714e-06)
BH2O.add_nuclide("B11",     3.2247e-05)
BH2O.add_nuclide("H1_H2O",  4.9456e-02)
BH2O.add_nuclide("O16",     RATIO_LW*2.4673e-02)
BH2O.add_nuclide("O17",     RATIO_LW*9.3734e-06)
BH2O.add_nuclide("O18",     RATIO_LW*4.9474e-05)
BH2O.add_nuclide("H2_D2O",  7.7035e-06)
BH2O.add_nuclide("O16_D2O", RATIO_HW*2.4673e-02)
BH2O.add_nuclide("O17_D2O", RATIO_HW*9.3734e-06)
BH2O.add_nuclide("O18_D2O", RATIO_HW*4.9474e-05)

# Nozzle / Support Plate Borated Water
RATIO_LW = (6.5512e-02)/((6.5512e-02)+(1.0204e-05))
RATIO_HW = 1. - RATIO_LW
NS_BH2O = ab.CEMaterial(name="Nozzel/Support Plate Borated Water")
NS_BH2O.add_nuclide("B10",     1.0559e-05)
NS_BH2O.add_nuclide("B11",     4.2716e-05)
NS_BH2O.add_nuclide("H1_H2O",  6.5512e-02)
NS_BH2O.add_nuclide("O16",     RATIO_LW*3.2683e-02)
NS_BH2O.add_nuclide("O17",     RATIO_LW*1.2416e-05)
NS_BH2O.add_nuclide("O18",     RATIO_LW*6.5535e-05)
NS_BH2O.add_nuclide("H2_D2O",  1.0204e-05)
NS_BH2O.add_nuclide("O16_D2O", RATIO_HW*3.2683e-02)
NS_BH2O.add_nuclide("O17_D2O", RATIO_HW*1.2416e-05)
NS_BH2O.add_nuclide("O18_D2O", RATIO_HW*6.5535e-05)

# Nozzle / Support Plate Stainless Steel
NS_SS = ab.CEMaterial(name="Nozzel/Support Plate Stainless Steel")
NS_SS.temperature = 560.
NS_SS.add_nuclide("Cr50", 3.5223e-04)
NS_SS.add_nuclide("Cr52", 6.7924e-03)
NS_SS.add_nuclide("Cr53", 7.7020e-04)
NS_SS.add_nuclide("Cr54", 1.9172e-04)
NS_SS.add_nuclide("Fe54", 1.5882e-03)
NS_SS.add_nuclide("Fe56", 2.4931e-02)
NS_SS.add_nuclide("Fe57", 5.7578e-04)
NS_SS.add_nuclide("Fe58", 7.6625e-05)
NS_SS.add_nuclide("Mn55", 8.0762e-04)
NS_SS.add_nuclide("Ni58", 2.5731e-03)
NS_SS.add_nuclide("Ni60", 9.9117e-04)
NS_SS.add_nuclide("Ni61", 4.3085e-05)
NS_SS.add_nuclide("Ni62", 1.3738e-04)
NS_SS.add_nuclide("Ni64", 3.4985e-05)
NS_SS.add_nuclide("Si28", 4.3711e-04)
NS_SS.add_nuclide("Si29", 2.2195e-05)
NS_SS.add_nuclide("Si30", 1.4631e-05)

CarbonSteel = ab.CEMaterial(name="Carbon Steel")
CarbonSteel.temperature = 560.
CarbonSteel.add_nuclide("Al27",  4.3523e-05)
CarbonSteel.add_nuclide("B10",   2.5833e-06)
CarbonSteel.add_nuclide("B11",   1.0450e-05)
CarbonSteel.add_nuclide("C12",   1.0442e-03)
CarbonSteel.add_nuclide("C13",   1.1697e-05)
CarbonSteel.add_nuclide("Ca40",  1.7043e-05)
CarbonSteel.add_nuclide("Ca42",  1.1375e-07)
CarbonSteel.add_nuclide("Ca43",  2.3734e-08)
CarbonSteel.add_nuclide("Ca44",  3.6673e-07)
CarbonSteel.add_nuclide("Ca46",  7.0322e-10)
CarbonSteel.add_nuclide("Ca48",  3.2875e-08)
CarbonSteel.add_nuclide("Cr50",  1.3738e-05)
CarbonSteel.add_nuclide("Cr52",  2.6493e-04)
CarbonSteel.add_nuclide("Cr53",  3.0041e-05)
CarbonSteel.add_nuclide("Cr54",  7.4778e-06)
CarbonSteel.add_nuclide("Cu63",  1.0223e-04)
CarbonSteel.add_nuclide("Cu65",  4.5608e-05)
CarbonSteel.add_nuclide("Fe54",  4.7437e-03)
CarbonSteel.add_nuclide("Fe56",  7.4465e-02)
CarbonSteel.add_nuclide("Fe57",  1.7197e-03)
CarbonSteel.add_nuclide("Fe58",  2.2886e-04)
CarbonSteel.add_nuclide("Mn55",  6.4126e-04)
CarbonSteel.add_nuclide("Mo100", 2.9814e-05)
CarbonSteel.add_nuclide("Mo92",  4.4822e-05)
CarbonSteel.add_nuclide("Mo94",  2.8110e-05)
CarbonSteel.add_nuclide("Mo95",  4.8567e-05)
CarbonSteel.add_nuclide("Mo96",  5.1015e-05)
CarbonSteel.add_nuclide("Mo97",  2.9319e-05)
CarbonSteel.add_nuclide("Mo98",  7.4327e-05)
CarbonSteel.add_nuclide("Nb93",  5.0559e-06)
CarbonSteel.add_nuclide("Ni58",  4.0862e-04)
CarbonSteel.add_nuclide("Ni60",  1.5740e-04)
CarbonSteel.add_nuclide("Ni61",  6.8420e-06)
CarbonSteel.add_nuclide("Ni62",  2.1815e-05)
CarbonSteel.add_nuclide("Ni64",  5.5557e-06)
CarbonSteel.add_nuclide("P31",   3.7913e-05)
CarbonSteel.add_nuclide("S32",   3.4808e-05)
CarbonSteel.add_nuclide("S33",   2.7420e-07)
CarbonSteel.add_nuclide("S34",   1.5368e-06)
CarbonSteel.add_nuclide("S36",   5.3398e-09)
CarbonSteel.add_nuclide("Si28",  6.1702e-04)
CarbonSteel.add_nuclide("Si29",  3.1330e-05)
CarbonSteel.add_nuclide("Si30",  2.0653e-05)
CarbonSteel.add_nuclide("Ti46",  1.2144e-06)
CarbonSteel.add_nuclide("Ti47",  1.0952e-06)
CarbonSteel.add_nuclide("Ti48",  1.0851e-05)
CarbonSteel.add_nuclide("Ti49",  7.9634e-07)
CarbonSteel.add_nuclide("Ti50",  7.6249e-07)
CarbonSteel.add_nuclide("V50",   1.1526e-07)
CarbonSteel.add_nuclide("V51",   4.5989e-05)

# ALL AXIAL SURFACES
LZC = -230. # Lattice Z Correction
highest_extent   = ab.ZPlane(460.0000+LZC, boundary_type='vacuum')
top_upper_nozzel = ab.ZPlane(431.8760+LZC)
bot_upper_nozzel = ab.ZPlane(423.0490+LZC)
top_bpra_rod     = ab.ZPlane(421.5320+LZC)
top_fuel_rod     = ab.ZPlane(419.7040+LZC)
top_fuel_rod_pln = ab.ZPlane(417.1640+LZC)
top_cr_pln       = ab.ZPlane(415.5580+LZC)
g8_top           = ab.ZPlane(415.1640+LZC)
g8_bot           = ab.ZPlane(411.8060+LZC)
bot_cr_pln       = ab.ZPlane(403.7780+LZC)
top_act_fuel     = ab.ZPlane(402.5080+LZC)
top_act_abs      = ab.ZPlane(401.2380+LZC)
cr_step_288      = ab.ZPlane(400.6380+LZC)
g7_top           = ab.ZPlane(364.7250+LZC)
g7_bot           = ab.ZPlane(359.0100+LZC)
g6_top           = ab.ZPlane(312.5280+LZC)
g6_bot           = ab.ZPlane(306.8130+LZC)
g5_top           = ab.ZPlane(260.3310+LZC)
g5_bot           = ab.ZPlane(254.6160+LZC)
g4_top           = ab.ZPlane(208.1340+LZC)
g4_bot           = ab.ZPlane(202.4190+LZC)
g3_top           = ab.ZPlane(155.9370+LZC)
g3_bot           = ab.ZPlane(150.2220+LZC)
g2_top           = ab.ZPlane(103.7400+LZC)
g2_bot           = ab.ZPlane(098.0250+LZC)
bot_lwr_abs      = ab.ZPlane(041.8280+LZC)
bot_act_abs      = ab.ZPlane(040.5580+LZC)
g1_top           = ab.ZPlane(040.5200+LZC)
cr_step_0        = ab.ZPlane(039.9580+LZC)
bot_bpra_rod     = ab.ZPlane(038.6600+LZC)
g1_bot           = ab.ZPlane(037.1621+LZC)
bot_act_fuel     = ab.ZPlane(036.7480+LZC)
bot_fuel_rod     = ab.ZPlane(035.0000+LZC)
bot_suprt_plt    = ab.ZPlane(020.0000+LZC)
lowest_extent    = ab.ZPlane(000.0000+LZC, boundary_type='vacuum')

# Fuel Pin / Cell Surfaces
spring_rad = ab.ZCylinder(0.06459)
fuel_rad   = ab.ZCylinder(0.39218)
he_rad     = ab.ZCylinder(0.40005)
clad_rad   = ab.ZCylinder(0.45720)

inc_grid_xmin = ab.XPlane(-0.61015)
inc_grid_xmax = ab.XPlane( 0.61015)
inc_grid_ymin = ab.YPlane(-0.61015)
inc_grid_ymax = ab.YPlane( 0.61015)

zirc_grid_xmin = ab.XPlane(-0.61049)
zirc_grid_xmax = ab.XPlane( 0.61049)
zirc_grid_ymin = ab.YPlane(-0.61049)
zirc_grid_ymax = ab.YPlane( 0.61049)

# Shared cells for all fuel pins
gap = ab.Cell(+fuel_rad & -he_rad, He)
clad = ab.Cell(+he_rad & -clad_rad, Zirc4)
mod_no_grid = ab.Cell(+clad_rad, BH2O)
mod_inc_grid = ab.Cell(+clad_rad & +inc_grid_xmin & -inc_grid_xmax & +inc_grid_ymin & -inc_grid_ymax, BH2O)
mod_zirc_grid = ab.Cell(+clad_rad & +zirc_grid_xmin & -zirc_grid_xmax & +zirc_grid_ymin & -zirc_grid_ymax, BH2O)
inc_grid = ab.Cell(-inc_grid_xmin | +inc_grid_xmax | -inc_grid_ymin | +inc_grid_ymax, Inconel)
zirc_grid = ab.Cell(-zirc_grid_xmin | +zirc_grid_xmax | -zirc_grid_ymin | +zirc_grid_ymax, Zirc4)

nozzel = ab.CellUniverse([ab.Cell(-clad_rad, NS_BH2O), ab.Cell(+clad_rad, NS_SS)])
zirc_pin = ab.CellUniverse([ab.Cell(-clad_rad, Zirc4), mod_no_grid])
mod = ab.CellUniverse([ab.Cell(fill=BH2O)]) # Infinite Moderator CellUniverse

# Common cells for tops of fuel pins
spring_inc = ab.Cell(-spring_rad, Inconel)
spring_he  = ab.Cell(+spring_rad & -he_rad, He)
spring_no_grid = ab.CellUniverse([spring_inc, spring_he, clad, mod_no_grid])
spring_inc_grid = ab.CellUniverse([spring_inc, spring_he, clad, mod_inc_grid, inc_grid])

top_mod = ab.Cell(+top_upper_nozzel, mod)
top_noz = ab.Cell(+bot_upper_nozzel & -top_upper_nozzel, nozzel)
mid_mod  = ab.Cell(+top_fuel_rod & -bot_upper_nozzel, mod)
zirc_top = ab.Cell(+top_fuel_rod_pln & -top_fuel_rod, zirc_pin)
sprg_top = ab.Cell(-top_fuel_rod_pln & +g8_top, spring_no_grid)
sprg_grd = ab.Cell(-g8_top & +g8_bot, spring_inc_grid)
sprg_bot = ab.Cell(-g8_bot & +top_act_fuel, spring_no_grid)

# Common cells for bottoms of fuel pins
zirc_bot = ab.Cell(-bot_act_fuel & +bot_fuel_rod, zirc_pin)
bot_noz  = ab.Cell(-bot_fuel_rod & +bot_suprt_plt, nozzel)
bot_mod  = ab.Cell(-bot_suprt_plt, mod)

def make_fuel_rod_uni(Fuel):
  fuel = ab.Cell(-fuel_rad, Fuel)

  fp_no_grid = ab.CellUniverse([fuel, gap, clad, mod_no_grid])
  fp_inc_grid = ab.CellUniverse([fuel, gap, clad, mod_inc_grid, inc_grid])
  fp_zirc_grid = ab.CellUniverse([fuel, gap, clad, mod_zirc_grid, zirc_grid])

  fp_ag7   = ab.Cell(-top_act_fuel & +g7_top, fp_no_grid) # Fuel Pin Above Grid 7
  fp_g7    = ab.Cell(-g7_top & +g7_bot, fp_zirc_grid)     # Fuel Pin Grid 7
  fp_ag6   = ab.Cell(-g7_bot & +g6_top, fp_no_grid)       # Fuel Pin Above Grid 6
  fp_g6    = ab.Cell(-g6_top & +g6_bot, fp_zirc_grid)     # Fuel Pin Grid 6
  fp_ag5   = ab.Cell(-g6_bot & +g5_top, fp_no_grid)       # Fuel Pin Above Grid 5
  fp_g5    = ab.Cell(-g5_top & +g5_bot, fp_zirc_grid)     # Fuel Pin Grid 5
  fp_ag4   = ab.Cell(-g5_bot & +g4_top, fp_no_grid)       # Fuel Pin Above Grid 4
  fp_g4    = ab.Cell(-g4_top & +g4_bot, fp_zirc_grid)     # Fuel Pin Grid 4
  fp_ag3   = ab.Cell(-g4_bot & +g3_top, fp_no_grid)       # Fuel Pin Above Grid 3
  fp_g3    = ab.Cell(-g3_top & +g3_bot, fp_zirc_grid)     # Fuel Pin Grid 3
  fp_ag2   = ab.Cell(-g3_bot & +g2_top, fp_no_grid)       # Fuel Pin Above Grid 2
  fp_g2    = ab.Cell(-g2_top & +g2_bot, fp_zirc_grid)     # Fuel Pin Grid 2
  fp_ag1   = ab.Cell(-g2_bot & +g1_top, fp_no_grid)       # Fuel Pin Above Grid 1
  fp_g1    = ab.Cell(-g1_top & +g1_bot, fp_inc_grid)      # Fuel Pin Grid 1
  fp_bg1   = ab.Cell(-g1_bot & +bot_act_fuel, fp_no_grid) # Fuel Pin Below Grid 1
  

  return ab.CellUniverse([top_mod, top_noz, mid_mod, zirc_top, sprg_top, sprg_grd,
                          sprg_bot, fp_ag7, fp_g7, fp_ag6, fp_g6, fp_ag5, fp_g5,
                          fp_ag4, fp_g4, fp_ag3, fp_g3, fp_ag2, fp_g2, fp_ag1,
                          fp_g1, fp_bg1, zirc_bot, bot_noz, bot_mod])

F16U = make_fuel_rod_uni(Fuel16)
F24U = make_fuel_rod_uni(Fuel24)
F31U = make_fuel_rod_uni(Fuel31)

#===============================================================================
# Simulation
sources = [ab.Source(spatial=ab.Box(ab.Point(-29.17, -29.17, 0.05), ab.Point(29.17, 29.17, 96.), fissile_only=True),
                  direction=ab.Isotropic(),
                  energy=ab.Watt(0.977, 2.546),
                  weight=1.)
          ]

simulation = ab.PowerIterator(nparticles=10000, ngenerations=3100, nignored=100, sources=sources)

input = ab.Input(F16U, simulation)
input.to_file('beavrs_pin.yaml')


