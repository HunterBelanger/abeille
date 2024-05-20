from abeille import *
import numpy as np

#===============================================================================
# Materials
UO2 = Material("UO2")
UO2.add_nuclide("O16", 4.70902E-02)
UO2.add_nuclide("U235", 4.30565E-04)
UO2.add_nuclide("U238", 2.31145E-02)

Umetal = Material("UMetal")
Umetal.add_nuclide("U235", 4.53160E-04)
Umetal.add_nuclide("U238", 4.68003E-02)

UO2_Clad = Material("UO2 Cladding")
UO2_Clad.add_nuclide("Al27", 5.00614E-02)

Umetal_Clad = Material("UMetal Cladding")
Umetal_Clad.add_nuclide("Al27", 5.17799E-02)

Moderator = Material("Moderator")
Moderator.add_nuclide("H1_H2O", 6.67578E-02)
Moderator.add_nuclide("O16", 3.33789E-02)

Base_Plate = Material("Base Plate")
Base_Plate.add_nuclide("Al27", 6.02611E-02)

Cadmium = Material("Cadmium")
Cadmium.add_nuclide("Cd106", 0.01245*4.63334E-02)
Cadmium.add_nuclide("Cd108", 0.00888*4.63334E-02)
Cadmium.add_nuclide("Cd110", 0.12470*4.63334E-02)
Cadmium.add_nuclide("Cd111", 0.12795*4.63334E-02)
Cadmium.add_nuclide("Cd112", 0.24109*4.63334E-02)
Cadmium.add_nuclide("Cd113", 0.12227*4.63334E-02)
Cadmium.add_nuclide("Cd114", 0.28754*4.63334E-02)
Cadmium.add_nuclide("Cd116", 0.07512*4.63334E-02)

Filler_Gas = Material("Filler Gas")
Filler_Gas.add_nuclide("He4", 1.6400E-04)

Air = Material("Air")
Air.add_nuclide("N14", 4.1400E-05)
Air.add_nuclide("O16", 9.0700E-06)

Void = Material("Void")
Void.add_nuclide("He4", 1.E-16)

#===============================================================================
# Geometry

# Shared Pin Cell Surfaces
Pin_BG_Bot = ZPlane(-.55 - 57.3)
Pin_BG_Cad_Bot = ZPlane(-0.05 - 57.3)
Pin_BG_Cad_Top = ZPlane(0.0 - 57.3)
Pin_BG_Top = ZPlane(0.5 - 57.3)
Pin_TG_Bot = ZPlane(100.5 - 57.3)
Pin_TG_Cad_Bot = ZPlane(101.0 - 57.3)
Pin_TG_Cad_Top = ZPlane(101.05 - 57.3)
Pin_TG_Top = ZPlane(102.55 - 57.3)
Pin_Water_Top = ZPlane(96.51 - 57.3) # Measures from bottom of active fuel

# UO2 Pin Surfaces
UO2_Fuel_Rad = ZCylinder(0.526)
UO2_Clad_Rad = ZCylinder(0.630)
UO2_Fuel_Bot = ZPlane(  0.0 - 57.3)
UO2_Fuel_Top = ZPlane(100.0 - 57.3)
UO2_Cap_Bot  = ZPlane(100.5 - 57.3)
UO2_Cap_Top  = ZPlane(117.3 - 57.3)

# U-Metal Pin Surfaces
Umetal_Fuel_Rad = ZCylinder(0.8500)
Umetal_Clad_Rad = ZCylinder(0.965)
Umetal_Fuel_Bot = ZPlane(  0.00 - 57.3)
Umetal_Fuel_Top = ZPlane(100.00 - 57.3)
Umetal_Cap_Bot  = ZPlane(101.47 - 57.3)
Umetal_Cap_Top  = ZPlane(117.30 - 57.3)

# Bottom Plate
BP_N  = YPlane(36.)
BP_NE = Plane(1., 1., 0., 54.)
BP_E  = XPlane(36.)
BP_SE = Plane(-1., 1., 0., -54.)
BP_S  = YPlane(-36.)
BP_SW = Plane(1., 1., 0., -54.)
BP_W  = XPlane(-36.)
BP_NW = Plane(-1., 1., 0., 54.)
BP_Top = ZPlane(-2.7)
BP_Bot = ZPlane(-5.7) # Also the reactor tank bottom
BP_Hexagon = -BP_N & +BP_S & +BP_W & -BP_E & -BP_NE & +BP_SW & -BP_NW & +BP_SE
Bottom_Plate = Cell(BP_Hexagon & +BP_Bot & -BP_Top, Base_Plate)

# Bottom Grid
BG_N  = YPlane(38.)
BG_NE = Plane(1., 1., 0., 57.)
BG_E  = XPlane(38.)
BG_SE = Plane(-1., 1., 0., -57.)
BG_S  = YPlane(-38.)
BG_SW = Plane(1., 1., 0., -57.)
BG_W  = XPlane(-38.)
BG_NW = Plane(-1., 1., 0., 57.)
BG_Bot = ZPlane(-.55)
BG_Cad_Bot = ZPlane(-0.05)
BG_Cad_Top = ZPlane(0.0)
BG_Top = ZPlane(0.5)
BG_Hexagon = -BG_N & +BG_S & +BG_W & -BG_E & -BG_NE & +BG_SW & -BG_NW & +BG_SE
BG_Lower = Cell(BG_Hexagon & +BG_Bot & -BG_Cad_Bot, Base_Plate)
BG_Cadmium = Cell(BG_Hexagon & +BG_Cad_Bot & -BG_Cad_Top, Cadmium)
BG_Upper = Cell(BG_Hexagon & +BG_Cad_Top & -BG_Top, Base_Plate)

# Top Grid
TG_N  = YPlane(42.)
TG_NE = Plane(1., 1., 0., 63.)
TG_E  = XPlane(42.)
TG_SE = Plane(-1., 1., 0., -63.)
TG_S  = YPlane(-42.)
TG_SW = Plane(1., 1., 0., -63.)
TG_W  = XPlane(-42.)
TG_NW = Plane(-1., 1., 0., 63.)
TG_Bot = ZPlane(100.5)
TG_Cad_Bot = ZPlane(101.0)
TG_Cad_Top = ZPlane(101.05)
TG_Top = ZPlane(102.55)
TG_Hexagon = -TG_N & +TG_S & +TG_W & -TG_E & -TG_NE & +TG_SW & -TG_NW & +TG_SE
TG_Lower = Cell(TG_Hexagon & +TG_Bot & -TG_Cad_Bot, Base_Plate)
TG_Cadmium = Cell(TG_Hexagon & +TG_Cad_Bot & -TG_Cad_Top, Cadmium)
TG_Upper = Cell(TG_Hexagon & +TG_Cad_Top & -TG_Top, Base_Plate)

# Other Needed surfaces
Tank_iRad = ZCylinder(65.)
Tank_oRad = ZCylinder(66.2, boundary_type='vacuum')
Water_Top = ZPlane(96.51) # Measures from bottom of active fuel
Tank_Top = ZPlane(102.3)
Air_Top = ZPlane(120.0, boundary_type='vacuum')
Reactor_Bottom = ZPlane(-9.7, boundary_type='vacuum')

Tank_Water = Cell(-Tank_iRad & +BP_Bot & -Water_Top & ~Bottom_Plate.region &
                  ~BG_Lower.region & ~BG_Cadmium.region & ~BG_Upper.region,
                  Moderator, 'tank water')

Tank_Air = Cell(-Tank_iRad & +Water_Top & -Air_Top &
                ~TG_Lower.region & ~TG_Cadmium.region & ~TG_Upper.region,
                Air, 'tank air')

Tank_Wall = Cell(+Tank_iRad & -Tank_oRad & +Reactor_Bottom & -Tank_Top, Base_Plate)

Air_Tank_Wall = Cell(-Tank_oRad & +Tank_iRad & +Tank_Top & -Air_Top, Air)

Tank_Bottom = Cell(-Tank_oRad & -BP_Bot & +Reactor_Bottom, Base_Plate)

Tank_Universe = CellUniverse([Bottom_Plate, BG_Lower, BG_Cadmium, BG_Upper,
                              TG_Lower, TG_Cadmium, TG_Upper, Tank_Water,
                              Tank_Air, Air_Tank_Wall, Tank_Bottom, Tank_Wall])

# UO2 Pin Cell
UO2_Al_Bot = Cell(-UO2_Clad_Rad & - UO2_Fuel_Bot, UO2_Clad)
UO2_Fuel = Cell(-UO2_Fuel_Rad & +UO2_Fuel_Bot & -UO2_Fuel_Top, UO2)
UO2_He_Gap  = Cell(-UO2_Fuel_Rad & +UO2_Fuel_Top & -UO2_Cap_Bot, Filler_Gas)
UO2_Cladding = Cell(+UO2_Fuel_Rad & -UO2_Clad_Rad, UO2_Clad)
UO2_Al_Top = Cell(-UO2_Clad_Rad & +UO2_Cap_Bot & -UO2_Cap_Top, UO2_Clad)
UO2_Air_Top_Pin = Cell(-UO2_Clad_Rad & +UO2_Cap_Top, Air)
UO2_Low_Water = Cell(+UO2_Clad_Rad & -Pin_BG_Bot, Moderator)
UO2_BG_Bot = Cell(+UO2_Clad_Rad & +Pin_BG_Bot & -Pin_BG_Cad_Bot, Base_Plate)
UO2_BG_Cad = Cell(+UO2_Clad_Rad & +Pin_BG_Cad_Bot & -Pin_BG_Cad_Top, Cadmium)
UO2_BG_Top = Cell(+UO2_Clad_Rad & +Pin_BG_Cad_Top & -Pin_BG_Top, Base_Plate)
UO2_Mod = Cell(+UO2_Clad_Rad & +Pin_BG_Top & -Pin_Water_Top, Moderator)
UO2_Low_Air = Cell(+UO2_Clad_Rad & +Pin_Water_Top & -Pin_TG_Bot, Air)
UO2_TG_Bot = Cell(+UO2_Clad_Rad & +Pin_TG_Bot & -Pin_TG_Cad_Bot, Base_Plate)
UO2_TG_Cad = Cell(+UO2_Clad_Rad & +Pin_TG_Cad_Bot & -Pin_TG_Cad_Top, Cadmium)
UO2_Tg_Top = Cell(+UO2_Clad_Rad & +Pin_TG_Cad_Top & -Pin_TG_Top, Base_Plate)
UO2_Air_Top = Cell(+UO2_Clad_Rad & +Pin_TG_Top, Air)

U = CellUniverse([UO2_Al_Bot, UO2_Fuel, UO2_He_Gap,
                  UO2_Cladding, UO2_Al_Top, UO2_Air_Top_Pin,
                  UO2_Low_Water, UO2_BG_Bot, UO2_BG_Cad,
                  UO2_BG_Top, UO2_Mod, UO2_Low_Air, UO2_TG_Bot,
                  UO2_TG_Cad, UO2_Tg_Top, UO2_Air_Top],
                 "UO2 Cell")


# U-Metal Pin Cell
Umetal_Al_Bot = Cell(-Umetal_Clad_Rad & - Umetal_Fuel_Bot, Umetal_Clad)
Umetal_Fuel = Cell(-Umetal_Fuel_Rad & +Umetal_Fuel_Bot & -Umetal_Fuel_Top, Umetal)
Umetal_He_Gap  = Cell(-Umetal_Fuel_Rad & +Umetal_Fuel_Top & -Umetal_Cap_Bot, Filler_Gas)
Umetal_Cladding = Cell(+Umetal_Fuel_Rad & -Umetal_Clad_Rad, Umetal_Clad)
Umetal_Al_Top = Cell(-Umetal_Clad_Rad & +Umetal_Cap_Bot & -Umetal_Cap_Top, Umetal_Clad)
Umetal_Air_Top_Pin = Cell(-Umetal_Clad_Rad & +Umetal_Cap_Top, Air)
Umetal_Low_Water = Cell(+Umetal_Clad_Rad & -Pin_BG_Bot, Moderator)
Umetal_BG_Bot = Cell(+Umetal_Clad_Rad & +Pin_BG_Bot & -Pin_BG_Cad_Bot, Base_Plate)
Umetal_BG_Cad = Cell(+Umetal_Clad_Rad & +Pin_BG_Cad_Bot & -Pin_BG_Cad_Top, Cadmium)
Umetal_BG_Top = Cell(+Umetal_Clad_Rad & +Pin_BG_Cad_Top & -Pin_BG_Top, Base_Plate)
Umetal_Mod = Cell(+Umetal_Clad_Rad & +Pin_BG_Top & -Pin_Water_Top, Moderator)
Umetal_Low_Air = Cell(+Umetal_Clad_Rad & +Pin_Water_Top & -Pin_TG_Bot, Air)
Umetal_TG_Bot = Cell(+Umetal_Clad_Rad & +Pin_TG_Bot & -Pin_TG_Cad_Bot, Base_Plate)
Umetal_TG_Cad = Cell(+Umetal_Clad_Rad & +Pin_TG_Cad_Bot & -Pin_TG_Cad_Top, Cadmium)
Umetal_Tg_Top = Cell(+Umetal_Clad_Rad & +Pin_TG_Cad_Top & -Pin_TG_Top, Base_Plate)
Umetal_Air_Top = Cell(+Umetal_Clad_Rad & +Pin_TG_Top, Air)

M = CellUniverse([Umetal_Al_Bot, Umetal_Fuel,
                  Umetal_He_Gap, Umetal_Cladding, Umetal_Al_Top,
                  Umetal_Air_Top_Pin, Umetal_Low_Water,
                  Umetal_BG_Bot, Umetal_BG_Cad, Umetal_BG_Top,
                  Umetal_Mod, Umetal_Low_Air, Umetal_TG_Bot,
                  Umetal_TG_Cad, Umetal_Tg_Top, Umetal_Air_Top],
                 "U-Metal Cell")

n = None

Outter_Lattice = RectLattice(shape=(20,20,1), pitch=(2.917, 2.917, 120.),
                             origin=(0.,0.,57.3), outer_universe=Tank_Universe)
Outter_Lattice.universes = [n, n, n, n, n, n, n, M, M, M, M, M, M, n, n, n, n, n, n, n,
                            n, n, n, n, n, M, M, M, M, M, M, M, M, M, M, n, n, n, n, n, 
                            n, n, n, M, M, M, M, M, M, M, M, M, M, M, M, M, M, n, n, n, 
                            n, n, M, M, M, M, M, M, n, n, n, n, M, M, M, M, M, M, n, n, 
                            n, n, M, M, M, M, n, n, n, n, n, n, n, n, M, M, M, M, n, n, 
                            n, M, M, M, M, M, n, n, n, n, n, n, n, n, M, M, M, M, M, n,
                            n, M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M, n,
                            M, M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M, M,
                            M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M,
                            M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M,
                            M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M,
                            M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M,
                            M, M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M, M,
                            n, M, M, M, n, n, n, n, n, n, n, n, n, n, n, n, M, M, M, n,
                            n, M, M, M, M, M, n, n, n, n, n, n, n, n, M, M, M, M, M, n,
                            n, n, M, M, M, M, n, n, n, n, n, n, n, n, M, M, M, M, n, n,
                            n, n, M, M, M, M, M, M, n, n, n, n, M, M, M, M, M, M, n, n, 
                            n, n, n, M, M, M, M, M, M, M, M, M, M, M, M, M, M, n, n, n,
                            n, n, n, n, n, M, M, M, M, M, M, M, M, M, M, n, n, n, n, n,
                            n, n, n, n, n, n, n, M, M, M, M, M, M, n, n, n, n, n, n, n ]


Inner_Lattice = RectLattice(shape=(22, 22, 1), pitch=(1.837, 1.837, 120.),
                            origin=(0., 0.,57.3), outer_universe= Outter_Lattice)
Inner_Lattice.universes = [n, n, n, n, n, n, n, n, U, U, U, U, U, U, n, n, n, n, n, n, n, n,
                           n, n, n, n, n, n, n, n, U, U, U, U, U, U, n, n, n, n, n, n, n, n,
                           n, n, n, n, n, U, U, U, U, U, U, U, U, U, U, U, U, n, n, n, n, n,
                           n, n, n, n, n, U, U, U, U, U, U, U, U, U, U, U, U, n, n, n, n, n,
                           n, n, n, n, n, U, U, U, U, U, U, U, U, U, U, U, U, n, n, n, n, n,
                           n, n, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, n, n,
                           n, n, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, n, n,
                           n, n, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, n, n,
                           U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
                           U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
                           U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
                           U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
                           U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
                           U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
                           n, n, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, n, n,
                           n, n, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, n, n,
                           n, n, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, n, n,
                           n, n, n, n, n, U, U, U, U, U, U, U, U, U, U, U, U, n, n, n, n, n,
                           n, n, n, n, n, U, U, U, U, U, U, U, U, U, U, U, U, n, n, n, n, n,
                           n, n, n, n, n, U, U, U, U, U, U, U, U, U, U, U, U, n, n, n, n, n,
                           n, n, n, n, n, n, n, n, U, U, U, U, U, U, n, n, n, n, n, n, n, n,
                           n, n, n, n, n, n, n, n, U, U, U, U, U, U, n, n, n, n, n, n, n, n ]

#===============================================================================
# Tallies
#tallies = []
#
#tallies.append(Tally(low=Point(-65., -65., 0.), hi=Point(65., 65., 96.5),
#                     shape=(500, 500, 100), energy_bounds=[1.E-11, 0.625E-6, 20.],
#                     quantity='flux', name='flux', estimator='track-length'))
#
#tallies.append(Tally(low=Point(-65., -65., 0.), hi=Point(65., 65., 96.5),
#                     shape=(1, 1, 1), quantity='flux', name='flux-spectrum', estimator='track-length',
#                     energy_bounds=list(np.logspace(np.log10(1.E-11), np.log10(20.), 2000))))

#===============================================================================
# Settings
settings = Settings()

#===============================================================================
# Simulation
sources = [Source(spatial=Box(Point(-29.17, -29.17, 0.05), Point(29.17, 29.17, 96.), fissile_only=True),
                  direction=Isotropic(),
                  energy=Watt(0.977, 2.546),
                  weight=1.)
          ]

entropy = Entropy(Point(-29.17, -29.17, -2.7), Point(29.17, 29.17, 100.), (5,5,5))

simulation = PowerIterator(nparticles=100000, ngenerations=3100, nignored=100, sources=sources)
simulation.entropy = entropy

input = Input(Inner_Lattice, settings, simulation)
input.to_file('crocus.yaml')
