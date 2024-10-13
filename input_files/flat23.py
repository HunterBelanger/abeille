# FLATTOP-23 - U233-MET-FAST-006
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille      : 1.000246 +/- 0.000222
# OpenMC       : 0.99973  +/- 0.00017

from abeille import *

# Materials
fuel = CEMaterial('U233 Fuel')
fuel.add_nuclide('U233', 4.6710E-2)
fuel.add_nuclide('U234', 5.8772E-4)
fuel.add_nuclide('U235', 1.4158E-5)
fuel.add_nuclide('U238', 2.7959E-4)

ref = CEMaterial('U Reflector')
ref.add_nuclide('U235', 3.5050E-4)
ref.add_nuclide('U238', 4.7719E-2)

# Geometry
fuel_rad = Sphere(r=4.2058)
ref_rad = Sphere(r=24.1194, boundary_type='vacuum')

fuel_cell = Cell(-fuel_rad, fuel)
ref_cell = Cell(+fuel_rad & -ref_rad, ref)
root_uni = CellUniverse([fuel_cell, ref_cell])

# Simulation
sources = [Source(Point(0., 0., 0.), Isotropic(), MonoEnergetic(0.7), 1.)]
simulation = KeffPowerIterator(10000, 2100, 100, sources)

# Write input file
inpt = Input(root_uni, simulation)
inpt.to_file('flat23.yaml')