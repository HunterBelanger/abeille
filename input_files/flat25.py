# FLATTOP-25 - HEU-MET-FAST-028
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille      : 1.001045 +/- 0.000205
# OpenMC       : 1.00085  +/- 0.00017

from abeille import *

# Materials
fuel = CEMaterial('U Fuel')
fuel.add_nuclide('U234', 4.8869E-4)
fuel.add_nuclide('U235', 4.4482E-2)
fuel.add_nuclide('U238', 2.7038E-3)

ref = CEMaterial('U Reflector')
ref.add_nuclide('U234', 2.6438E-6)
ref.add_nuclide('U235', 3.4610E-4)
ref.add_nuclide('U238', 4.7721E-2)

# Geometry
fuel_rad = Sphere(r=6.1156)
ref_rad = Sphere(r=24.1242, boundary_type='vacuum')

fuel_cell = Cell(-fuel_rad, fuel)
ref_cell = Cell(+fuel_rad & -ref_rad, ref)
root_uni = CellUniverse([fuel_cell, ref_cell])

# Simulation
sources = [Source(Point(0., 0., 0.), Isotropic(), MonoEnergetic(0.7), 1.)]
simulation = KeffPowerIterator(10000, 2100, 100, sources)

# Write input file
inpt = Input(root_uni, simulation)
inpt.to_file('flat25.yaml')