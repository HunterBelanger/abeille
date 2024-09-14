# Godiva - HEU-MET-FAST-001 Case 2
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille      : 0.999988 +/- 0.000208
# OpenMC       : 1.00041  +/- 0.00015

from abeille import *

# Material
fuel = CEMaterial('HEU Metal')
fuel.add_nuclide('U234', 4.9184E-4)
fuel.add_nuclide('U235', 4.4994E-2)
fuel.add_nuclide('U238', 2.4984E-3)

# Geometry
surf = Sphere(r=8.7407, boundary_type='vacuum')
sphr = Cell(-surf, fuel)
root_uni = CellUniverse([sphr])

# Simulation
sources = [Source(Point(0., 0., 0.), Isotropic(), MonoEnergetic(0.7), 1.)]
simulation = PowerIterator(10000, 2100, 100, sources)

# Write input file
inpt = Input(root_uni, simulation)
inpt.to_file('godiva.yaml')