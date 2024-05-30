# Jezebel-233 - U233-MET-FAST-001
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille-dev  : 1.000337 +/- 0.000235
# OpenMC       : 1.00074  +/- 0.00014

from abeille import *

# Material
fuel = CEMaterial('U Metal')
fuel.add_nuclide('U233', 4.6712E-2)
fuel.add_nuclide('U234', 5.9026E-4)
fuel.add_nuclide('U235', 1.4281E-5)
fuel.add_nuclide('U238', 2.8561E-4)

# Geometry
surf = Sphere(r=5.9838, boundary_type='vacuum')
sphr = Cell(-surf, fuel)
root_uni = CellUniverse([sphr])

# Simulation
sources = [Source(Point(0., 0., 0.), Isotropic(), MonoEnergetic(0.7), 1.)]
simulation = PowerIterator(10000, 2100, 100, sources)

# Write input file
inpt = Input(root_uni, simulation)
inpt.to_file('jez233.yaml')