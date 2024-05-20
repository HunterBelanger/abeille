# Jezebel-233 - U233-MET-FAST-001
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille-dev  : 1.000337 +/- 0.000235
# OpenMC       : 1.00074  +/- 0.00014

from abeille import *

# Material
fuel = Material('U Metal')
fuel.add_nuclide('U233', 4.6712E-2)
fuel.add_nuclide('U234', 5.9026E-4)
fuel.add_nuclide('U235', 1.4281E-5)
fuel.add_nuclide('U238', 2.8561E-4)

# Geometry
surf = Sphere(r=5.9838, boundary_type='vacuum')
sphr = Cell(-surf, fuel)
root_uni = CellUniverse([sphr])

# Sources
sources = [Source(Point(0., 0., 0.), Isotropic(), MonoEnergetic(0.7), 1.)]

# Settings
settings = Settings()
settings.nparticles = 10000
settings.ngenerations = 2100
settings.nignored = 100
settings.use_dbrc = False

# Write input file
inpt = Input(root_uni, sources, settings)
inpt.to_file('jez233.yaml')