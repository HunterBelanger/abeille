# Jezebel - PU-MET-FAST-001
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille      : 0.999343 +/- 0.000234
# OpenMC       : 0.99983  +/- 0.00015

from abeille import *

# Material
fuel = Material('Pu Metal')
fuel.add_nuclide('Pu239', 3.7047E-2)
fuel.add_nuclide('Pu240', 1.7512E-3)
fuel.add_nuclide('Pu241', 1.1674E-4)
fuel.add_nuclide('Ga69',  8.2663E-4)
fuel.add_nuclide('Ga71',  5.4857E-4)

# Geometry
surf = Sphere(r=6.3849, boundary_type='vacuum')
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