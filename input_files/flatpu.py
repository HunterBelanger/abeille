# FLATTOP-Pu - PU-MET-FAST-006
#
#                keff         std
# ENDF/B-VIII.0
#-------------------------------------
# Abeille      : 0.998517 +/- 0.000219
# OpenMC       : 0.99886  +/- 0.00018

from abeille import *

# Materials
fuel = Material('Pu Fuel')
fuel.add_nuclide('Pu239', 3.6697E-2)
fuel.add_nuclide('Pu240', 1.8700E-3)
fuel.add_nuclide('Pu241', 1.1639E-4)
fuel.add_nuclide('Ga69',  8.8692E-5)
fuel.add_nuclide('Ga71',  5.8858E-5)

ref = Material('U Reflector')
ref.add_nuclide('U234', 2.6438E-6)
ref.add_nuclide('U235', 3.4610E-4)
ref.add_nuclide('U238', 4.7721E-2)

# Geometry
fuel_rad = Sphere(r=4.5332)
ref_rad = Sphere(r=24.1420, boundary_type='vacuum')

fuel_cell = Cell(-fuel_rad, fuel)
ref_cell = Cell(+fuel_rad & -ref_rad, ref)
root_uni = CellUniverse([fuel_cell, ref_cell])

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
inpt.to_file('flat23.yaml')