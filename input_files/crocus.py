import yaml
import numpy as np

if __name__ == "__main__":
  # Fuel and Cladding Radii
  fuel_UO2_rad = 0.526
  clad_UO2_rad = 0.630
  fuel_met_rad = 0.850
  clad_met_rad = 0.965

  # Fuel Axial Divisions
  H = 96.51
  fuel_bottom = 0.
  fuel_top = fuel_bottom + 100.0
  pin_Al_bottom = fuel_bottom - 2.7
  UO2_He_top = fuel_top + 0.5
  UO2_top_Al_top = UO2_He_top + 16.8
  Met_He_top = fuel_top + 1.47
  Met_top_Al_top = Met_He_top + 1.583
  
  tank_top = 102.3
  air_top = 120.

  # Other Axial Divisions
  top_Cd_top =  101.05 # Top of Cd layer in top grid
  top_Cd_bot =  101.0  # Bottom of the Cd Layer in top grid
  top_Al_top =  102.55 # Top of top grid plate
  top_Al_bot =  100.5  # Bottom of top grid plate

  bot_Cd_top = fuel_bottom # Top of Cd layer in bottom plate (same level as fuel bottom !)
  bot_Cd_bot = fuel_bottom - 0.05 #  Bottom of the Cd Layer
  bot_Al_top = bot_Cd_top + 0.5 # Top of bottom grid plate
  bot_Al_bot = fuel_bottom - 0.55 # Bottom of bottom grid plate

  base_plate_top = fuel_bottom - 2.7
  base_plate_bottom = fuel_bottom - 5.7
  reactor_vessel_bottom = base_plate_bottom - 4.0 # TODO Check this value !!!
  

  #==============================================================================
  # Surfaces
  surfaces = []

  # Pin Surfaces in Cell Universe
  surfaces.append({"id": 1, "type": "zcylinder", "x0": 0., "y0": 0., "r": fuel_UO2_rad, "name": "fuel_UO2_rad"})
  surfaces.append({"id": 2, "type": "zcylinder", "x0": 0., "y0": 0., "r": clad_UO2_rad, "name": "clad_UO2_rad"})
  surfaces.append({"id": 3, "type": "zcylinder", "x0": 0., "y0": 0., "r": fuel_met_rad, "name": "fuel_met_rad"})
  surfaces.append({"id": 4, "type": "zcylinder", "x0": 0., "y0": 0., "r": clad_met_rad, "name": "clad_met_rad"})
  
  surfaces.append({"id": 5, "type": "zplane", "z0": -57.85, "name": "bottom_grid_bottom"})
  surfaces.append({"id": 6, "type": "zplane", "z0": -57.35, "name": "bottom_Cd_bottom"})
  surfaces.append({"id": 7, "type": "zplane", "z0": -57.3, "name": "fuel_bottom and Cd top"})
  surfaces.append({"id": 8, "type": "zplane", "z0": -56.8, "name": "bottom_grid_top"})
  surfaces.append({"id": 9, "type": "zplane", "z0": -57.3 + H, "name": "water_top"})
  surfaces.append({"id":10, "type": "zplane", "z0": 42.7, "name": "fuel_top"})
  surfaces.append({"id":11, "type": "zplane", "z0": 43.2, "name": "UO2_He_top"})
  surfaces.append({"id":12, "type": "zplane", "z0": 44.17, "name": "UMet_He_top"})
  surfaces.append({"id":13, "type": "zplane", "z0": 43.2, "name": "top_grid_bottom"})
  surfaces.append({"id":14, "type": "zplane", "z0": 43.7, "name": "top_Cd_bottom"})
  surfaces.append({"id":15, "type": "zplane", "z0": 43.75, "name": "top_Cd_top"})
  surfaces.append({"id":16, "type": "zplane", "z0": 45.25, "name": "top_grid_top"})

  # Bottom Plate
  surfaces.append({"id": 20, "type": "yplane", "y0": 36., "name": "N-bot-plate"})
  surfaces.append({"id": 21, "type": "plane", "A": 1., "B": 1., "C": 0., "D": 54., "name": "NE-bot-plate"})
  surfaces.append({"id": 22, "type": "xplane", "x0": 36., "name": "E-bot-plate"})
  surfaces.append({"id": 23, "type": "plane", "A": -1., "B": 1., "C": 0., "D": -54., "name": "SE-bot-plate"})
  surfaces.append({"id": 24, "type": "yplane", "y0": -36., "name": "S-bot-plate"})
  surfaces.append({"id": 25, "type": "plane", "A": 1., "B": 1., "C": 0., "D": -54., "name": "SW-bot-plate"})
  surfaces.append({"id": 26, "type": "xplane", "x0": -36., "name": "W-bot-plate"})
  surfaces.append({"id": 27, "type": "plane", "A": -1., "B": 1., "C": 0., "D": 54., "name": "NW-bot-plate"})
  surfaces.append({"id": 28, "type": "zplane", "z0": base_plate_bottom, "name": "base_plate_bottom"})
  surfaces.append({"id": 29, "type": "zplane", "z0": base_plate_top, "name": "base_plate_top"})

  # Bottom Grid
  surfaces.append({"id": 30, "type": "yplane", "y0": 38., "name": "N-bot-grid"})
  surfaces.append({"id": 31, "type": "plane", "A": 1., "B": 1., "C": 0., "D": 57., "name": "NE-bot-grid"})
  surfaces.append({"id": 32, "type": "xplane", "x0": 38., "name": "E-bot-grid"})
  surfaces.append({"id": 33, "type": "plane", "A": -1., "B": 1., "C": 0., "D": -57., "name": "SE-bot-grid"})
  surfaces.append({"id": 34, "type": "yplane", "y0": -38., "name": "S-bot-grid"})
  surfaces.append({"id": 35, "type": "plane", "A": 1., "B": 1., "C": 0., "D": -57., "name": "SW-bot-grid"})
  surfaces.append({"id": 36, "type": "xplane", "x0": -38., "name": "W-bot-grid"})
  surfaces.append({"id": 37, "type": "plane", "A": -1., "B": 1., "C": 0., "D": 57., "name": "NW-bot-grid"})
  surfaces.append({"id": 38, "type": "zplane", "z0": bot_Al_bot, "name": "bot_Al_bot"})
  surfaces.append({"id": 39, "type": "zplane", "z0": bot_Al_top, "name": "bot_Al_top"})
  surfaces.append({"id": 41, "type": "zplane", "z0": bot_Cd_bot, "name": "bot_Cd_bot"})
  surfaces.append({"id": 42, "type": "zplane", "z0": bot_Cd_top, "name": "bot_Cd_top"})

  # Top Grid
  surfaces.append({"id": 50, "type": "yplane", "y0": 42., "name": "N-top-grid"})
  surfaces.append({"id": 51, "type": "plane", "A": 1., "B": 1., "C": 0., "D": 63., "name": "NE-top-grid"})
  surfaces.append({"id": 52, "type": "xplane", "x0": 42., "name": "E-bot-grid"})
  surfaces.append({"id": 53, "type": "plane", "A": -1., "B": 1., "C": 0., "D": -63., "name": "SE-top-grid"})
  surfaces.append({"id": 54, "type": "yplane", "y0": -42., "name": "S-bot-grid"})
  surfaces.append({"id": 55, "type": "plane", "A": 1., "B": 1., "C": 0., "D": -63., "name": "SW-top-grid"})
  surfaces.append({"id": 56, "type": "xplane", "x0": -42., "name": "W-bot-grid"})
  surfaces.append({"id": 57, "type": "plane", "A": -1., "B": 1., "C": 0., "D": 63., "name": "NW-top-grid"})
  surfaces.append({"id": 58, "type": "zplane", "z0": top_Al_bot, "name": "top_Al_bot"})
  surfaces.append({"id": 59, "type": "zplane", "z0": top_Al_top, "name": "top_Al_top"})
  surfaces.append({"id": 61, "type": "zplane", "z0": top_Cd_bot, "name": "top_Cd_bot"})
  surfaces.append({"id": 62, "type": "zplane", "z0": top_Cd_top, "name": "top_Cd_top"})

  surfaces.append({"id": 100, "type": "zcylinder", "x0": 0., "y0": 0., "r": 65.0, "name": "inner-tank-edge"})
  surfaces.append({"id": 101, "type": "zcylinder", "x0": 0., "y0": 0., "r": 66.2, "boundary": "vacuum", "name": "outer-tank-edge"})
  surfaces.append({"id": 102, "type": "zplane", "z0": H, "name": "water_top"})
  surfaces.append({"id": 103, "type": "zplane", "z0": tank_top, "name": "tank-top"})
  surfaces.append({"id": 104, "type": "zplane", "z0": air_top, "boundary": "vacuum", "name": "air_top"})
  surfaces.append({"id": 105, "type": "zplane", "z0": reactor_vessel_bottom, "boundary": "vacuum", "name": "vessel_bottom"})

  #==============================================================================
  # Materials
  materials = []
  
  Fuel_UO2 = 1
  materials.append({"id": Fuel_UO2, "name": "Fuel UO2", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "O16"   , "fraction": 4.70902E-02},
                                    {"nuclide": "U235"  , "fraction": 4.30565E-04},
                                    {"nuclide": "U238"  , "fraction": 2.31145E-02}]})

  Fuel_Metal = 2
  materials.append({"id": Fuel_Metal, "name": "Fuel Metal", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "U235"  , "fraction": 4.53160E-04},
                                    {"nuclide": "U238"  , "fraction": 4.68003E-02}]})
  
  Clad_UO2 = 3
  materials.append({"id": Clad_UO2, "name": "Clad UO2", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "Al27", "fraction": 5.00614E-02}]})

  Clad_Metal = 4
  materials.append({"id": Clad_Metal, "name": "Clad Metal", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "Al27", "fraction": 5.17799E-02}]})
  
  Moderator = 5
  materials.append({"id": Moderator, "name": "Moderator", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "H1_H2O", "fraction": 6.67578E-02},
                                    {"nuclide": "O16"   , "fraction": 3.33789E-02}]})
  Base_Plate = 6
  materials.append({"id": Base_Plate, "name": "Base Plate", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "Al27", "fraction": 6.02611E-02}]})

  Cadmium = 7
  materials.append({"id": Cadmium, "name": "Cadmium", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "Cd106", "fraction": 0.01245*4.63334E-02},
                                    {"nuclide": "Cd108", "fraction": 0.00888*4.63334E-02},
                                    {"nuclide": "Cd110", "fraction": 0.1247*4.63334E-02},
                                    {"nuclide": "Cd111", "fraction": 0.12795*4.63334E-02},
                                    {"nuclide": "Cd112", "fraction": 0.24109*4.63334E-02},
                                    {"nuclide": "Cd113", "fraction": 0.12227*4.63334E-02},
                                    {"nuclide": "Cd114", "fraction": 0.28754*4.63334E-02},
                                    {"nuclide": "Cd116", "fraction": 0.07512*4.63334E-02}]})

  Filler_Gas = 8
  materials.append({"id": Filler_Gas, "name": "Filler Gas", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "He4", "fraction": 1.6400E-04}]})

  Air = 9
  materials.append({"id": Air, "name": "Air", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "N14", "fraction": 4.1400E-05},
                                    {"nuclide": "O16", "fraction": 9.0700E-06}]})

  Void = 10
  materials.append({"id": Void, "name": "Void", "temperature": 293.6,
                    "density-units": "sum", "fractions": "atoms",
                    "composition": [{"nuclide": "He4", "fraction": 1.E-16}]})

  #=============================================================================
  # Cells
  cells = []
  
  # Lower Aluminum plate
  cells.append({"id": 10, "region": "-20 & -21 & -22 & +23 & +24 & +25 & +26 & -27 & +28 & -29",
                "material": Base_Plate})
  
  # Lower portion of Bottom Al grid plate
  cells.append({"id": 11, "region": "-30 & -31 & -32 & +33 & +34 & +35 & +36 & -37 & +38 & -41",
                "material": Base_Plate})
  
  # Lower Cd layer
  cells.append({"id": 12, "region": "-30 & -31 & -32 & +33 & +34 & +35 & +36 & -37 & +41 & -42",
                "material": Cadmium})
  
  # Upper portion of Bottom Al grid plate
  cells.append({"id": 13, "region": "-30 & -31 & -32 & +33 & +34 & +35 & +36 & -37 & +42 & -39",
                "material": Base_Plate})

  # Lower portion of Top Al Grid plate
  cells.append({"id": 14, "region": "-50 & -51 & -52 & +53 & +54 & +55 & +56 & -57 & +58 & -61",
               "material": Base_Plate})

  # Cadmium of Top Al Grid plate
  cells.append({"id": 15, "region": "-50 & -51 & -52 & +53 & +54 & +55 & +56 & -57 & +61 & -62",
               "material": Cadmium})

  # Upper portion of Top Al Grid plate
  cells.append({"id": 16, "region": "-50 & -51 & -52 & +53 & +54 & +55 & +56 & -57 & +62 & -59",
               "material": Base_Plate})

  # Tank water
  cells.append({"id": 17, "region" : "-100 & +28 & -102 & ~ "
                                     "(-20 & -21 & -22 & +23 & +24 & +25 & +26 & -27 & +28 & -29) & ~ "
                                     "(-30 & -31 & -32 & +33 & +34 & +35 & +36 & -37 & +38 & -39)",
                "material": Moderator})

  # Tank Wall
  cells.append({"id": 18, "region": "+100 & -101 & +28 & -103", "material": Base_Plate})

  # Air in tank and above
  cells.append({"id": 19, "region": "-100 & +102 & -104 & ~ "
                                    "(-50 & -51 & -52 & +53 & +54 & +55 & +56 & -57 & +58 & -59)",
                "material": Air})
  
  # Air above tank wall
  cells.append({"id": 20, "region": "+100 & -101 & +103 & -104", "material": Air})

  # Tank bottom
  cells.append({"id": 21, "region": "-101 & +105 & -28", "material": Base_Plate})

  #-----------------------------------------------------------------------------
  # UO2 Pin Cell
  
  # Al bottom
  cells.append({"id": 30, "region": "-2 & -7", "material": Clad_UO2})
  # Fuel
  cells.append({"id": 31, "region": "-1 & +7 & -10", "material": Fuel_UO2})
  # Clad
  cells.append({"id": 32, "region": "+1  & -2 & +7 & -11", "material": Clad_UO2})
  # He
  cells.append({"id": 33, "region": "-1 & +10 & -11", "material": Filler_Gas})
  # Al top
  cells.append({"id": 34, "region": "-2 & +11", "material": Clad_UO2})

  # Lower Water
  cells.append({"id": 35, "region": "+2 & -5", "material": Moderator})
  # Bottom Al of bottom grid
  cells.append({"id": 36, "region": "+2 & +5 & -6", "material": Base_Plate})
  # Bottom Cd 
  cells.append({"id": 37, "region": "+2 & +6 & -7", "material": Cadmium})
  # Top Al of bottom grid
  cells.append({"id": 38, "region": "+2 & +7 & -8", "material": Base_Plate})
  # Main Water
  cells.append({"id": 39, "region": "+2 & +8 & -9", "material": Moderator})
  # Air
  cells.append({"id": 40, "region": "+2 & +9 & -13", "material": Air})
  # Bottom Al of top grid
  cells.append({"id": 41, "region": "+2 & +13 & -14", "material": Base_Plate})
  # Cd top grid
  cells.append({"id": 42, "region": "+2 & +14 & -15", "material": Cadmium})
  # Top Al of top grid
  cells.append({"id": 43, "region": "+2 & +15 & -16", "material": Base_Plate})
  # Top Air 
  cells.append({"id": 44, "region": "+2 & +16", "material": Air})

  #-----------------------------------------------------------------------------
  # U Metal Pin Cell
  
  # Al bottom
  cells.append({"id": 50, "region": "-4 & -7", "material": Clad_Metal})
  # Fuel
  cells.append({"id": 51, "region": "-3 & +7 & -10", "material": Fuel_Metal})
  # Clad
  cells.append({"id": 52, "region": "+3  & -4 & +7 & -12", "material": Clad_Metal})
  # He
  cells.append({"id": 53, "region": "-3 & +10 & -12", "material": Filler_Gas})
  # Al top
  cells.append({"id": 54, "region": "-4 & +12", "material": Clad_Metal})

  # Lower Water
  cells.append({"id": 55, "region": "+4 & -5", "material": Moderator})
  # Bottom Al of bottom grid
  cells.append({"id": 56, "region": "+4 & +5 & -6", "material": Base_Plate})
  # Bottom Cd 
  cells.append({"id": 57, "region": "+4 & +6 & -7", "material": Cadmium})
  # Top Al of bottom grid
  cells.append({"id": 58, "region": "+4 & +7 & -8", "material": Base_Plate})
  # Main Water
  cells.append({"id": 59, "region": "+4 & +8 & -9", "material": Moderator})
  # Air
  cells.append({"id": 60, "region": "+4 & +9 & -13", "material": Air})
  # Bottom Al of top grid
  cells.append({"id": 61, "region": "+4 & +13 & -14", "material": Base_Plate})
  # Cd top grid
  cells.append({"id": 62, "region": "+4 & +14 & -15", "material": Cadmium})
  # Top Al of top grid
  cells.append({"id": 63, "region": "+4 & +15 & -16", "material": Base_Plate})
  # Top Air 
  cells.append({"id": 64, "region": "+4 & +16", "material": Air})

  #==============================================================================
  # Universes
  LC = +57.3
  universes = []
  
  # UO2 Pin Cell Universe
  universes.append({"id": 1, "cells": [30,31,32,33,34,35,36,37,38,39,40,41,42,43,44]})

  # U Metal Pin Cell Universe
  universes.append({"id": 2, "cells": [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64]})

  # UO2 Lattice
  universes.append({"id": 3, "type": "rectlinear", "shape": [22, 22, 1],
    "pitch": [1.837, 1.837, 120.], "outer": 4,
    "origin": [0., 0., LC], "name": "UO2 Lattice",
    "universes": [-1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,

                  -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,
                  
                  -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,

                  -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,

                   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  
                   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,

                   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,

                   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  
                   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,

                   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,

                  -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,
                  
                  -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,

                  -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,
                  
                  -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,
                  
                  -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1 ]})

  # U Metal Lattice
  universes.append({"id": 4, "type": "rectlinear", "shape": [20, 20, 1],
    "pitch": [2.917, 2.917, 120.], "outer": 5,
    "origin": [0., 0., LC], "name": "U Metal Lattice",
    "universes": [-1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, -1, -1, -1, -1, -1, 
                  
                  -1, -1, -1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, -1, -1, -1, 

                  -1, -1,  2,  2,  2,  2,  2,  2, -1, -1, -1, -1,  2,  2,  2,  2,  2,  2, -1, -1, 
                  
                  -1, -1,  2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2, -1, -1, 

                  -1,  2,  2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2,  2, -1,

                  -1,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2, -1,

                   2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2,

                   2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,
                  
                   2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,

                   2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,
                  
                   2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,
                   
                   2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2,

                  -1,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2, -1,
                  
                  -1,  2,  2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2,  2, -1,

                  -1, -1,  2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2, -1, -1,
                  
                  -1, -1,  2,  2,  2,  2,  2,  2, -1, -1, -1, -1,  2,  2,  2,  2,  2,  2, -1, -1, 
                  
                  -1, -1, -1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, -1, -1, -1,

                  -1, -1, -1, -1, -1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, -1, -1, -1, -1, -1,

                  -1, -1, -1, -1, -1, -1, -1,  2,  2,  2,  2,  2,  2, -1, -1, -1, -1, -1, -1, -1 ]})

  # Base plate, lower grid, and water, and tank
  universes.append({"id": 5, "cells": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]})

  root_universe = 3

  #=============================================================================
  # Sources
  sources = []
  
  sources.append({})
  sources[0]["spatial"] = {"type": "box", "low": [-29.17, -29.17, 0.05],
                                           "hi": [29.17, 29.17, H-0.05]}
  sources[0]["direction"] = {"type": "isotropic"}
  sources[0]["energy"] = {"type": "watt", "a": 0.977, "b": 2.546} # U235 Watt spectrum from MCNP manual
  sources[0]["fissile-only"] = True
  sources[0]["weight"] = 1.

  #=============================================================================
  # Tallies
  tallies = []

  tallies.append({})
  tallies[0]['shape'] = [500, 500, 100]
  tallies[0]['low'] = [-65., -65., 0.]
  tallies[0]['hi'] = [65., 65., H]
  tallies[0]['name'] = "flux"
  tallies[0]['quantity'] = "flux"
  tallies[0]['estimator'] = "track-length"
  tallies[0]['energy-bounds'] = [1.E-11, 0.625E-6, 20.]

  tallies.append({})
  tallies[1]['shape'] = [1, 1, 1]
  tallies[1]['low'] = [-65., -65., 0.]
  tallies[1]['hi'] = [65., 65., H]
  tallies[1]['name'] = "flux_spectrum"
  tallies[1]['quantity'] = "flux"
  tallies[1]['estimator'] = "track-length"
  Ebounds = []
  Ebounds_array = list(np.logspace(np.log10(1.E-11), np.log10(20.), 2000))
  for Ebound in Ebounds_array:
    Ebounds.append(float(Ebound))
  tallies[1]['energy-bounds'] = Ebounds


  #=============================================================================
  # Entropy
  entropy = {}
  entropy['shape'] = [20,20,50]
  entropy['low'] = [-29.17, -29.17, -2.7]
  entropy['hi'] = [29.17, 29.17, 100.]


  #=============================================================================
  # Settings
  settings = {}
  settings['nparticles'] = 10000
  settings['ngenerations'] = 3200
  settings['nignored'] = 200
  settings['simulation'] = 'k-eigenvalue'
  settings['transport'] = 'surface-tracking'


  #=============================================================================
  # Plots
  plots = []
  
  plots.append({"type": "slice", "basis": "xz", "resolution": [5000, 5000],
    "origin": [0., 0.5*1.837, base_plate_bottom - 0.5*(base_plate_bottom - air_top)], "dimensions": [1.1*2*66.2, 1.1*(base_plate_bottom - air_top)],
    "color": "material", "name": "crocus_core_long"})

  plots.append({"type": "slice", "basis": "xy", "resolution": [5000, 5000],
    "origin": [0., 0., 15.], "dimensions": [133., 133.],
    "color": "material", "name": "crocus_core"})


  #=============================================================================
  # Write input file
  crocus = {}

  crocus['materials'] = materials
  crocus['surfaces'] = surfaces
  crocus['cells'] = cells
  crocus['universes'] = universes
  crocus['root-universe'] = root_universe
  crocus['sources'] = sources
  crocus['tallies'] = tallies
  crocus['entropy'] = entropy
  crocus['settings'] = settings
  crocus['plots'] = plots

  with open('crocus.yaml', 'w') as input_file:
    yaml.dump(crocus, input_file)
