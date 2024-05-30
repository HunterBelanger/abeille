from abeille import *
import numpy as np

# Results from Abeille
# -----------------------------------
# | kcol    = 1.165292 +/- 0.000046 |
# | ktrk    = 1.165348 +/- 0.000061 |
# -----------------------------------

# This is for normalizing the fission production so that the system will be critical
keff = 1.
#keff = 1.165292
norm = 1. / keff

# Same for all materials
decay_constants = [1.247E-02, 2.829E-02, 4.252E-02, 1.330E-01, 2.925E-01, 6.665E-01, 1.635E+00, 3.555E+00]

nu = norm * np.array([2.78145,     2.47443,     2.43383,     2.43380,     2.43380,     2.43380,     2.43380])
UO2 = MGMaterial("UO2")
UO2.color = [255,127,14]
UO2.Et = [1.77949E-01, 3.29805E-01, 4.80388E-01, 5.54367E-01, 3.11801E-01, 3.95168E-01, 5.64406E-01]
UO2.Ea = [8.02480E-03, 3.71740E-03, 2.67690E-02, 9.62360E-02, 3.00200E-02, 1.11260E-01, 2.82780E-01]
UO2.Ef = [7.21206E-03, 8.19301E-04, 6.45320E-03, 1.85648E-02, 1.78084E-02, 8.30348E-02, 2.16004E-01]
UO2.chi = [[5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.,          0.,          0.]]
UO2.Es = [[1.27537E-01, 4.23780E-02, 9.43740E-06, 5.51630E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
          [0.00000E+00, 3.24456E-01, 1.63140E-03, 3.14270E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
          [0.00000E+00, 0.00000E+00, 4.50940E-01, 2.67920E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.52565E-01, 5.56640E-03, 0.00000E+00, 0.00000E+00],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.25250E-04, 2.71401E-01, 1.02550E-02, 1.00210E-08],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.29680E-03, 2.65802E-01, 1.68090E-02],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.54580E-03, 2.73080E-01]]
UO2.delayed_probabilities = np.array([2.13333E-04, 1.04514E-03, 6.03969E-04, 1.33963E-03 , 2.29386E-03, 7.05174E-04 , 6.00381E-04, 2.07736E-04])
UO2.delayed_constants = decay_constants
UO2.delayed_chi = [[0.00075, 0.98512, 0.01413, 0., 0., 0., 0.],
                   [0.03049, 0.96907, 0.00044, 0., 0., 0., 0.],
                   [0.00457, 0.97401, 0.02142, 0., 0., 0., 0.],
                   [0.02002, 0.97271, 0.00727, 0., 0., 0., 0.],
                   [0.05601, 0.93818, 0.00581, 0., 0., 0., 0.],
                   [0.06098, 0.93444, 0.00458, 0., 0., 0., 0.],
                   [0.10635, 0.88298, 0.01067, 0., 0., 0., 0.],
                   [0.09346, 0.90260, 0.00394, 0., 0., 0., 0.]]
P_delayed = np.sum(UO2.delayed_probabilities)
UO2.nu_prompt = (1. - P_delayed) * nu
UO2.nu_delayed = P_delayed * nu
UO2.group_speeds = [2.23466E+09, 5.07347E+08, 3.86595E+07, 5.13931E+06, 1.67734E+06, 7.28603E+05, 2.92902E+05]


nu = norm * np.array([2.85209,     2.89099,     2.85486,     2.86073,     2.85447,     2.86415,     2.86780])
M4 = MGMaterial("MOX 4.3%")
M4.color = [0, 0, 0]
M4.Et = [1.78731E-01, 3.30849E-01, 4.83772E-01, 5.66922E-01, 4.26227E-01, 6.78997E-01, 6.82852E-01]
M4.Ea = [8.43390E-03, 3.75770E-03, 2.79700E-02, 1.04210E-01, 1.39940E-01, 4.09180E-01, 4.09350E-01]
M4.Ef = [7.62704E-03, 8.76898E-04, 5.69835E-03, 2.28872E-02, 1.07635E-02, 2.32757E-01, 2.48968E-01]
M4.chi = [[5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.,          0.,          0.]]
M4.Es = [[1.28876E-01, 4.14130E-02, 8.22900E-06, 5.04050E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 3.25452E-01, 1.63950E-03, 1.59820E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 4.53188E-01, 2.61420E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.57173E-01, 5.53940E-03, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.60460E-04, 2.76814E-01, 9.31270E-03, 9.16560E-09],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.00510E-03, 2.52962E-01, 1.48500E-02],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.49480E-03, 2.65007E-01]]
M4.delayed_probabilities = np.array([7.82484E-05, 6.40534E-04, 2.27884E-04, 5.78624E-04, 9.97539E-04, 4.33265E-04, 3.22355E-04, 1.23882E-04])
M4.delayed_constants = decay_constants
M4.delayed_chi = [[0.00075, 0.98512, 0.01413, 0., 0., 0., 0.],
                  [0.03069, 0.96887, 0.00044, 0., 0., 0., 0.],
                  [0.00607, 0.97276, 0.02117, 0., 0., 0., 0.],
                  [0.01887, 0.97282, 0.00831, 0., 0., 0., 0.],
                  [0.04990, 0.94419, 0.00591, 0., 0., 0., 0.],
                  [0.05524, 0.93984, 0.00492, 0., 0., 0., 0.],
                  [0.10140, 0.88508, 0.01351, 0., 0., 0., 0.],
                  [0.08055, 0.91408, 0.00537, 0., 0., 0., 0.]]
P_delayed = np.sum(M4.delayed_probabilities)
M4.nu_prompt = (1. - P_delayed) * nu
M4.nu_delayed = P_delayed * nu
M4.group_speeds = [2.23473E+09, 5.07114E+08, 3.88385E+07, 5.16295E+06, 1.75719E+06, 7.68973E+05, 2.94764E+05]


nu = norm * np.array([2.88498,     2.91079,     2.86574,     2.87063,     2.86714,     2.86658,     2.87539])
M7 = MGMaterial("MOX 7%")
M7.color = [214,39,40]
M7.Et = [1.81323E-01, 3.34368E-01, 4.93785E-01, 5.91216E-01, 4.74198E-01, 8.33601E-01, 8.53603E-01]
M7.Ea = [9.06570E-03, 4.29670E-03, 3.28810E-02, 1.22030E-01, 1.82980E-01, 5.68460E-01, 5.85210E-01]
M7.Ef = [8.25446E-03, 1.32565E-03, 8.42156E-03, 3.28730E-02, 1.59636E-02, 3.23794E-01, 3.62803E-01]
M7.chi = [[5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.,          0.,          0.]]
M7.Es = [[1.30457E-01, 4.17920E-02, 8.51050E-06, 5.13290E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 3.28428E-01, 1.64360E-03, 2.20170E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 4.58371E-01, 2.53310E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.63709E-01, 5.47660E-03, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.76190E-04, 2.82313E-01, 8.72890E-03, 9.00160E-09],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.27600E-03, 2.49751E-01, 1.31140E-02],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.86450E-03, 2.59529E-01]]
M7.delayed_probabilities = np.array([7.65120E-05, 6.34833E-04, 2.23483E-04, 5.68882E-04, 9.81163E-04, 4.29227E-04, 3.18971E-04, 1.21830E-04])
M7.delayed_constants = decay_constants
M7.delayed_chi = [[0.00075, 0.98512, 0.01413, 0., 0., 0., 0.],
                  [0.03069, 0.96887, 0.00044, 0., 0., 0., 0.],
                  [0.00612, 0.97272, 0.02116, 0., 0., 0., 0.],
                  [0.01883, 0.97283, 0.00834, 0., 0., 0., 0.],
                  [0.04968, 0.94440, 0.00592, 0., 0., 0., 0.],
                  [0.05506, 0.94002, 0.00492, 0., 0., 0., 0.],
                  [0.10115, 0.88527, 0.01358, 0., 0., 0., 0.],
                  [0.08021, 0.91438, 0.00541, 0., 0., 0., 0.]]
P_delayed = np.sum(M7.delayed_probabilities)
M7.nu_prompt = (1. - P_delayed) * nu
M7.nu_delayed = P_delayed * nu
M7.group_speeds = [2.23479E+09, 5.07355E+08, 3.91436E+07, 5.18647E+06, 1.78072E+06, 7.84470E+05, 3.02310E+05]


nu = norm * np.array([2.90426,     2.91795,     2.86986,     2.87491,     2.87175,     2.86752,     2.87808]) 
M8 = MGMaterial("MOX 8.7%")
M8.color = [150,150,150]
M8.Et = [1.83045E-01, 3.36705E-01, 5.00507E-01, 6.06174E-01, 5.02754E-01, 9.21028E-01, 9.55231E-01]
M8.Ea = [9.48620E-03, 4.65560E-03, 3.62400E-02, 1.32720E-01, 2.08400E-01, 6.58700E-01, 6.90170E-01]
M8.Ef = [8.67209E-03, 1.62426E-03, 1.02716E-02, 3.90447E-02, 1.92576E-02, 3.74888E-01, 4.30599E-01]
M8.chi = [[5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.,          0.,          0.]]
M8.Es = [[1.31504E-01, 4.20460E-02, 8.69720E-06, 5.19380E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 3.30403E-01, 1.64630E-03, 2.60060E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 4.61792E-01, 2.47490E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.68021E-01, 5.43300E-03, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.85970E-04, 2.85771E-01, 8.39730E-03, 8.92800E-09],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.39160E-03, 2.47614E-01, 1.23220E-02],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.96810E-03, 2.56093E-01]]
M8.delayed_probabilities = np.array([7.58799E-05, 6.33750E-04, 2.22271E-04, 5.66810E-04, 9.77854E-04, 4.29965E-04, 3.19265E-04, 1.21188E-04])
M8.delayed_constants = decay_constants
M8.delayed_chi = [[0.00075, 0.98512, 0.01413, 0., 0., 0., 0.],
                  [0.03069, 0.96887, 0.00044, 0., 0., 0., 0.],
                  [0.00614, 0.97270, 0.02116, 0., 0., 0., 0.],
                  [0.01880, 0.97284, 0.00836, 0., 0., 0., 0.],
                  [0.04960, 0.94448, 0.00592, 0., 0., 0., 0.],
                  [0.05496, 0.94012, 0.00492, 0., 0., 0., 0.],
                  [0.10101, 0.88540, 0.01359, 0., 0., 0., 0.],
                  [0.08003, 0.91454, 0.00543, 0., 0., 0., 0.]]
P_delayed = np.sum(M8.delayed_probabilities)
M8.nu_prompt = (1. - P_delayed) * nu
M8.nu_delayed = P_delayed * nu
M8.group_speeds = [2.23483E+09, 5.07520E+08, 3.93259E+07, 5.20109E+06, 1.79321E+06, 7.91377E+05, 3.05435E+05]


FC = MGMaterial("Fission Chamber")
FC.color = [44,160,44]
FC.Et = [1.26032E-01, 2.93160E-01, 2.84250E-01, 2.81020E-01, 3.34460E-01, 5.65640E-01, 1.17214E+00]
FC.Ea = [5.11320E-04, 7.58130E-05, 3.16430E-04, 1.16750E-03, 3.39770E-03, 9.18860E-03, 2.32440E-02]
FC.Ef = [4.79002E-09, 5.82564E-09, 4.63719E-07, 5.24406E-06, 1.45390E-07, 7.14972E-07, 2.08041E-06]
FC.nu = norm * np.array([2.76283,     2.46239,     2.43380,     2.43380,     2.43380,     2.43380,     2.43380])
FC.chi = [[5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.,          0.,          0.]]
FC.Es = [[6.61659E-02, 5.90700E-02, 2.83340E-04, 1.46220E-06, 2.06420E-08, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 2.40377E-01, 5.24350E-02, 2.49900E-04, 1.92390E-05, 2.98750E-06, 4.21400E-07],
         [0.00000E+00, 0.00000E+00, 1.83425E-01, 9.22880E-02, 6.93650E-03, 1.07900E-03, 2.05430E-04],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.90769E-02, 1.69990E-01, 2.58600E-02, 4.92560E-03],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 3.73400E-05, 9.97570E-02, 2.06790E-01, 2.44780E-02],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 9.17420E-04, 3.16774E-01, 2.38760E-01],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.97930E-02, 1.09910E+00]]
FC.group_speeds = [2.24885E+09, 5.12300E+08, 3.75477E+07, 5.02783E+06, 1.66563E+06, 6.70396E+05, 2.51392E+05]


GT = MGMaterial("Guide Tube")
GT.color = [240, 240, 240]
GT.Et = [1.26032E-01, 2.93160E-01, 2.84240E-01, 2.80960E-01, 3.34440E-01, 5.65640E-01, 1.17215E+00]
GT.Ea = [5.11320E-04, 7.58010E-05, 3.15720E-04, 1.15820E-03, 3.39750E-03, 9.18780E-03, 2.32420E-02]
GT.Ef = [0.,          0.,          0.,          0.,          0.,          0.,          0.]
GT.Es = [[6.61659E-02, 5.90700E-02, 2.83340E-04, 1.46220E-06, 2.06420E-08, 0.00000E+00, 0.00000E+00],
         [0.00000E+00, 2.40377E-01, 5.24350E-02, 2.49900E-04, 1.92390E-05, 2.98750E-06, 4.21400E-07],
         [0.00000E+00, 0.00000E+00, 1.83297E-01, 9.23970E-02, 6.94460E-03, 1.08030E-03, 2.05670E-04],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.88511E-02, 1.70140E-01, 2.58810E-02, 4.92970E-03],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 3.73330E-05, 9.97372E-02, 2.06790E-01, 2.44780E-02],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 9.17260E-04, 3.16765E-01, 2.38770E-01],
         [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.97920E-02, 1.09912E+00]]
GT.group_speeds = [2.21473E+09, 4.54712E+08, 4.22099E+07, 5.36964E+06, 1.71422E+06, 7.63783E+05, 2.93629E+05]


WTR = MGMaterial("Water")
WTR.color = [31,119,180]
WTR.Et = [1.59206E-01, 4.12970E-01, 5.90310E-01, 5.84350E-01, 7.18000E-01, 1.25445E+00, 2.65038E+00]
WTR.Ea = [6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02]
WTR.Ef = [0.,          0.,          0.,          0.,          0.,          0.,          0.]
WTR.Es = [[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
          [0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
          [0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21570E-03, 6.99913E-01, 5.37320E-01],
          [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]]
WTR.group_speeds = [2.23517E+09, 4.98880E+08, 3.84974E+07, 5.12639E+06, 1.67542E+06, 7.26031E+05, 2.81629E+05]


CR = MGMaterial("Control Rod")
CR.color = [0, 0, 0]
CR.Et = [2.16768E-01, 4.80098E-01, 8.86369E-01, 9.70009E-01, 9.10482E-01, 1.13775E+00, 1.84048E+00]
CR.Ea = [1.70490E-03, 8.36224E-03, 8.37901E-02, 3.97797E-01, 6.98763E-01, 9.29508E-01, 1.17836E+00]
CR.Ef = [0.,          0.,          0.,          0.,          0.,          0.,          0.]
CR. Es = [[1.7056E-01, 4.4401E-02, 9.8367E-05, 1.2779E-07, 0.0000E+00, 0.0000E+00, 0.0000E+00,],
          [0.0000E+00, 4.7105E-01, 6.8548E-04, 3.9140E-10, 0.0000E+00, 0.0000E+00, 0.0000E+00,],
          [0.0000E+00, 0.0000E+00, 8.0186E-01, 7.2013E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,],
          [0.0000E+00, 0.0000E+00, 0.0000E+00, 5.7075E-01, 1.4602E-03, 0.0000E+00, 0.0000E+00,],
          [0.0000E+00, 0.0000E+00, 0.0000E+00, 6.5556E-05, 2.0784E-01, 3.8149E-03, 3.6976E-09,],
          [0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.0243E-03, 2.0247E-01, 4.7529E-03,],
          [0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 3.5304E-03, 6.5860E-01,]]
CR.group_speeds = [2.18553E+09, 4.21522E+08, 8.76487E+07, 7.47375E+06, 2.28533E+06, 1.01738E+06, 4.11374E+05]


pin_rad = 0.54
pin_pitch = 1.26
fuel_len = 128.52
cr_rest_pos = .5*fuel_len
cr_insertion = 0.
geom_height = 21.42 + fuel_len + 21.42
asmbly_pitch = 17 * pin_pitch
core_bottom = 0.

# Surfaces for fuel pins
pin_cyl = ZCylinder(pin_rad)
pin_bot = ZPlane(-0.5*fuel_len)
pin_top = ZPlane(0.5*fuel_len)

# Create fuel pin cells
UO2_pin = Cell(-pin_cyl & +pin_bot & -pin_top, UO2)
M4_pin = Cell(-pin_cyl & +pin_bot & -pin_top, M4)
M7_pin = Cell(-pin_cyl & +pin_bot & -pin_top, M7)
M8_pin = Cell(-pin_cyl & +pin_bot & -pin_top, M8)
FC_pin = Cell(-pin_cyl & +pin_bot & -pin_top, FC)

# Create moderator cells for water around, above, and bellow pins
Mod_pin = Cell(+pin_cyl & +pin_bot & -pin_top, WTR)
Mod_bot = Cell(-pin_bot, WTR)
Mod_top = Cell(+pin_top, WTR)

# Create cells for control rods
cr_bottom = ZPlane(cr_rest_pos - cr_insertion)
GT_pin = Cell(-pin_cyl & +pin_bot & -cr_bottom, GT)
CR_pin = Cell(-pin_cyl & +cr_bottom, CR)
CR_mod = Cell(+pin_cyl & +cr_bottom, WTR)

# Create universes for the pins
U2u = CellUniverse([UO2_pin, Mod_pin, Mod_bot, Mod_top])
M4u = CellUniverse([M4_pin, Mod_pin, Mod_bot, Mod_top])
M7u = CellUniverse([M7_pin, Mod_pin, Mod_bot, Mod_top])
M8u = CellUniverse([M8_pin, Mod_pin, Mod_bot, Mod_top])
FCu = CellUniverse([FC_pin, Mod_pin, Mod_bot, Mod_top])
GTu = CellUniverse([GT_pin, Mod_pin, Mod_bot, CR_pin, CR_mod])

# Create fuel assemblies
UO2_asmbly = RectLattice(shape=(17, 17, 1), pitch=(pin_pitch, pin_pitch, geom_height), origin=[0., 0., 0.])
UO2_asmbly.universes = [U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, GTu, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, GTu, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, FCu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, GTu, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, GTu, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, GTu, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u,
                        U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u, U2u]

MOX_asmbly = RectLattice(shape=(17, 17, 1), pitch=(pin_pitch, pin_pitch, geom_height), origin=[0., 0., 0.])
MOX_asmbly.universes = [M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u,
                        M4u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M4u,
                        M4u, M7u, M7u, M7u, M7u, GTu, M7u, M7u, GTu, M7u, M7u, GTu, M7u, M7u, M7u, M7u, M4u,
                        M4u, M7u, M7u, GTu, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, GTu, M7u, M7u, M4u,
                        M4u, M7u, M7u, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, M7u, M7u, M4u,
                        M4u, M7u, GTu, M8u, M8u, GTu, M8u, M8u, GTu, M8u, M8u, GTu, M8u, M8u, GTu, M7u, M4u,
                        M4u, M7u, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, M7u, M4u,
                        M4u, M7u, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, M7u, M4u,
                        M4u, M7u, GTu, M8u, M8u, GTu, M8u, M8u, FCu, M8u, M8u, GTu, M8u, M8u, GTu, M7u, M4u,
                        M4u, M7u, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, M7u, M4u,
                        M4u, M7u, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, M7u, M4u,
                        M4u, M7u, GTu, M8u, M8u, GTu, M8u, M8u, GTu, M8u, M8u, GTu, M8u, M8u, GTu, M7u, M4u,
                        M4u, M7u, M7u, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, M7u, M7u, M4u,
                        M4u, M7u, M7u, GTu, M7u, M8u, M8u, M8u, M8u, M8u, M8u, M8u, M7u, GTu, M7u, M7u, M4u,
                        M4u, M7u, M7u, M7u, M7u, GTu, M7u, M7u, GTu, M7u, M7u, GTu, M7u, M7u, M7u, M7u, M4u,
                        M4u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M7u, M4u,
                        M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u, M4u]

# Make empty water cell that consumes all space
WTR_asmbly_cell = Cell(fill=WTR)
WTR_asmbly = CellUniverse([WTR_asmbly_cell])

# Make cell with boundary conditions, that will be filled with core
xl = XPlane(-0.5*3*asmbly_pitch, boundary_type="reflective")
xh = XPlane( 0.5*3*asmbly_pitch, boundary_type="vacuum")
yl = YPlane(-0.5*3*asmbly_pitch, boundary_type="reflective")
yh = YPlane( 0.5*3*asmbly_pitch, boundary_type="vacuum")
zl = ZPlane(core_bottom, boundary_type="vacuum")
zh = ZPlane(core_bottom + geom_height, boundary_type="vacuum")
outer = CellUniverse([Cell(+xl & -xh & +yl & -yh & +zl & -zh, WTR)])

# Make the core assembly
core = RectLattice(shape=(3,3,1), pitch=(asmbly_pitch, asmbly_pitch, geom_height), origin=(0., 0., core_bottom + 0.5*geom_height))
core.outer_universe = outer
core.universes = [UO2_asmbly, MOX_asmbly, WTR_asmbly,
                  MOX_asmbly, UO2_asmbly, WTR_asmbly,
                  WTR_asmbly, WTR_asmbly, WTR_asmbly]


sources = [Source(spatial=Box(Point(xl.x0, yl.y0, core_bottom+21.42), Point(xl.x0 + 2.*asmbly_pitch, yl.y0 + 2.*asmbly_pitch, core_bottom+21.42+128.52), fissile_only=True),
                  direction=Isotropic(),
                  energy=MonoEnergetic(0.5),
                  weight=1.)
          ]

entropy = Entropy(Point(xl.x0, yl.y0, core_bottom+21.42), Point(xl.x0 + 2.*asmbly_pitch, yl.y0 + 2.*asmbly_pitch, core_bottom+21.42+128.52), (8,8,8))

simulation = PowerIterator(nparticles=100000, ngenerations=3500, nignored=500, sources=sources)
simulation.entropy = entropy

settings = Settings()
settings.energy_mode = 'multi-group'
settings.ngroups = 7
settings.energy_bounds = [0., 1., 2., 3., 4., 5., 6., 7.]

#input = Input(geom, simulation)
input = Input(core, simulation, settings=settings)
input.to_file('c5g7_td.yaml')
