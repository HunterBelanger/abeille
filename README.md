# MGMC - A Multi-Group Monte Carlo Transport Code
[![arXiv](https://img.shields.io/badge/arXiv-2103.13891-b31b1b.svg?style=flat)](https://arxiv.org/abs/2103.13891)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4585368.svg)](https://doi.org/10.5281/zenodo.4585368)
[![License](https://img.shields.io/badge/License-CeCILL%20v2.1-brightgreen)](http://www.cecill.info)

```
 ███╗   ███╗ ██████╗ ███╗   ███╗ ██████╗
 ████╗ ████║██╔════╝ ████╗ ████║██╔════╝
 ██╔████╔██║██║  ███╗██╔████╔██║██║
 ██║╚██╔╝██║██║   ██║██║╚██╔╝██║██║
 ██║ ╚═╝ ██║╚██████╔╝██║ ╚═╝ ██║╚██████╗
 ╚═╝     ╚═╝ ╚═════╝ ╚═╝     ╚═╝ ╚═════╝

 Multi-Group Monte Carlo Transport Code
```

MGMC is a 3D multi-group Monte Carlo transport code which solves the Boltzmann
neutron transport equation for fixed-source, k-eigenvalue, and neutron noise
problems. 

All problem parameters such as geometry, materials, and simulation settings are
controlled through a single YAML input file. Several examples of such files are
provided in the `input_file` directory. A problem geometry may be constructed
using combinations of surface half-spaces to create cells. Like in other Monte
Carlo codes, Universes and Lattices may be used to construct more complex
geometries such as for fuel assemblies found in nuclear reactors. Material
cross sections may be provided with any number of energy groups, so long as all
materials in the problem have the same number of groups.

Three different particle tracking methods are available to users:
surface-tracking, delta-tracking [1], and carter-tracking [2]. Either of these
may be selected in the settings portion of the input file.  When using
carter-tracking, if the majorant cross section is under-estimated, the particle
population will not be stable unless particle weight cancellation is used. A
mesh may be defined over which 3D regional cancellation is used to annihilate
positive and negative weight.

Tallies for the scalar flux and reaction rates on rectilinear mesh are
available, using a collision estimator or a track-length estimator. The tallies
are saved in `.npy` files, for easy plotting in Python. Plots of the geometry
may also be specified in the input file, and are generated when the program is
run with the `--plot` flag. Shared memory parallelism is implemented with
OpenMP, and is turned on by default. Distributed memory parallelism with MPI
can be turned on at compiled time with the `-DMGMC_USE_MPI=ON` opiton
when running cmake.

In addition to solving standard k-eigenvalue and fixed-source problems, MGMC
is also able to solve neutron noise problems for macroscopic cross section
oscillations and for vibrations of flat surfaces. The basic methods to perform
noise transport were developed by Dr Amélie Rouchon durring her PhD [3,4].

[1] J. Leppänen, “On the use of delta-tracking and the collision flux estimator
in the Serpent 2 Monte Carlo particle transport code,” Ann Nucl Energy, vol.
105, pp. 161–167, 2017, doi:
[10.1016/j.anucene.2017.03.006](https://dx.doi.org/10.1016/j.anucene.2017.03.006).

[2] L. L. Carter, E. D. Cashwell, and W. M. Taylor, “Monte Carlo Sampling with
Continuously Varying Cross Sections Along Flight Paths,” Nucl Sci Eng, vol. 48,
no. 4, pp. 403–411, 1972, doi:
[10.13182/nse72-1](https://dx.doi.org/10.13182/nse72-1). 

[3] A. Rouchon, “Analyse et développement d’outils numériques déterministes et
stochastiques résolvant les équations du bruit neutronique et applications aux
réacteurs thermiques et rapides,” 2016.

[4] A. Rouchon, A. Zoia, and R. Sanchez, “A new Monte Carlo method for neutron
noise calculations in the frequency domain,” Ann Nucl Energy, vol. 102,
pp. 465–475, 2017, doi: 10.1016/j.anucene.2016.11.035. 

## Papers Using MGMC

H. Belanger, D. Mancusi, and A. Zoia, “Variance Reduction Techniques for Monte
Carlo Neutron Noise Simulations,” May 2022, PHYSOR 2022.

H. Belanger, C. Larmier, D. Mancusi, and A. Zoia, “Optimization of Particle
Tracking Methods for Stochastic Media,” May 2022, PHYSOR 2022.

H. Belanger, D. Mancusi, and A. Zoia, “Exact weight cancellation in Monte Carlo
eigenvalue transport problems,” Phys. Rev. E, vol. 104, no. 1, p. 015306, 2021,
doi: 10.1103/physreve.104.015306. 

H. Belanger, D. Mancusi, and A. Zoia, “Solving Eigenvalue Transport Problems
with Negative Weights and Regional Cancellation,” Oct. 2021, M&C 2021,
doi: 10.13182/m&c21-33615. 

## Install
To build MGMC, a linux system with a C++17 compliant compiler is required
(gcc >= 7 or clang >= 6 works), along with cmake >= 3.11. A few third-party
compile-time dependencies ([yaml-cpp](https://github.com/jbeder/yaml-cpp),
[docopt](http://docopt.org/), [pcg](https://www.pcg-random.org),
[ndarray](https://github.com/HunterBelanger/ndarray)) are downloaded
and compiled automatically by CMake during the build.

To build the mgmc executable, the following commands can be used:
```
$ git clone https://github.com/HunterBelanger/mgmc.git
$ cd mgmc 
$ cmake -E make_directory build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
$ cmake --build .
```
This will produce an executable called `mgmc`. You can run the provided example
with the following command:
```
$ ./mgmc -i ref_sqr_c5g7.yaml
```

## Contributors
MGMC was developed by Hunter Belanger in the framework of his Ph.D. thesis
at the French Alternative Energies and Atomic Energy Commission (CEA).
