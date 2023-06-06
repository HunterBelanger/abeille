# Abeille - A Monte Carlo Transport Code

Abeille is a 3D multi-group Monte Carlo transport code which solves the Boltzmann
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
can be turned on at compiled time with the `-DABEILLE_USE_MPI=ON` opiton
when running cmake.

In addition to solving standard k-eigenvalue and fixed-source problems, Abeille
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

## Install
To build Abeille, a linux system with a C++20 compliant compiler is required
(gcc >= 11 or clang >= 15 works), along with cmake >= 3.11. A few third-party
compile-time dependencies ([papillon-ndl](https://github.com/HunterBelanger/papillon-ndl),
[yaml-cpp](https://github.com/jbeder/yaml-cpp),
[docopt](http://docopt.org/), [pcg](https://www.pcg-random.org),
[ImApp](https://github.com/HunterBelanger/ImApp),
[ndarray](https://github.com/HunterBelanger/ndarray)) are downloaded
and compiled automatically by CMake during the build.

To build the Abeille executable, the following commands can be used:
```
$ git clone https://github.com/HunterBelanger/abeille.git
$ cd abeille
$ cmake -E make_directory build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
$ cmake --build .
```
This will produce an executable called `abeille`. You can run the provided example
with the following command:
```
$ ./abeille -i c5g7.yaml
```

## Contributors
Abeille is based on MGMC, which was developed by Hunter Belanger in the
framework of his Ph.D. thesis at the French Alternative Energies and
Atomic Energy Commission (CEA).
