.. Abeille documentation master file, created by
   sphinx-quickstart on Sat Sep  9 18:31:49 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Abeille : A Monte Carlo Transport Code
======================================

Abeille is a 3D Monte Carlo transport code, which solves the Boltzmann transport
equation for neutrons. It is capable of solving fixed-source, k-eigenvalue, and
neutron noise problems with both continuous energy and multi-group physics.
ACE files are used for continuous energy data, like other popular transport
codes. Input files are written in easy to read yaml, which is very friendly for
new users. The geometry may be constructed using combinations of surface
half-spaces to create cells. Like in other Monte Carlo codes, Universes and
Lattices may be used to construct more complex geometries such as for fuel
assemblies found in nuclear reactors. An interactive graphical interface is
built in, for plotting the problem geometry for verification. Problem results
are written to and HDF5 file at the end of the simulation.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   install
   nuclear_data

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
