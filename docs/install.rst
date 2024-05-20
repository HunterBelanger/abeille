############
Installation
############
Here, you will find all instructions for compiling Abeille from source.

*************
Prerequisites
*************

Compiler
========
Abeille requires a C++20 compatible compiler. If you are on Windows 10 or 11,
you should install Visual Studio 2019 or later, and this will give you a
sufficient version of the MSVC compiler.

If you are using a Linux/Unix like operating system, you will need either
GCC >= 11 or Clang >= 15. On Ubuntu 22.04 or later and Fedora 36 or later, you
should have a sufficient version of GCC.

CMake
=====
Abeille (and its dependencies) uses CMake as the build system. It is therefore
necessary for the build process.

On Windows, you can download the installer for the latest version of `CMake here <https://cmake.org/download/>`_.

On a linux operating system, you can likely used one of the following commands
to install CMake:

.. code-block:: bash

  sudo dnf install cmake  # Fedora

.. code-block:: bash

  sudo apt install cmake # Ubuntu

HDF5
====
Abeille writes all problem results to an HDF5 file, so you will need to install
the HDF5 library before you can build Abeille. Follow the instructions for your
particular operating system.

If you are building Abeille on Windows, you can download the installer for the
latest version of `HDF5 here <https://www.hdfgroup.org/downloads/hdf5>`_.

On a linux operating system, you can likely used one of the following commands
to install HDF5:

.. code-block:: bash

  sudo dnf install hdf5-devel  # Fedora

.. code-block:: bash

  sudo apt install libhdf5-dev # Ubuntu

MPI (Optional)
==============

If you want to use Abeille on a cluster or supercomputer with distributed memory
parallelism, then you will need to install an MPI library. This dependency is
not necessary if you will not be compiling Abeille with the ABEILLE_USE_MPI
option.

On Windows, you can use the Microsoft MPI distribution `MS-MPI <https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi#ms-mpi-downloads>`_.
Make sure that you install both the stand alone library **and** the SDK !

On linux operating systems, it is recommended to use OpenMPI. This can be
installed with one of the following commands:

.. code-block:: bash

  sudo dnf install openmpi-devel  # Fedora

.. code-block:: bash

  sudo apt install libopenmpi-dev # Ubuntu

Graphics Libraries (Optional)
=============================

Abeille has a built-in graphical interface for making slice plots of the input
geometry. This is very useful for verifying the geometry and materials for your
problem, and identifying areas where particles could get lost durring tracking.

If you are building Abeille on Windows, there are no additional dependencies
which are required to build the GUI plotting interface. It will therefore be
built by default.

When building Abeille on a Linux operating system, the header files for several
differnt parts of the graphics library will likely need to be installed. The
following commands will install all necessary dependencies on Ubuntu or Fedora
based distributions.

.. code-block:: bash
  
  # Fedora
  sudo dnf install  libX11-devel libXrandr-devel libXinerama-devel libXcursor-devel libXi-devel libGL-devel

.. code-block:: bash

  # Ubuntu
  sudo apt install libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl-dev

*****
Build
*****

First, navigate to the location in your file system where you would like to keep
the Abeille sources. Once there, you can get the source files by downloading
them directly from `GitHub <https://github.com/HunterBelanger/abeille>`_, or
using git, with the following command:

.. code-block:: bash

  git clone https://github.com/HunterBelanger/abeille.git

Once this is done, there should be a folder called `abeille`. Move into this
directory with 

.. code-block:: bash

  cd abeille

You can then create and move into a build directory with

.. code-block:: bash
  
  mkdir build
  cd build

At this point, it is time to configure the build opitons. Below, you will find
a list of the specific Abeille build options.

ABEILLE_USE_OMP
  This option is used to compile Abeille with OpenMP for shared memory
  parallelism. This is turned ON by default.

ABEILLE_USE_LTO
  This option is used to build Abeille with link time optimizations, if
  supported by your compiler. This is turned ON by default.

ABEILLE_USE_MPI
  This option builds Abeille with MPI support for distributed memory
  parallelism. Using MPI requires an MPI library be installed, in addition to
  the header files (see `MPI section <#mpi-optional>`_ above). This is turned
  OFF by default.

ABEILLE_GUI_PLOT
  This option builds the graphical geometry plotter with Abeille. This might
  require additional dependencies, depending on your operating system
  (see `Graphics Libraries section <#graphics-libraries-optional>`_ above). It
  is turned OFF by default on Linux systems and turned ON by default on Windows.

To configure your build, you must run the CMake command, with your desired
configuration options. For example, if you desire to compile Abeille with MPI
support and the GUI plotter, then you should run

.. code-block:: bash

  cmake -DABEILLE_USE_MPI=ON -DABEILLE_GUI_PLOT=ON ..

If this command results in errors, make sure you have installed all necessary
dependencies for your desired build configuration. Once the configuration is
complete, you can compile Abeille with

.. code-block:: bash

  cmake --build .

Compliation Optimizations
=========================

How you set the level of compilation optimizations for the build depends on the
build system being used.

If using classic make files on Linux, then you can specify the build type by
using the CMAKE_BUILD_TYPE option with the initial cmake command. For example,
if you want to make a Release build with MPI, then you would run

.. code-block:: bash

  cmake -DCMAKE_BUILD_TYPE=Release -DABEILLE_USE_MPI=ON ..

The available build types are Debug, RelWithDebInfo, and Release.

On Windows, the build type must be specified at the moment of compilation by
passing the --config option with the build type. For example, making a debug
build on Windows would then require:

.. code-block:: bash

  cmake --build . --config=Debug
