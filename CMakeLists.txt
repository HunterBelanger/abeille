cmake_minimum_required(VERSION 3.16...4.0)
project(abeille CXX)

# These two are for LTO enforcement when used
cmake_policy(SET CMP0069 NEW) 
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

option(ABEILLE_USE_OMP "Compile Abeille with OpenMP for shared memory parallelism" ON)
option(ABEILLE_USE_LTO "Compile with Link-Time Optimization for all components" ON)
option(ABEILLE_USE_MPI "Compile Abeille with MPI for distributed memory parallelism" OFF)

# If building Abeille on Windows, the default setting is to build with the GUI
# plotter, as no other dependencies are necessary. On Linux, additional libraries
# are required to build the GUI, so we default this to OFF.
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  option(ABEILLE_GUI_PLOT "Compile Abeille with support for the GUI plotter" ON)
else()
  option(ABEILLE_GUI_PLOT "Compile Abeille with support for the GUI plotter" OFF)
endif()

# Check for using LTO
if(ABEILLE_USE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)

  if(LTO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "Using LTO for build")
  else()
    message(WARNING "LTO is not supported: ${LTO_ERROR}")
    message(WARNING "Will continue build WITHOUT LTO")
  endif()
endif()

# Get FetchContent for downloading dependencies
include(FetchContent)

#===============================================================================
# Get YAML-CPP version v0.9.0
message(STATUS "Downloading yaml-cpp v0.9.0")
FetchContent_Declare(yaml-cpp
  URL            https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.9.0.tar.gz
  EXCLUDE_FROM_ALL
)
set(YAML_CPP_BUILD_SHARED_LIBS OFF CACHE BOOL "Build yaml-cpp shared library")
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Enable testing")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Enable parse tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Enable contrib stuff in library")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Enable generation of install target")
FetchContent_MakeAvailable(yaml-cpp)

#===============================================================================
# Get PapillonNDL develop
message(STATUS "Downloading PapillonNDL feature/sab branch")
FetchContent_Declare(PapillonNDL
  URL      https://github.com/HunterBelanger/papillon-ndl/archive/f1e3416ccbb1b8afd436bf379e821fd7b44dd903.tar.gz
  EXCLUDE_FROM_ALL
)
set(PNDL_SHARED OFF CACHE BOOL "Build PapillonNDL as a shared library")
set(PNDL_INSTALL OFF CACHE BOOL "Install the PapillonNDL library and header files")
set(PNDL_PYTHON OFF CACHE BOOL "Enable Python interface to PapillonNDL")
FetchContent_MakeAvailable(PapillonNDL)

#===============================================================================
# Get docopt
message(STATUS "Downloading docopt.cpp v0.6.3")
FetchContent_Declare(docopt.cpp
  URL            https://github.com/docopt/docopt.cpp/archive/refs/tags/v0.6.3.tar.gz
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(docopt.cpp)

#===============================================================================
# Get XTL
message(STATUS "Downloading xtl v0.7.7")
FetchContent_Declare(xtl
  URL      https://github.com/xtensor-stack/xtl/archive/refs/tags/0.7.7.tar.gz
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(xtl)

#===============================================================================
# Get XSIMD
message(STATUS "Downloading xsimd v12.1.1")
FetchContent_Declare(xsimd
  URL      https://github.com/xtensor-stack/xsimd/archive/refs/tags/12.1.1.tar.gz
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(xsimd)

#===============================================================================
# Get XTENSOR
message(STATUS "Downloading xtensor v0.25.0")
FetchContent_Declare(xtensor
  URL https://github.com/xtensor-stack/xtensor/archive/refs/tags/0.25.0.tar.gz
  EXCLUDE_FROM_ALL
)
set(XTENSOR_USE_XSIMD ON CACHE BOOL "Enable SIMD acceleration.")
FetchContent_MakeAvailable(xtensor)

#===============================================================================
# Get ImApp
if(ABEILLE_GUI_PLOT)
  message(STATUS "Downloading ImApp")
  FetchContent_Declare(ImApp
    URL      https://github.com/HunterBelanger/ImApp/archive/7bcee00689165944ca608266ce31ee9ea488ef3f.tar.gz
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(ImApp)
endif()

#===============================================================================
# Get PCG-C++
message(STATUS "Downloading PCG-C++")
FetchContent_Declare(PCG_CXX
  URL            https://github.com/HunterBelanger/pcg-cpp/archive/eb816d2f01862399b67a5de68e3d423c9c73d2e2.tar.gz
  EXCLUDE_FROM_ALL
)
set(PCG_CPP_INSTALL OFF CACHE BOOL "Install PCG-CPP")
FetchContent_MakeAvailable(PCG_CXX)

#===============================================================================
# Get Boost
message(STATUS "Downloading Boost 1.90.0")
set(BOOST_INCLUDE_LIBRARIES unordered container)
set(BOOST_ENABLE_CMAKE ON)
FetchContent_Declare(Boost
  URL            https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-cmake.tar.gz
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(Boost)

#===============================================================================
# Get HDF5
message(STATUS "Downloading HDF5 v2.0.0")
set(HDF5_EXTERNALLY_CONFIGURED  1)
set(HDF5_ENABLE_SZIP_SUPPORT  OFF)
set(HDF5_ENABLE_ZLIB_SUPPORT OFF)
set(HDF5_BUILD_EXAMPLES       OFF)
set(HDF5_BUILD_TOOLS          OFF)
set(HDF5_BUILD_UTILS          OFF)
set(HDF5_BUILD_HL_LIB         OFF)
set(BUILD_SHARED_LIBS         OFF)
FetchContent_Declare(hdf5
  URL      https://github.com/HDFGroup/hdf5/releases/download/2.0.0/hdf5-2.0.0.tar.gz
  URL_HASH MD5=801597eccdf74d7304444fdc18652b75
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(hdf5)

#===============================================================================
# Get HighFive
message(STATUS "Downloading HighFive")
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_TEST_SPAN OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_FIND_HDF5 OFF)
FetchContent_Declare(HighFive
  URL      https://github.com/highfive-devs/highfive/archive/6ea1c7b8b38db5a58acf33f011f442e0c28f0aa7.tar.gz # Dec. 16 2025
  URL_HASH MD5=6fdc23148460812548f67b07019b209c
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(HighFive)

#===============================================================================
# Get SOBOL
add_subdirectory("${CMAKE_SOURCE_DIR}/vendor/sobol")

#===============================================================================
include(sourcelist.cmake) # For ABEILLE_SOURCE_FILES list

# Get the Git hash for define
execute_process(COMMAND git describe --dirty --always
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                OUTPUT_VARIABLE ABEILLE_GIT_HASH_RAW
                OUTPUT_STRIP_TRAILING_WHITESPACE
               )

add_executable(abeille src/main.cpp ${ABEILLE_SOURCE_FILES})
target_include_directories(abeille PRIVATE include)
target_compile_features(abeille PRIVATE cxx_std_17)
target_link_libraries(abeille PUBLIC docopt_s yaml-cpp xtl xtensor PapillonNDL::PapillonNDL PCG::PCG_CXX sobol Boost::unordered Boost::container HighFive hdf5-static)
target_compile_definitions(abeille PUBLIC ABEILLE_GIT_HASH=\"${ABEILLE_GIT_HASH_RAW}\")
target_compile_definitions(abeille PUBLIC ABEILLE_COMPILER_NAME=\"${CMAKE_CXX_COMPILER_ID}\")
target_compile_definitions(abeille PUBLIC ABEILLE_COMPILER_VERSION=\"${CMAKE_CXX_COMPILER_VERSION}\")

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Compile options for Windows
  target_compile_options(abeille PRIVATE /W3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
  target_compile_options(abeille PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
  target_compile_options(abeille PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
  target_compile_options(abeille PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
endif()

# Find OpenMP if desired
if(ABEILLE_USE_OMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Compile options for Windows
      target_compile_options(abeille PRIVATE /openmp:llvm) 
    else()
      target_link_libraries(abeille PUBLIC OpenMP::OpenMP_CXX)
    endif()
    target_compile_definitions(abeille PUBLIC ABEILLE_USE_OMP)
  endif()
endif()

# Find MPI if desired
if(ABEILLE_USE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_CXX_FOUND)
    target_link_libraries(abeille PUBLIC MPI::MPI_CXX)
    target_compile_definitions(abeille PUBLIC ABEILLE_USE_MPI)
  endif()
endif()

if(ABEILLE_GUI_PLOT)
  target_link_libraries(abeille PUBLIC ImApp::ImApp)
  target_compile_definitions(abeille PUBLIC ABEILLE_GUI_PLOT)
endif()

#===============================================================================
# Allow for install of executable
install(TARGETS abeille
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
