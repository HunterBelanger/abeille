cmake_minimum_required(VERSION 3.11)
project(abeille CXX)

# These two are for LTO enforcement when used
cmake_policy(SET CMP0069 NEW) 
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

option(ABEILLE_USE_OMP "Compile Abeille with OpenMP for shared memory parallelism" ON)
option(ABEILLE_USE_LTO "Compile with Link-Time Optimization for all components" ON)
option(ABEILLE_USE_MPI "Compile Abeille with MPI for distributed memory parallelism" OFF)
option(ABEILLE_GUI_PLOT "Compile Abeille with support for the GUI plotter" ON)

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
# Get YAML-CPP version 0.8.0
message(STATUS "Downloading yaml-cpp 0.8.0")
FetchContent_Declare(yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG        0.8.0
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
  GIT_REPOSITORY https://github.com/HunterBelanger/papillon-ndl.git
  GIT_TAG        feature/sab
)
set(PNDL_SHARED OFF CACHE BOOL "Build PapillonNDL as a shared library")
set(PNDL_INSTALL OFF CACHE BOOL "Install the PapillonNDL library and header files")
set(PNDL_PYTHON OFF CACHE BOOL "Enable Python interface to PapillonNDL")
FetchContent_MakeAvailable(PapillonNDL)

#===============================================================================
# Get docopt
message(STATUS "Downloading docopt.cpp v0.6.3")
FetchContent_Declare(docopt.cpp
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git
  GIT_TAG        v0.6.3
)
if(NOT docopt.cpp_POPULATED)
  FetchContent_Populate(docopt.cpp)
  add_subdirectory(${docopt.cpp_SOURCE_DIR} ${docopt.cpp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#===============================================================================
# Get NDArray
message(STATUS "Downloading NDArray")
FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray.git
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

#===============================================================================
# Get ImApp
if(ABEILLE_GUI_PLOT)
  message(STATUS "Downloading ImApp")
  FetchContent_Declare(ImApp
    GIT_REPOSITORY https://github.com/HunterBelanger/ImApp.git
    GIT_TAG        master
  )
  FetchContent_MakeAvailable(ImApp)
endif()

#===============================================================================
# Get PCG-C++
message(STATUS "Downloading PCG-C++")
FetchContent_Declare(PCG_CXX
  GIT_REPOSITORY https://github.com/HunterBelanger/pcg-cpp.git
  GIT_TAG        feature/cmake
)
set(PCG_CPP_INSTALL OFF CACHE BOOL "Install PCG-CPP")
FetchContent_MakeAvailable(PCG_CXX)

#===============================================================================
# Get Boost
message(STATUS "Downloading Boost 1.81.0")
set(BOOST_INCLUDE_LIBRARIES unordered)
set(BOOST_ENABLE_CMAKE ON)
FetchContent_Declare(Boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG        boost-1.81.0
)
FetchContent_MakeAvailable(Boost)

#===============================================================================
# Get HighFive
message(STATUS "Downloading HighFive")
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_STATIC_HDF5 ON)
set(HighFive_FIND_QUIETLY ON)
FetchContent_Declare(HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
  GIT_TAG        master
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
target_link_libraries(abeille PUBLIC docopt_s yaml-cpp PapillonNDL::PapillonNDL NDArray::NDArray PCG::PCG_CXX sobol Boost::unordered HighFive)
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
    target_link_libraries(abeille PUBLIC OpenMP::OpenMP_CXX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Compile options for Windows
      target_compile_options(abeille PRIVATE /openmp:llvm) 
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
