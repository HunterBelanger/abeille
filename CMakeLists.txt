cmake_minimum_required(VERSION 3.11)
project(abeille CXX)

option(ABEILLE_USE_OMP "Compile Abeille with OpenMP for shared memory parallelism" ON)
option(ABEILLE_USE_MPI "Compile Abeille with MPI for distributed memory parallelism" OFF)
option(ABEILLE_GUI_PLOT "Compile Abeille with support for the GUI plotter" ON)

# Get FetchContent for downloading dependencies
include(FetchContent)

#===============================================================================
# Get YAML-CPP version 0.7.0
message(STATUS "Downloading yaml-cpp 0.7.0")
FetchContent_Declare(yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG        yaml-cpp-0.7.0
)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Enable testing")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Enable parse tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Enable contrib stuff in library")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Enable generation of install target")
FetchContent_MakeAvailable(yaml-cpp)

#===============================================================================
# Get PapillonNDL version 0.3.0
message(STATUS "Downloading PapillonNDL develop")
FetchContent_Declare(PapillonNDL
  GIT_REPOSITORY https://github.com/HunterBelanger/papillon-ndl.git
  GIT_TAG        develop
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
target_compile_options(abeille PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
target_link_libraries(abeille PUBLIC docopt_s yaml-cpp PapillonNDL::PapillonNDL NDArray::NDArray PCG::PCG_CXX sobol)
target_compile_definitions(abeille PUBLIC ABEILLE_GIT_HASH=\"${ABEILLE_GIT_HASH_RAW}\")

# Find OpenMP if desired
if(ABEILLE_USE_OMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(abeille PUBLIC OpenMP::OpenMP_CXX)
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
# Copy input files to build directory so that the user can try them right
# after the build.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input_files/c5g7.yaml"
               "${CMAKE_CURRENT_BINARY_DIR}/c5g7.yaml")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input_files/ref_sqr_c5g7.yaml"
               "${CMAKE_CURRENT_BINARY_DIR}/ref_sqr_c5g7.yaml")

# Allow for install of executable
install(TARGETS abeille
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
