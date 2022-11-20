cmake_minimum_required(VERSION 3.11)

project(sobol
  DESCRIPTION "Library to sample sobol sequences."
  LANGUAGES CXX
)

include(CMakePackageConfigHelpers)


# Create PapillonNDL library target
add_library(sobol STATIC "${CMAKE_CURRENT_SOURCE_DIR}/src/sobol.cpp")

# Add public interface includes directory
target_include_directories(sobol PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Comile options for Windows
  target_compile_options(sobol PRIVATE /W4)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
  target_compile_options(sobol PRIVATE -W -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
  target_compile_options(sobol PRIVATE -W -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
  target_compile_options(sobol PRIVATE -W -Wall -Wextra -Wpedantic)
endif()
