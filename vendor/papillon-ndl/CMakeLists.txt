cmake_minimum_required(VERSION 3.11)

project(PapillonNDL
  VERSION "0.2.1"
  DESCRIPTION "Library to interact with ACE nuclear data files cross sections."
  LANGUAGES CXX
)

include(CMakePackageConfigHelpers)

# List of source files for PapillonNDL
set(PNDL_SOURCE_LIST src/ace.cpp
                     src/energy_grid.cpp
                     src/cross_section.cpp
)

# Create PapillonNDL library target
add_library(PapillonNDL STATIC ${PNDL_SOURCE_LIST})

# Add alias to make more friendly with FetchConent
add_library(PapillonNDL::PapillonNDL ALIAS PapillonNDL)

# Add public interface includes directory
target_include_directories(PapillonNDL PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Require C++17 standard
target_compile_features(PapillonNDL PUBLIC cxx_std_17)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Comile options for Windows
  target_compile_options(PapillonNDL PRIVATE /W4)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
  target_compile_options(PapillonNDL PRIVATE -W -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
  target_compile_options(PapillonNDL PRIVATE -W -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
  target_compile_options(PapillonNDL PRIVATE -W -Wall -Wextra -Wpedantic)
endif()
