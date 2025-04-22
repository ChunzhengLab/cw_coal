# ==============================================================================
# CMake configuration for cw_coal library and tools
# ==============================================================================

cmake_minimum_required(VERSION 3.10)
project(cw_coal LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Standard install dirs & RPATH (macOS support)
# ------------------------------------------------------------------------------
include(GNUInstallDirs)
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(CMAKE_BUILD_RPATH   "@loader_path/../lib")
set(CMAKE_INSTALL_RPATH "@loader_path/../lib")

# ------------------------------------------------------------------------------
# Default install prefix: if not specified, use parent directory of cw_coal/
# ------------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/.."
      CACHE PATH "Install prefix (override with -DCMAKE_INSTALL_PREFIX)" FORCE)
endif()

# ------------------------------------------------------------------------------
# Compiler settings
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Find dependencies
# ------------------------------------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

# ------------------------------------------------------------------------------
# Find nanoflann (mandatory)
# ------------------------------------------------------------------------------
find_path(NANOFLANN_INCLUDE_DIR nanoflann.hpp)

if(NANOFLANN_INCLUDE_DIR)
  message(STATUS "Found nanoflann in ${NANOFLANN_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "nanoflann is required but was not found. Please install nanoflann first.")
endif()

# ------------------------------------------------------------------------------
# Include paths (global)
# ------------------------------------------------------------------------------
set(PROJECT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
include_directories(
  ${PROJECT_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${NANOFLANN_INCLUDE_DIR}
)

# ------------------------------------------------------------------------------
# ROOT dictionary generation
# ------------------------------------------------------------------------------
set(DICT_HEADERS
  ${PROJECT_INCLUDE_DIR}/core/Particle.h
  ${PROJECT_INCLUDE_DIR}/core/Event.h
  ${PROJECT_INCLUDE_DIR}/Combiners.h
)
set(LINKDEF_FILE ${PROJECT_INCLUDE_DIR}/LinkDef.h)
include(${ROOT_DIR}/ROOTMacros.cmake)
ROOT_GENERATE_DICTIONARY(
  G__cw_coalDict
  ${DICT_HEADERS}
  LINKDEF ${LINKDEF_FILE}
  MODULE cw_coal
)
# Install PCM for ROOT introspection
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/libcw_coal_rdict.pcm"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# ------------------------------------------------------------------------------
# Source collections
# ------------------------------------------------------------------------------
file(GLOB_RECURSE SRC_CORE     src/core/*.cxx)
file(GLOB_RECURSE SRC_IO       src/io/*.cxx)
file(GLOB_RECURSE SRC_COMBINER src/combiner/*.cxx)
file(GLOB_RECURSE SRC_ANA      src/ana/*.cxx)
set(SOURCES
  ${SRC_CORE}
  ${SRC_IO}
  ${SRC_COMBINER}
  ${SRC_ANA}
  ${CMAKE_CURRENT_BINARY_DIR}/G__cw_coalDict.cxx
)

# ------------------------------------------------------------------------------
# Build shared library: cw_coal
# ------------------------------------------------------------------------------
add_library(cw_coal SHARED ${SOURCES})
target_compile_definitions(cw_coal
  PUBLIC DATA_INSTALL_DIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
)
target_include_directories(cw_coal PUBLIC
  ${PROJECT_INCLUDE_DIR}
  ${PROJECT_INCLUDE_DIR}/core
  ${PROJECT_INCLUDE_DIR}/io
  ${PROJECT_INCLUDE_DIR}/combiner
)
target_link_libraries(cw_coal PRIVATE ${ROOT_LIBRARIES})
# Install the cw_coal shared library
install(TARGETS cw_coal
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}   # install .dylib on macOS
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}   # install .so on Linux
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}   # install static .a if any
)

# ------------------------------------------------------------------------------
# cwcoal CLI tool
# ------------------------------------------------------------------------------
add_executable(cwcoal src/app/cwcoal.cpp)
target_link_libraries(cwcoal PRIVATE cw_coal ${ROOT_LIBRARIES})
target_compile_definitions(cwcoal
  PRIVATE DATA_INSTALL_DIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
)
install(TARGETS cwcoal RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ------------------------------------------------------------------------------
# Tests: utilities for cw_coal
# ------------------------------------------------------------------------------
function(add_cw_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE cw_coal ${ROOT_LIBRARIES})
  target_compile_definitions(${name}
    PRIVATE
      DATA_INSTALL_DIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
  )
  set_target_properties(${name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test
    BUILD_RPATH "@loader_path/../lib"
    INSTALL_RPATH "@loader_path/../lib"
  )
  install(TARGETS ${name} RUNTIME DESTINATION test)
endfunction()
add_cw_test(test_combiner        test/test_combiner.cpp)
add_cw_test(test_eventReaderAMPT test/test_eventReaderAMPT.cpp)

# ------------------------------------------------------------------------------
# Installation of headers, reference data, scripts, and analysis macros
# ------------------------------------------------------------------------------
install(
  DIRECTORY ${PROJECT_INCLUDE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/refdata/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cw_coal/refdata
  FILES_MATCHING PATTERN "*.root"
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cw_coal/scripts
  FILES_MATCHING PATTERN "*.C" PATTERN "*.py"
)
