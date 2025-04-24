# ==============================================================================
# CMake configuration for cw_coal library and tools
# ==============================================================================

cmake_minimum_required(VERSION 3.10)
project(cw_coal LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Standard install dirs & RPATH (Linux/macOS support)
# ------------------------------------------------------------------------------
include(GNUInstallDirs)

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
  set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_BUILD_RPATH   "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
  set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_BUILD_RPATH   "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# ------------------------------------------------------------------------------
# Default install prefix: if not overridden, install into parent of cw_coal/
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
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist EG)

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
# Include paths
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

ROOT_GENERATE_DICTIONARY(
  G__cw_coalDict
  ${DICT_HEADERS}
  LINKDEF ${LINKDEF_FILE}
  OPTIONS "-I${PROJECT_INCLUDE_DIR}"
)

# ------------------------------------------------------------------------------
# Install dictionary files so ROOT can load them
# ------------------------------------------------------------------------------
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}cw_coalDict_rdict.pcm"
  "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}cw_coalDict.rootmap"
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

# Install the cw_coal shared library into lib/
install(TARGETS cw_coal
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}   # .dylib on macOS
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}   # .so on Linux
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}   # .a if generated
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
# analysis CLI tool
# ------------------------------------------------------------------------------
add_executable(analysis src/app/analysis.cpp)
target_link_libraries(analysis PRIVATE cw_coal ${ROOT_LIBRARIES})
target_compile_definitions(analysis PRIVATE
  DATA_INSTALL_DIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
)
set_target_properties(analysis PROPERTIES
  BUILD_RPATH "${CMAKE_BUILD_RPATH}"
  INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
)
install(TARGETS analysis RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ------------------------------------------------------------------------------
# Tests (optional small binaries)
# ------------------------------------------------------------------------------
function(add_cw_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE cw_coal ${ROOT_LIBRARIES})
  target_compile_definitions(${name}
    PRIVATE DATA_INSTALL_DIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
  )

  # Ensure test uses the *same* install libdir (lib64 or lib) as the main app
  if(APPLE)
    set(testrpath "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  else()
    set(testrpath "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  endif()

  set_target_properties(${name} PROPERTIES
    BUILD_RPATH   "${testrpath}"
    INSTALL_RPATH "${testrpath}"
  )

  # install tests to <prefix>/test alongside bin/ and lib/
  install(TARGETS ${name} RUNTIME DESTINATION test)
endfunction()

add_cw_test(test_combiner        test/test_combiner.cpp)
add_cw_test(test_eventReaderAMPT test/test_eventReaderAMPT.cpp)
add_cw_test(test_pidInference    test/test_pidInference.cpp)

# ------------------------------------------------------------------------------
# Install headers, reference data, scripts
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

# ------------------------------------------------------------------------------
# Install condor directory contents into bin/condor
# ------------------------------------------------------------------------------
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/condor/
  DESTINATION ${CMAKE_INSTALL_BINDIR}/condor
)