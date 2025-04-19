# ==============================================================================
# CMake configuration for cw_coal library and tools
# ==============================================================================

# Declare project name and languages
cmake_minimum_required(VERSION 3.10)
project(cw_coal)

# Include GNUInstallDirs to use standard directory variables:
#   CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_DATADIR
include(GNUInstallDirs)
# Enable macOS RPATH so installed executables retain runtime paths
set(CMAKE_MACOSX_RPATH ON)
# Include link directories (e.g., ROOT library path) in RPATH on install
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
# Use loader_path so executables look in sibling lib/ at both build and install time
set(CMAKE_BUILD_RPATH "@loader_path/../lib")
set(CMAKE_INSTALL_RPATH "@loader_path/../lib")

# If user did not set CMAKE_INSTALL_PREFIX, default install into parent directory (CWCoalProject)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/.." CACHE PATH "Install path prefix" FORCE)
endif()


# Use C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile commands support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Locate ROOT libraries and headers
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

# Modern include path setup
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

execute_process(COMMAND brew --prefix nanoflann OUTPUT_VARIABLE NANOFANN_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(${NANOFANN_PREFIX}/include)

# Generate ROOT dictionary for classes stored in ROOT files
# Dictionary generation
set(DICT_HEADERS
  ${PROJECT_INCLUDE_DIR}/core/Particle.h
  ${PROJECT_INCLUDE_DIR}/core/Event.h
  ${PROJECT_INCLUDE_DIR}/Combiners.h
)
set(LINKDEF_FILE ${PROJECT_INCLUDE_DIR}/LinkDef.h)

include(${ROOT_DIR}/ROOTMacros.cmake)
ROOT_GENERATE_DICTIONARY(G__cw_coalDict
  ${DICT_HEADERS}
  LINKDEF ${LINKDEF_FILE}
  MODULE cw_coal
)

# Install ROOT PCM file for cw_coal dictionary
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/libcw_coal_rdict.pcm"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Collect source files for core, I/O, and combiner modules
# Source files
file(GLOB_RECURSE SRC_CORE src/core/*.cxx)
file(GLOB_RECURSE SRC_IO src/io/*.cxx)
file(GLOB_RECURSE SRC_COMBINER src/combiner/*.cxx)
file(GLOB_RECURSE SRC_ANA src/ana/*.cxx)
set(SOURCES ${SRC_CORE} ${SRC_IO} ${SRC_COMBINER} ${SRC_ANA} ${CMAKE_CURRENT_BINARY_DIR}/G__cw_coalDict.cxx)

# Build shared library 'cw_coal' from collected sources
add_library(cw_coal SHARED ${SOURCES})
# Make DATA_INSTALL_DIR available in code for locating refdata at runtime
target_compile_definitions(cw_coal
  PUBLIC DATA_INSTALL_DIR="${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
)
target_include_directories(cw_coal PUBLIC
  ${PROJECT_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${PROJECT_INCLUDE_DIR}/core
  ${PROJECT_INCLUDE_DIR}/io
  ${PROJECT_INCLUDE_DIR}/combiner
)
target_link_libraries(cw_coal ${ROOT_LIBRARIES})

# Build command-line tool 'cwcoal' linking against cw_coal library
# Add executable for app
add_executable(cwcoal src/app/cwcoal.cpp)
target_link_libraries(cwcoal PRIVATE cw_coal ${ROOT_LIBRARIES})
# Pass DATA_INSTALL_DIR to cwcoal executable for locating refdata
target_compile_definitions(cwcoal
  PRIVATE DATA_INSTALL_DIR="${CMAKE_INSTALL_DATADIR}/cw_coal/refdata"
)
install(TARGETS cwcoal RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

## Build test_combiner and output into build/test directory
add_executable(test_combiner test/test_combiner.cpp)
target_link_libraries(test_combiner PRIVATE cw_coal ${ROOT_LIBRARIES})
# Provide install path macros to test_combiner
target_compile_definitions(test_combiner
  PRIVATE
    DATA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}"
    DATA_INSTALL_DATADIR="${CMAKE_INSTALL_DATADIR}"
)
set_target_properties(test_combiner PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test
  BUILD_RPATH "@loader_path/../lib"
  INSTALL_RPATH "@loader_path/../lib"
)
install(
  TARGETS test_combiner
  RUNTIME DESTINATION test
)

# Build test_eventReaderAMPT and output into build/test directory
add_executable(test_eventReaderAMPT test/test_eventReaderAMPT.cpp)
target_link_libraries(test_eventReaderAMPT PRIVATE cw_coal ${ROOT_LIBRARIES})
# Provide install path macros to test_eventReaderAMPT
target_compile_definitions(test_eventReaderAMPT
  PRIVATE
    DATA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}"
    DATA_INSTALL_DATADIR="${CMAKE_INSTALL_DATADIR}"
)
set_target_properties(test_eventReaderAMPT PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test
  BUILD_RPATH "@loader_path/../lib"
  INSTALL_RPATH "@loader_path/../lib"
)
install(
  TARGETS test_eventReaderAMPT
  RUNTIME DESTINATION test
)

# Install the cw_coal shared library
install(TARGETS cw_coal
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install reference ROOT files (refdata)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/refdata/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cw_coal/refdata
  FILES_MATCHING PATTERN "*.root"
)

# Install all public header files (core, io, combiner, ana)
install(
  DIRECTORY ${PROJECT_INCLUDE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
)

# Install analysis scripts
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cw_coal/scripts
  FILES_MATCHING PATTERN "*.C" PATTERN "*.py"
)

# Install analysis macros for ROOT
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/src/ana/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cw_coal/ana
)
