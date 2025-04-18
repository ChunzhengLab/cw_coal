cmake_minimum_required(VERSION 3.5)
project(cw_coal)

include_directories(${PROJECT_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile commands support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

# Modern include path setup
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# Include ROOT and local headers
include_directories(${ROOT_INCLUDE_DIRS}
                    ${PROJECT_INCLUDE_DIR}
                    ${PROJECT_INCLUDE_DIR}/core
                    ${PROJECT_INCLUDE_DIR}/io
                    ${PROJECT_INCLUDE_DIR}/combiner)
execute_process(COMMAND brew --prefix nanoflann OUTPUT_VARIABLE NANOFANN_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(${NANOFANN_PREFIX}/include)

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

# Source files
file(GLOB_RECURSE SRC_CORE src/core/*.cxx)
file(GLOB_RECURSE SRC_IO src/io/*.cxx)
file(GLOB_RECURSE SRC_COMBINER src/combiner/*.cxx)
set(SOURCES ${SRC_CORE} ${SRC_IO} ${SRC_COMBINER} ${CMAKE_CURRENT_BINARY_DIR}/G__cw_coalDict.cxx)

# Create shared library
add_library(cw_coal SHARED ${SOURCES})
target_include_directories(cw_coal PUBLIC
  ${PROJECT_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${PROJECT_INCLUDE_DIR}/core
  ${PROJECT_INCLUDE_DIR}/io
  ${PROJECT_INCLUDE_DIR}/combiner
)
target_link_libraries(cw_coal ${ROOT_LIBRARIES})

# Example analyzer macro (ROOT interpreter usage)
# Nothing to compile for analyze_tree.C

# Test target
add_executable(test_combiner test/test_combiner.cpp)
target_link_libraries(test_combiner cw_coal ${ROOT_LIBRARIES})

# Install support
install(TARGETS cw_coal DESTINATION lib)
install(FILES ${DICT_HEADERS} DESTINATION include)