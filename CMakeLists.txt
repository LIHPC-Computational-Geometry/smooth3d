cmake_minimum_required(VERSION 3.14.0)

include(CMakePrintHelpers)
include (cmake/version.cmake)

project(smooth3D 
  LANGUAGES C CXX
  VERSION ${SMOOTH3D_VERSION})

cmake_print_variables(smooth3D_VERSION)

#==============================================================================


option(BUILD_SHARED_LIBS    ON)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/smooth3D)

