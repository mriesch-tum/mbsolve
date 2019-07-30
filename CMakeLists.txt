cmake_minimum_required(VERSION 3.6)
project(mbsolve)

set(MBSOLVE_VERSION_MAJOR 0)
set(MBSOLVE_VERSION_MINOR 2)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

# CXX_STANDARD is required
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# specify architecture
set(ARCH "HOST" CACHE STRING "specify architecture")
set_property(CACHE ARCH PROPERTY STRINGS "HOST" "AVX2" "SSE3")

# set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was "
    "specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build."
    FORCE)
  # set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# specify components
option(WITH_HDF5 "with(out) HDF5 writer" ON)
option(WITH_CPU "with(out) CPU solver" ON)
option(WITH_CUDA "with(out) CUDA solver" OFF)
option(WITH_DOXYGEN "with(out) Doxygen documentation" ON)
option(WITH_PYTHON "with(out) Python interface" ON)
option(WITH_CPPTOOL "with(out) C++ tool" ON)

# Eigen library
find_package(Eigen3 REQUIRED)

# Python interface
if(WITH_PYTHON)
  find_package(SWIG 2.0)
  find_package(PythonLibs)
  if(SWIG_FOUND AND PythonLibs_FOUND)
    include(${SWIG_USE_FILE})
    set(MBSOLVE_PYTHON TRUE)
  endif(SWIG_FOUND AND PythonLibs_FOUND)
endif(WITH_PYTHON)

# mbsolve library
add_subdirectory(mbsolve-lib)

# HDF support?
if(WITH_HDF5)
  find_package(HDF5 COMPONENTS C CXX)
  if(HDF5_FOUND)
    set(MBSOLVE_HDF5 TRUE)
    add_subdirectory(writer-hdf5)
  endif(HDF5_FOUND)
endif(WITH_HDF5)

# CUDA support?
if(WITH_CUDA)
  find_package(CUDA 7.0)
  if(CUDA_FOUND)
    set(MBSOLVE_CUDA TRUE)
    add_subdirectory(solver-cuda)
  endif(CUDA_FOUND)
endif(WITH_CUDA)

# CPU support?
if(WITH_CPU)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(MBSOLVE_CPU TRUE)
    add_subdirectory(solver-cpu)
  endif(OPENMP_FOUND)
endif(WITH_CPU)

# C++ tool
if(WITH_CPPTOOL)
  find_package(Boost COMPONENTS program_options system timer)
  if(Boost_FOUND)
    add_subdirectory(mbsolve-tool)
  endif(Boost_FOUND)
endif(WITH_CPPTOOL)

# Doxygen documentation
if(WITH_DOXYGEN)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM)
  endif(DOXYGEN_FOUND)
endif(WITH_DOXYGEN)
