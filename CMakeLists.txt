cmake_minimum_required(VERSION 3.6)
project(mbsolve)

set(MBSOLVE_VERSION_MAJOR 0)
set(MBSOLVE_VERSION_MINOR 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

# CXX_STANDARD is required
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# specify architecture
set(ARCH "HOST" CACHE STRING "specify architecture")
set_property(CACHE ARCH PROPERTY STRINGS "HOST" "PHI" "AVX2" "SSE3")

# specify components
option(WITH_MATLAB "with(out) MATLAB writer" ON)
option(WITH_CUDA "with(out) CUDA solver" ON)
option(WITH_OPENMP "with(out) OPENMP solver" ON)
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

# MATLAB support?
if(WITH_MATLAB)
  find_package(MATLAB)
  if(MATLAB_FOUND)
    set(MBSOLVE_MATLAB TRUE)
    add_subdirectory(writer-matlab)
  endif(MATLAB_FOUND)
endif(WITH_MATLAB)

# CUDA support?
if(WITH_CUDA)
  find_package(CUDA 7.0)
  if(CUDA_FOUND)
    set(MBSOLVE_CUDA TRUE)
    add_subdirectory(solver-cuda)
  endif(CUDA_FOUND)
endif(WITH_CUDA)

# OpenMP support?
if(WITH_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(MBSOLVE_OPENMP TRUE)
    add_subdirectory(solver-openmp)
  endif(OPENMP_FOUND)
endif(WITH_OPENMP)

# Generic support
add_subdirectory(solver-generic)

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
