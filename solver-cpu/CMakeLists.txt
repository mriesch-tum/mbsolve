#
# mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
#
# Copyright (c) 2016, Computational Photonics Group, Technical University of
# Munich.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

cmake_minimum_required(VERSION 3.6)
project(solver-cpu)

string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

add_library(solver-cpu
  src/solver_cpu_fdtd.cpp
  src/solver_cpu_fdtd_red.cpp)
set_target_properties(solver-cpu PROPERTIES
  CXX_STANDARD 11
  VERSION ${VERSION_CORE}
  SOVERSION ${VERSION_MAJOR})

if(ARCH STREQUAL "AVX2")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-cpu PRIVATE "-mavx2")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-cpu PRIVATE "-xAVX2")
  else()
  endif()
elseif(ARCH STREQUAL "SSE3")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-cpu PRIVATE "-msse3")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-cpu PRIVATE "-xSSE3")
  else()
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-cpu PRIVATE "-march=native")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-cpu PRIVATE "-xHost")
  else()
  endif()
endif()

target_compile_options(solver-cpu PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(solver-cpu PRIVATE ${mbsolve-lib_INCLUDE_DIR})
target_include_directories(solver-cpu PUBLIC "include")
target_link_libraries(solver-cpu Eigen3::Eigen)
target_link_libraries(solver-cpu ${OpenMP_CXX_FLAGS})
target_link_libraries(solver-cpu mbsolve-lib)

if(MBSOLVE_PYTHON)
  include_directories("include" ${PYTHON_INCLUDE_DIRS}
    ${mbsolve-lib_INCLUDE_DIR})
  set_source_files_properties(swig/solvercpu.i PROPERTIES CPLUSPLUS ON)
  set(CMAKE_SWIG_OUTDIR ${PYTHON_IFACE_PATH})
  swig_add_module(solvercpu python swig/solvercpu.i)
  set_target_properties(${SWIG_MODULE_solvercpu_REAL_NAME} PROPERTIES
    CXX_STANDARD 11
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_IFACE_PATH})
  set(SWIG_MODULE_solvercpu_REAL_NAME ${SWIG_MODULE_solvercpu_REAL_NAME}
    PARENT_SCOPE)
  target_link_libraries(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
    mbsolve-lib solver-cpu)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND APPLE)
    target_link_libraries(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
      "-undefined dynamic_lookup")
  endif()
  if(WIN32)
    target_link_libraries(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
      ${PYTHON_LIBRARIES})
  endif(WIN32)
endif(MBSOLVE_PYTHON)

# install library
install(TARGETS solver-cpu
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# install public header files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
