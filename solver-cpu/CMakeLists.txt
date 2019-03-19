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

cmake_minimum_required(VERSION 3.9)
project(solver-cpu)

string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

set(SRC_BASE_CPU ${CMAKE_CURRENT_SOURCE_DIR})
add_library(solver-cpu
  ${SRC_BASE_CPU}/include/mbsolve/solver-cpu/internal/common_cpu.hpp
  ${SRC_BASE_CPU}/include/mbsolve/solver-cpu/solver_cpu_fdtd.hpp
  ${SRC_BASE_CPU}/include/mbsolve/solver-cpu/solver_cpu_fdtd_red.hpp
  ${SRC_BASE_CPU}/include/mbsolve/solver-cpu/solver_cpu_fdtd.tpp
  ${SRC_BASE_CPU}/include/mbsolve/solver-cpu/solver_cpu_fdtd_red.tpp
  ${SRC_BASE_CPU}/include/mbsolve/solver-cpu/solver_cpu_loader.hpp
  ${SRC_BASE_CPU}/src/solver_cpu_fdtd_noop.cpp
  ${SRC_BASE_CPU}/src/solver_cpu_fdtd_red_noop.cpp
  ${SRC_BASE_CPU}/src/solver_cpu_fdtd_cvr_rodr.cpp
  ${SRC_BASE_CPU}/src/solver_cpu_fdtd_red_cvr_rodr.cpp
  ${SRC_BASE_CPU}/src/solver_cpu_fdtd_reg_cayley.cpp
  ${SRC_BASE_CPU}/src/solver_cpu_fdtd_red_reg_cayley.cpp
  ${SRC_BASE_CPU}/src/solver_cpu_loader.cpp)
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

# MSVC does not support OpenMP pragma "simd"
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  target_compile_options(solver-cpu PRIVATE "-DUSE_OMP_SIMD=0")
else()
  target_compile_options(solver-cpu PRIVATE "-DUSE_OMP_SIMD=1")
endif()

target_compile_options(solver-cpu PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(solver-cpu PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(solver-cpu
  PRIVATE Eigen3::Eigen ${OpenMP_CXX_LIBRARIES}
  PUBLIC mbsolve-lib)

if(MBSOLVE_PYTHON)
  set_property(SOURCE swig/solvercpu.i PROPERTY CPLUSPLUS ON)
  set(CMAKE_SWIG_OUTDIR ${PYTHON_IFACE_PATH})
  swig_add_library(solvercpu
    LANGUAGE python
    SOURCES swig/solvercpu.i)
  set_target_properties(${SWIG_MODULE_solvercpu_REAL_NAME} PROPERTIES
    CXX_STANDARD 11
    ARCHIVE_OUTPUT_DIRECTORY $<1:${PYTHON_IFACE_PATH}>
    LIBRARY_OUTPUT_DIRECTORY $<1:${PYTHON_IFACE_PATH}>
    RUNTIME_OUTPUT_DIRECTORY $<1:${PYTHON_IFACE_PATH}>)
  target_include_directories(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
    ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
    mbsolve-lib solver-cpu)
  if(APPLE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      target_link_libraries(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
        "-undefined dynamic_lookup")
    endif()
  endif()
  if(WIN32)
    target_link_libraries(${SWIG_MODULE_solvercpu_REAL_NAME} PRIVATE
      ${PYTHON_LIBRARIES})
  endif()
  set(PYTHON_MOD_STRING
    "${PYTHON_MOD_STRING}
            '$<TARGET_FILE_NAME:${SWIG_MODULE_solvercpu_REAL_NAME}>',"
    PARENT_SCOPE)
endif()

# install library
install(TARGETS solver-cpu
  EXPORT ${PROJECT}-solver-cpu-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# install CMake targets
install(EXPORT ${PROJECT}-solver-cpu-targets
  NAMESPACE ${PROJECT}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT})

# install public header files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
