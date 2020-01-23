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
project(mbsolve-lib)

# library
add_library(mbsolve-lib
  src/device.cpp
  src/material.cpp
  src/record.cpp
  src/result.cpp
  src/scenario.cpp
  src/solver.cpp
  src/writer.cpp)
set_target_properties(mbsolve-lib PROPERTIES
  CXX_STANDARD 11
  VERSION ${VERSION_CORE}
  SOVERSION ${VERSION_MAJOR})
target_include_directories(mbsolve-lib PUBLIC "include")

if(MBSOLVE_PYTHON)
  include_directories(${PYTHON_INCLUDE_DIRS})
  set_source_files_properties(swig/lib.i PROPERTIES CPLUSPLUS ON)
  set(CMAKE_SWIG_OUTDIR ${PYTHON_IFACE_PATH})
  swig_add_module(lib python swig/lib.i)
  set_target_properties(${SWIG_MODULE_lib_REAL_NAME} PROPERTIES
    CXX_STANDARD 11
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_IFACE_PATH})
  set(SWIG_MODULE_lib_REAL_NAME ${SWIG_MODULE_lib_REAL_NAME} PARENT_SCOPE)
  target_link_libraries(${SWIG_MODULE_lib_REAL_NAME} PRIVATE mbsolve-lib)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND APPLE)
    target_link_libraries(${SWIG_MODULE_lib_REAL_NAME} PRIVATE
      "-undefined dynamic_lookup")
  endif()
  if(WIN32)
    target_link_libraries(${SWIG_MODULE_lib_REAL_NAME} PRIVATE
      ${PYTHON_LIBRARIES})
  endif(WIN32)
endif(MBSOLVE_PYTHON)

# install library
install(TARGETS mbsolve-lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# install public header files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
