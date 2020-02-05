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
project(writer-hdf5)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS} include/)

add_library(writer-hdf5 src/writer_hdf5.cpp)
set_target_properties(writer-hdf5 PROPERTIES
  CXX_STANDARD 11
  VERSION ${VERSION_CORE}
  SOVERSION ${VERSION_MAJOR})
target_include_directories(writer-hdf5 PUBLIC "include")
target_link_libraries(writer-hdf5 mbsolve-lib)
target_link_libraries(writer-hdf5 ${HDF5_LIBRARIES})

if(MBSOLVE_PYTHON)
  include_directories("include" ${PYTHON_INCLUDE_DIRS}
    ${mbsolve-lib_INCLUDE_DIR})
  set_source_files_properties(swig/writerhdf5.i PROPERTIES CPLUSPLUS ON)
  set(CMAKE_SWIG_OUTDIR ${PYTHON_IFACE_PATH})
  swig_add_module(writerhdf5 python swig/writerhdf5.i)
  set_target_properties(${SWIG_MODULE_writerhdf5_REAL_NAME} PROPERTIES
    CXX_STANDARD 11
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_IFACE_PATH}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PYTHON_IFACE_PATH}
    LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${PYTHON_IFACE_PATH}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PYTHON_IFACE_PATH}
    LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${PYTHON_IFACE_PATH})
  set(SWIG_MODULE_writerhdf5_REAL_NAME ${SWIG_MODULE_writerhdf5_REAL_NAME}
    PARENT_SCOPE)
  target_link_libraries(${SWIG_MODULE_writerhdf5_REAL_NAME} PRIVATE
    mbsolve-lib writer-hdf5)
  if(APPLE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      target_link_libraries(${SWIG_MODULE_writerhdf5_REAL_NAME} PRIVATE
        "-undefined dynamic_lookup")
    endif()
  endif()
  if(WIN32)
    target_link_libraries(${SWIG_MODULE_writerhdf5_REAL_NAME} PRIVATE
      ${PYTHON_LIBRARIES})
  endif()
endif(MBSOLVE_PYTHON)

# install library
install(TARGETS writer-hdf5
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# install public header files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
