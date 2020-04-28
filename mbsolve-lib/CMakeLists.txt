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
project(mbsolve-lib)

# library
add_library(mbsolve-lib
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/device.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/material.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/qm_description.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/record.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/result.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/scenario.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/source.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/types.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mbsolve/lib/writer.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/device.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/material.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/record.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/scenario.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/writer.cpp)
set_target_properties(mbsolve-lib PROPERTIES
  CXX_STANDARD 11
  VERSION ${VERSION_CORE}
  SOVERSION ${VERSION_MAJOR})
target_include_directories(mbsolve-lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
if(MBSOLVE_TESTS)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(mbsolve-lib PRIVATE -fprofile-arcs
      -ftest-coverage)
    target_link_libraries(mbsolve-lib -lgcov)
  endif()
endif()

if(MBSOLVE_PYTHON)
  set_property(SOURCE swig/lib.i PROPERTY CPLUSPLUS ON)
  set(CMAKE_SWIG_OUTDIR ${PYTHON_IFACE_PATH})
  swig_add_library(lib
    LANGUAGE python
    SOURCES swig/lib.i)
  set_target_properties(${SWIG_MODULE_lib_REAL_NAME} PROPERTIES
    CXX_STANDARD 11
    ARCHIVE_OUTPUT_DIRECTORY $<1:${PYTHON_IFACE_PATH}>
    LIBRARY_OUTPUT_DIRECTORY $<1:${PYTHON_IFACE_PATH}>
    RUNTIME_OUTPUT_DIRECTORY $<1:${PYTHON_IFACE_PATH}>)
  target_include_directories(${SWIG_MODULE_lib_REAL_NAME} PRIVATE
    ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${SWIG_MODULE_lib_REAL_NAME} PRIVATE mbsolve-lib)
  if(APPLE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      target_link_libraries(${SWIG_MODULE_lib_REAL_NAME} PRIVATE
        "-undefined dynamic_lookup")
    endif()
  endif()
  if(WIN32)
    target_link_libraries(${SWIG_MODULE_lib_REAL_NAME} PRIVATE
      ${PYTHON_LIBRARIES})
  endif()
  set(PYTHON_MOD_STRING
    "${PYTHON_MOD_STRING}
            '$<TARGET_FILE_NAME:${SWIG_MODULE_lib_REAL_NAME}>',"
    PARENT_SCOPE)
endif()

# Catch2 unit tests
if(MBSOLVE_TESTS)
  add_executable(mbsolve-lib-test
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests_common_fdtd.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests_device.cpp)
  target_link_libraries(mbsolve-lib-test
    PRIVATE Catch2::Catch2 mbsolve-lib)
  set_target_properties(mbsolve-lib-test PROPERTIES CXX_STANDARD 11)
  catch_discover_tests(mbsolve-lib-test)
endif()

# install library
install(TARGETS mbsolve-lib
  EXPORT ${PROJECT}-lib-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# install CMake targets
install(EXPORT ${PROJECT}-lib-targets
  NAMESPACE ${PROJECT}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT})

# install public header files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
