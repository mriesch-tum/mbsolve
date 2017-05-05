cmake_minimum_required(VERSION 3.6)
project(mbsolve-lib)

include_directories(include/)

add_library(mbsolve-lib SHARED
  #  src/writer.cpp
  src/device.cpp
  src/material.cpp
  src/solver.cpp
  src/Quantity.cpp)
set_property(TARGET mbsolve-lib PROPERTY CXX_STANDARD 11)

target_include_directories(mbsolve-lib PUBLIC include/)

if(MBSOLVE_PYTHON)
  include_directories(${PYTHON_INCLUDE_PATH})
  set_source_files_properties(swig/pymbsolvelib.i PROPERTIES CPLUSPLUS ON)
  # set_source_files_properties(swig/pymbsolvelib.i PROPERTIES SWIG_FLAGS
  #   "-includeall")
  swig_add_module(pymbsolvelib python swig/pymbsolvelib.i)
  set_property(TARGET ${SWIG_MODULE_pymbsolvelib_REAL_NAME} PROPERTY
    CXX_STANDARD 11)
  swig_link_libraries(pymbsolvelib ${PYTHON_LIBRARIES} mbsolve-lib)
endif(MBSOLVE_PYTHON)
