cmake_minimum_required(VERSION 3.6)
project(solver-generic)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(solver-generic SHARED src/solver_generic.cpp)
set_property(TARGET solver-generic PROPERTY CXX_STANDARD 11)

target_include_directories(solver-generic PUBLIC include/)
target_link_libraries(solver-generic mbsolve-lib)

if(MBSOLVE_PYTHON)
  include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
  include_directories(${PYTHON_INCLUDE_PATH})
  set_source_files_properties(swig/pysolvergeneric.i PROPERTIES CPLUSPLUS ON)
#  set_source_files_properties(swig/pysolvergeneric.i PROPERTIES SWIG_FLAGS
#    "-includeall")
  swig_add_module(pysolvergeneric python swig/pysolvergeneric.i)
  set_property(TARGET ${SWIG_MODULE_pysolvergeneric_REAL_NAME} PROPERTY
    CXX_STANDARD 11)
  swig_link_libraries(pysolvergeneric ${PYTHON_LIBRARIES} mbsolve-lib
    solver-generic)
endif(MBSOLVE_PYTHON)
