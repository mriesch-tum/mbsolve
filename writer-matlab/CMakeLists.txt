cmake_minimum_required(VERSION 3.6)
project(writer-matlab)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR} include/)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(writer-matlab SHARED src/writer_matlab.cpp)
set_property(TARGET writer-matlab PROPERTY CXX_STANDARD 11)

target_include_directories(writer-matlab PUBLIC include/)
target_link_libraries(writer-matlab mbsolve-lib)
target_link_libraries(writer-matlab ${MATLAB_LIBRARIES})

if(MBSOLVE_PYTHON)
  include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
  include_directories(${PYTHON_INCLUDE_PATH})
  set_source_files_properties(swig/pywritermatlab.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(pywritermatlab python swig/pywritermatlab.i)
  set_property(TARGET ${SWIG_MODULE_pywritermatlab_REAL_NAME} PROPERTY
    CXX_STANDARD 11)
  swig_link_libraries(pywritermatlab ${PYTHON_LIBRARIES} mbsolve-lib
    writer-matlab)
endif(MBSOLVE_PYTHON)
