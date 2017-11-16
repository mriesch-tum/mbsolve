cmake_minimum_required(VERSION 3.6)
project(writer-hdf5)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS} include/)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(writer-hdf5 SHARED src/writer_hdf5.cpp)
set_property(TARGET writer-hdf5 PROPERTY CXX_STANDARD 11)

target_include_directories(writer-hdf5 PUBLIC include/)
target_link_libraries(writer-hdf5 mbsolve-lib)
target_link_libraries(writer-hdf5 ${HDF5_LIBRARIES})

if(MBSOLVE_PYTHON)
  include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
  include_directories(${PYTHON_INCLUDE_PATH})
  set_source_files_properties(swig/pywriterhdf5.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(pywriterhdf5 python swig/pywriterhdf5.i)
  set_property(TARGET ${SWIG_MODULE_pywriterhdf5_REAL_NAME} PROPERTY
    CXX_STANDARD 11)
  swig_link_libraries(pywriterhdf5 ${PYTHON_LIBRARIES} mbsolve-lib
    writer-hdf5)
endif(MBSOLVE_PYTHON)
