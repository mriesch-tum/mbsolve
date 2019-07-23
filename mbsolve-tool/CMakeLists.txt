cmake_minimum_required(VERSION 3.6)
project(mbsolve-tool)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(mbsolve-tool src/mbsolve-tool.cpp)
set_property(TARGET mbsolve-tool PROPERTY CXX_STANDARD 11)

target_link_libraries(mbsolve-tool ${Boost_LIBRARIES} mbsolve-lib
  solver-generic)

if(MBSOLVE_CPU)
  target_link_libraries(mbsolve-tool solver-cpu)
endif(MBSOLVE_CPU)

if(MBSOLVE_CUDA)
  target_link_libraries(mbsolve-tool solver-cuda)
endif(MBSOLVE_CUDA)

if(MBSOLVE_MATLAB)
  target_link_libraries(mbsolve-tool writer-matlab)
endif(MBSOLVE_MATLAB)

if(MBSOLVE_HDF5)
  target_link_libraries(mbsolve-tool writer-hdf5)
endif(MBSOLVE_HDF5)
