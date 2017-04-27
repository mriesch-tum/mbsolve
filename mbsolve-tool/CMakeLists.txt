cmake_minimum_required(VERSION 3.6)
project(mbsolve-tool)

find_package(Boost REQUIRED COMPONENTS program_options system timer)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

add_executable(mbsolve-tool src/mbsolve-tool.cpp)
set_property(TARGET mbsolve-tool PROPERTY CXX_STANDARD 11)

target_link_libraries(mbsolve-tool ${Boost_LIBRARIES} mbsolve-lib
  solver-generic)

if(MBSOLVE_CUDA)
  target_link_libraries(mbsolve-tool solver-cuda)
endif(MBSOLVE_CUDA)

if(MBSOLVE_OPENMP)
  target_link_libraries(mbsolve-tool solver-openmp)
endif(MBSOLVE_OPENMP)

if(MBSOLVE_MATLAB)
  target_link_libraries(mbsolve-tool writer-matlab)
endif(MBSOLVE_MATLAB)
