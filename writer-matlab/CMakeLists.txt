cmake_minimum_required(VERSION 2.6)
project(writer-matlab)

find_package(MATLAB REQUIRED)

include_directories(${MBSOLVE_INCLUDE_DIR}
  ${WRITERMATLAB_INCLUDE_DIR}
  ${MATLAB_INCLUDE_DIR})

add_library(writer-matlab SHARED src/WriterMATLAB.cpp)
target_link_libraries(writer-matlab ${MATLAB_LIBRARIES})
