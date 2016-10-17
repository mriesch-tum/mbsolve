cmake_minimum_required(VERSION 2.6)
project(writer-matlab)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR} include/)

add_library(writer-matlab SHARED src/WriterMATLAB.cpp)
target_include_directories(writer-matlab PUBLIC include/)
target_link_libraries(writer-matlab mbsolve-lib)
target_link_libraries(writer-matlab ${MATLAB_LIBRARIES})
