cmake_minimum_required(VERSION 3.6)
project(writer-matlab)

find_package(Boost REQUIRED COMPONENTS program_options system timer)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR} include/  ${Boost_INCLUDE_DIR})

add_library(writer-matlab SHARED src/WriterMATLAB.cpp)
target_include_directories(writer-matlab PUBLIC include/)
target_link_libraries(writer-matlab mbsolve-lib)
target_link_libraries(writer-matlab ${MATLAB_LIBRARIES})
