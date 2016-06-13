cmake_minimum_required(VERSION 2.6)
project(mbsolve)

find_package(Boost REQUIRED COMPONENTS program_options)

include_directories(${Boost_INCLUDE_DIR})

add_executable(mbsolve mbsolve.cpp)
target_link_libraries(mbsolve ${Boost_LIBRARIES})
