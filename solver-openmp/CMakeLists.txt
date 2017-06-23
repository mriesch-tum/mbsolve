cmake_minimum_required(VERSION 3.6)
project(solver-openmp)

#string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

add_library(solver-openmp SHARED src/solver_openmp_2lvl_pc.cpp
  src/solver_openmp_2lvl_pc_red.cpp
)

set_property(TARGET solver-openmp PROPERTY CXX_STANDARD 11)

if(ARCH STREQUAL "AVX2")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-openmp PRIVATE "-mavx2")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-openmp PRIVATE "-xAVX2")
  else()
  endif()
elseif(ARCH STREQUAL "PHI")
  target_compile_definitions(solver-openmp PRIVATE XEON_PHI_OFFLOAD=1)
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-openmp PRIVATE "-march=native")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-openmp PRIVATE "-xHost")
  else()
  endif()
endif()

target_compile_options(solver-openmp PRIVATE ${OpenMP_CXX_FLAGS})

target_include_directories(solver-openmp PUBLIC include/)
target_link_libraries(solver-openmp ${OpenMP_CXX_FLAGS})
target_link_libraries(solver-openmp mbsolve-lib)

if(MBSOLVE_PYTHON)
  include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
  include_directories(${PYTHON_INCLUDE_PATH})
  set_source_files_properties(swig/pysolveropenmp.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(pysolveropenmp python swig/pysolveropenmp.i)
  set_property(TARGET ${SWIG_MODULE_pysolveropenmp_REAL_NAME} PROPERTY
    CXX_STANDARD 11)
  swig_link_libraries(pysolveropenmp ${PYTHON_LIBRARIES} mbsolve-lib
    solver-openmp)
endif(MBSOLVE_PYTHON)
