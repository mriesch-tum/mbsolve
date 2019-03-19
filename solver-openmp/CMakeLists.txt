cmake_minimum_required(VERSION 3.6)
project(solver-openmp)

# specify matrix exponential method
set(MAT_EXP_METHOD "RODR" CACHE STRING "specify matrix exponential method")
set_property(CACHE MAT_EXP_METHOD PROPERTY STRINGS "PADE" "DIAG" "RODR")

string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(${EIGEN3_INCLUDE_DIR})

# solvers that support Xeon Phi offloading
set(SOLVER_OPENMP_SRCS
  src/solver_openmp_2lvl_pc.cpp
  src/solver_openmp_2lvl_pc_red.cpp
)

# solvers without Xeon Phi support
if(NOT ARCH STREQUAL "PHI")
  list(APPEND SOLVER_OPENMP_SRCS
    src/solver_openmp_2lvl_os.cpp
    src/solver_openmp_clvl_rk.cpp
    src/solver_openmp_2lvl_os_red.cpp
    src/solver_openmp_clvl_os_red.cpp
    src/solver_openmp_fdtd.cpp
  )
endif()

add_library(solver-openmp SHARED ${SOLVER_OPENMP_SRCS})

set_property(TARGET solver-openmp PROPERTY CXX_STANDARD 11)

if(ARCH STREQUAL "AVX2")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-openmp PRIVATE "-mavx2")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-openmp PRIVATE "-xAVX2")
  else()
  endif()
elseif(ARCH STREQUAL "SSE3")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-openmp PRIVATE "-msse3")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-openmp PRIVATE "-xSSE3")
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

if(MAT_EXP_METHOD STREQUAL "DIAG")
  target_compile_definitions(solver-openmp PRIVATE EXP_METHOD=1)
elseif(MAT_EXP_METHOD STREQUAL "RODR")
  target_compile_definitions(solver-openmp PRIVATE EXP_METHOD=2)
else()
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
