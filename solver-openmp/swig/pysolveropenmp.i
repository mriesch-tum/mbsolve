%module pysolveropenmp
%{
#include "../include/solver_openmp_2lvl_pc.hpp"

using namespace mbsolve;
%}

%import(module="pymbsolvelib") "../../mbsolve-lib/include/mbsolve.hpp"

%include "exception.i"
%include "std_except.i"
%include "std_map.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_vector.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%shared_ptr(mbsolve::solver_openmp_2lvl_pc)

%include "../include/solver_openmp_2lvl_pc.hpp"
