%module solvercpu
%{
#include <solver_cpu_loader.hpp>
%}

%init
%{
    mbsolve::solver_cpu_loader cpu_load;
%}

%import(module="mbsolve.lib") "../../mbsolve-lib/include/mbsolve/lib.hpp"

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
