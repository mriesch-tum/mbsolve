%module writerhdf5
%{
#include "../include/writer_hdf5.hpp"

using namespace mbsolve;
%}

%import(module="mbsolve.lib") "../../mbsolve-lib/include/mbsolve.hpp"

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

%shared_ptr(mbsolve::device)
%shared_ptr(mbsolve::scenario)

%include "../include/writer_hdf5.hpp"
