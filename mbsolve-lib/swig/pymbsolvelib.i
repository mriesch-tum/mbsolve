%module pymbsolvelib
%{
#include "../include/mbsolve.hpp"
%}

%include "exception.i"
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
%shared_ptr(mbsolve::record)
%shared_ptr(mbsolve::scenario)
%shared_ptr(mbsolve::sech_pulse)
%shared_ptr(mbsolve::sine_source)
%shared_ptr(mbsolve::source)
%template (vector_double) std::vector<double>;

%include "../include/types.hpp"
%include "../include/material.hpp"
%include "../include/device.hpp"
%include "../include/record.hpp"
%include "../include/result.hpp"
%include "../include/source.hpp"
%include "../include/scenario.hpp"
%include "../include/solver.hpp"
