%module writerhdf5
%{
#include <writer_hdf5.hpp>
%}

%init
%{
    mbsolve::writer_hdf5_loader hdf5_load;
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
