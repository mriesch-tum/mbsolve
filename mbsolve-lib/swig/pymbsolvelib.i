%module pymbsolvelib
%{
#include "../include/mbsolve.hpp"
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_shared_ptr.i"

%shared_ptr(mbsolve::record)
%template (vector_double) std::vector<double>;

%include "../include/types.hpp"
%include "../include/material.hpp"
%include "../include/device.hpp"
%include "../include/record.hpp"
%include "../include/result.hpp"
%include "../include/source.hpp"
%include "../include/scenario.hpp"
%include "../include/solver.hpp"
