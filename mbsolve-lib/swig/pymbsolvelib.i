%module pymbsolvelib
%{
#include "../include/mbsolve.hpp"
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_shared_ptr.i"

%include "../include/types.hpp"
%include "../include/material.hpp"
%include "../include/device.hpp"
