%module pymbsolvelib
%{
#include "../include/mbsolve.hpp"
%}

%include "std_string.i"

%include "../include/types.hpp"
%include "../include/material.hpp"
%include "../include/device.hpp"
