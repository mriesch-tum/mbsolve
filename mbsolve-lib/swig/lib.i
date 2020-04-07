%module lib
%{
#include <mbsolve/lib.hpp>
%}

%include "stl.i"
%include "std_complex.i"
%include "std_except.i"
%include "std_map.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_vector.i"

%exception %{
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch(...) {
  }
%}

%shared_ptr(mbsolve::device)
%shared_ptr(mbsolve::gaussian_pulse)
%shared_ptr(mbsolve::ic_density)
%shared_ptr(mbsolve::ic_density_const)
%shared_ptr(mbsolve::ic_field)
%shared_ptr(mbsolve::ic_field_const)
%shared_ptr(mbsolve::ic_field_random)
%shared_ptr(mbsolve::material)
%shared_ptr(mbsolve::qm_operator)
%shared_ptr(mbsolve::qm_superop)
%shared_ptr(mbsolve::qm_lindblad_relaxation)
%shared_ptr(mbsolve::qm_description)
%shared_ptr(mbsolve::qm_desc_2lvl)
%shared_ptr(mbsolve::record)
%shared_ptr(mbsolve::region)
%shared_ptr(mbsolve::result)
%shared_ptr(mbsolve::scenario)
%shared_ptr(mbsolve::sech_pulse)
%shared_ptr(mbsolve::sine_source)
%shared_ptr(mbsolve::single_cycle_pulse)
%shared_ptr(mbsolve::solver)
%shared_ptr(mbsolve::source)
%shared_ptr(mbsolve::writer)

%include "../include/mbsolve/lib/types.hpp"

%template(vector_string) std::vector<std::string>;
%template(vector_real_d) std::vector<double>;
%template(vector_real_f) std::vector<float>;
%template(vector_complex_d) std::vector<std::complex<double> >;
%template(vector_complex_f) std::vector<std::complex<float> >;
%template(vector_result) std::vector<std::shared_ptr<mbsolve::result> >;
%template(matrix_real_d) std::vector<std::vector<double> >;
%template(matrix_real_f) std::vector<std::vector<float> >;

%feature("python:cdefaultargs") mbsolve::material::material;
%feature("python:cdefaultargs") mbsolve::scenario::scenario;

%include "../include/mbsolve/lib/qm_description.hpp"
%include "../include/mbsolve/lib/material.hpp"
%include "../include/mbsolve/lib/device.hpp"
%include "../include/mbsolve/lib/record.hpp"
%include "../include/mbsolve/lib/result.hpp"
%include "../include/mbsolve/lib/source.hpp"
%include "../include/mbsolve/lib/scenario.hpp"
%include "../include/mbsolve/lib/solver.hpp"
%include "../include/mbsolve/lib/writer.hpp"
