%module lib
%{
#include "../include/mbsolve.hpp"
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

/* TODO:
 *
 * bad_alloc?
 * vector<real>, vector<complex<real> >
 * real type


 */

%shared_ptr(mbsolve::device)
%shared_ptr(mbsolve::qm_operator)
%shared_ptr(mbsolve::qm_superop)
%shared_ptr(mbsolve::qm_lindblad_relaxation)
%shared_ptr(mbsolve::qm_description)
%shared_ptr(mbsolve::qm_desc_2lvl)
%shared_ptr(mbsolve::qm_desc_3lvl)
%shared_ptr(mbsolve::qm_desc_nlvl)
%shared_ptr(mbsolve::material)
%shared_ptr(mbsolve::record)
%shared_ptr(mbsolve::region)
%shared_ptr(mbsolve::result)
%shared_ptr(mbsolve::scenario)
%shared_ptr(mbsolve::gaussian_pulse)
%shared_ptr(mbsolve::sech_pulse)
%shared_ptr(mbsolve::sine_source)
%shared_ptr(mbsolve::single_cycle_pulse)
%shared_ptr(mbsolve::solver)
%shared_ptr(mbsolve::source)
%shared_ptr(mbsolve::writer)

%include "../include/types.hpp"

%template(vector_real) std::vector<mbsolve::real>;
%template(vector_complex) std::vector<mbsolve::complex>;
%template(vector_result) std::vector<std::shared_ptr<mbsolve::result> >;
%template(matrix_real) std::vector<std::vector<mbsolve::real> >;

%include "../include/qm_description.hpp"
%include "../include/material.hpp"
%include "../include/device.hpp"
%include "../include/record.hpp"
%include "../include/result.hpp"
%include "../include/source.hpp"
%include "../include/scenario.hpp"
%include "../include/solver.hpp"
%include "../include/writer.hpp"
