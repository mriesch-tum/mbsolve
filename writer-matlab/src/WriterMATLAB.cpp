#include <WriterMATLAB.hpp>

#include <mat.h>

namespace mbsolve {

void WriterMATLAB::write(const std::string& file,
			 const std::vector<Result *>& results,
			 const Device& device, const Scenario& scenario)
{
    MATFile *pmat;

    pmat = matOpen(file.c_str(), "w");
    if (pmat == NULL) {
	throw std::exception();
    }

    mxArray *t = mxCreateDoubleScalar(scenario.t_e);

    //matPutVariable(pmat, varname, var);
    matPutVariable(pmat, "t_e", t);

    mxDestroyArray(t);
    matClose(pmat);

}

}
