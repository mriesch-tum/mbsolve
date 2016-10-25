#include <WriterMATLAB.hpp>

#include <boost/foreach.hpp>
#include <mat.h>

namespace mbsolve {

static WriterFactory<WriterMATLAB> factory("MATLAB");

void
WriterMATLAB::write(const std::string& file,
		    const std::vector<Result *>& results,
		    const Device& device, const Scenario& scenario) const
{
    MATFile *pmat;

    pmat = matOpen(file.c_str(), "w");
    if (pmat == NULL) {
	throw std::invalid_argument("File \"" + file + "\" not found");
    }

    mxArray *t = mxCreateDoubleScalar(scenario.SimEndTime);
    matPutVariable(pmat, "SimEndTime", t);
    mxDestroyArray(t);

    BOOST_FOREACH(mbsolve::Result *result, results) {
	mxArray *var = mxCreateDoubleMatrix(result->size(), 1, mxREAL);

	for (int i = 0; i < result->size(); i++) {
	    *(mxGetPr(var) + i) = result->at(i);
	}

	matPutVariable(pmat, result->name().c_str(), var);

	mxDestroyArray(var);
    }

    matClose(pmat);
}

std::string
WriterMATLAB::getExtension() const
{
    return std::string("mat");
}

}
