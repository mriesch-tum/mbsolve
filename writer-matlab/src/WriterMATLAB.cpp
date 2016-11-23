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

    mxArray *t;
    /* put scenario data */
    t = mxCreateDoubleScalar(scenario.SimEndTime);
    matPutVariable(pmat, "SimEndTime", t);
    mxDestroyArray(t);
    t = mxCreateDoubleScalar(scenario.TimeStepSize);
    matPutVariable(pmat, "TimeStepSize", t);
    mxDestroyArray(t);
    t = mxCreateDoubleScalar(scenario.GridPointSize);
    matPutVariable(pmat, "GridPointSize", t);
    mxDestroyArray(t);
    /* put device data */
    t = mxCreateDoubleScalar(device.XDim());
    matPutVariable(pmat, "XDim", t);
    mxDestroyArray(t);

    /* TODO: use updated (updated by solver) scenario? */
    BOOST_FOREACH(mbsolve::Result *result, results) {
	/* matlab array is created transposed in order to match order */
	mxArray *var = mxCreateDoubleMatrix(result->cols(), result->rows(),
					    mxREAL);

	std::copy(result->data(), result->data() + result->count(),
		  mxGetPr(var));

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
