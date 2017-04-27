/*
 * mbsolve: Framework for solving the Maxwell-Bloch/-Lioville equations
 *
 * Copyright (c) 2016, Computational Photonics Group, Technical University of
 * Munich.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

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
