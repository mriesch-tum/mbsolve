#ifndef DENSITYMATRIX_H
#define DENSITYMATRIX_H

#include <Types.hpp>

namespace mbsolve {

class DensityMatrixData
{
private:
    unsigned int m_numGridPoints;
    unsigned int m_numLevels;
    unsigned int m_numMultistep;
    unsigned int m_sizeData;

    bool m_aIsOld;
    unsigned int m_head;

    real * m_dmA;
    real * m_dmB;
    real * m_rhs;

public:
    DensityMatrixData(unsigned int numGridPoints,
		      unsigned int numLevels,
		      unsigned int numMultistep,
		      real *data = NULL);

    ~DensityMatrixData();

    unsigned int getSize() const;

    unsigned int getNumLevels() const;

    unsigned int getNumMultistep() const;

    void next();

    real *oldDM(unsigned int row, unsigned int col) const;

    real *newDM(unsigned int row, unsigned int col) const;

    real *rhs(unsigned int row, unsigned int col,
	      unsigned int rhsIdx) const;
};

}

#endif
