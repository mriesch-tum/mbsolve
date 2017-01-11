#include <DensityMatrixData.hpp>

namespace mbsolve {

DensityMatrixData::DensityMatrixData(unsigned int numGridPoints,
					     unsigned int numLevels,
					     unsigned int numMultistep,
					     real *data) :
    m_numGridPoints(numGridPoints), m_numLevels(numLevels),
    m_numMultistep(numMultistep), m_aIsOld(false), m_head(0),
    m_sizeData(sizeof(real) * numLevels * numLevels * (2 + numMultistep) *
	       numGridPoints),
    m_dmA(data), m_dmB(data + numLevels * numLevels * numGridPoints),
    m_rhs(data + numLevels * numLevels * 2 * numGridPoints)
{
    m_dmA = new real[numGridPoints * numLevels * numLevels];
    m_dmB = new real[numGridPoints * numLevels * numLevels];
    m_rhs = new real[numGridPoints * numLevels * numLevels * numMultistep];
}

DensityMatrixData::~DensityMatrixData()
{
    delete[] m_dmA;
    delete[] m_dmB;
    delete[] m_rhs;
}

unsigned int
DensityMatrixData::getSize() const
{
    return m_sizeData;
}

unsigned int
DensityMatrixData::getNumLevels() const
{
    return m_numLevels;
}

unsigned int
DensityMatrixData::getNumMultistep() const
{
    return m_numMultistep;
}

real *
DensityMatrixData::oldDM(unsigned int row, unsigned int col) const
{
    unsigned int i = (row * m_numLevels + col) * m_numGridPoints;
    return (m_aIsOld ? m_dmA : m_dmB) + i;
}

real *
DensityMatrixData::newDM(unsigned int row, unsigned int col) const
{
    unsigned int i = (row * m_numLevels + col) * m_numGridPoints;
    return (m_aIsOld ? m_dmB : m_dmA) + i;
}

real *
DensityMatrixData::rhs(unsigned int row, unsigned int col,
			   unsigned int rhsIdx) const
{
    unsigned int idx = (rhsIdx + m_head) % m_numMultistep;
    unsigned int base = (row * m_numLevels + col) * m_numMultistep;
    return m_rhs + (base + idx) * m_numGridPoints;
}

void
DensityMatrixData::next()
{
    m_aIsOld = !m_aIsOld;
    m_head = (m_head + 1) % m_numMultistep;
}

}
