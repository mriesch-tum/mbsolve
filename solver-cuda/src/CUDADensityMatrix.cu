#include <CUDACommon.hpp>
#include <CUDADensityMatrix.hpp>

namespace mbsolve {

__host__
CUDADensityMatrixData::CUDADensityMatrixData(unsigned int numGridPoints,
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
}

/*
__host__ __device__
CUDADensityMatrixData::CUDADensityMatrixData(const CUDADensityMatrixData&
					     other) :
    m_numLevels(other.m_numLevels), m_numMultistep(other.m_numMultistep),
    m_sizeData(other.m_sizeData), m_dmA(other.m_dmA), m_dmB(other.m_dmB),
    m_rhs(other.m_rhs), m_head(other.m_head), m_aIsOld(other.m_aIsOld)
{
}*/


__host__
CUDADensityMatrixData::~CUDADensityMatrixData()
{
}

__host__ unsigned int
CUDADensityMatrixData::getSize() const
{
    return m_sizeData;
}

__host__ __device__ unsigned int
CUDADensityMatrixData::getNumLevels() const
{
    return m_numLevels;
}

__host__ __device__ unsigned int
CUDADensityMatrixData::getNumMultistep() const
{
    return m_numMultistep;
}

__host__ __device__ real *
CUDADensityMatrixData::oldDM(unsigned int row, unsigned int col) const
{
    unsigned int i = (row * m_numLevels + col) * m_numGridPoints;
    return (m_aIsOld ? m_dmA : m_dmB) + i;
}

__host__ __device__ real *
CUDADensityMatrixData::newDM(unsigned int row, unsigned int col) const
{
    unsigned int i = (row * m_numLevels + col) * m_numGridPoints;
    return (m_aIsOld ? m_dmB : m_dmA) + i;
}

__host__ __device__ real *
CUDADensityMatrixData::rhs(unsigned int row, unsigned int col,
			   unsigned int rhsIdx) const
{
    unsigned int idx = (rhsIdx + m_head) % m_numMultistep;
    unsigned int base = (row * m_numLevels + col) * m_numMultistep;
    return m_rhs + (base + idx) * m_numGridPoints;
}

__host__ __device__ void
CUDADensityMatrixData::next()
{
    m_aIsOld = !m_aIsOld;
    m_head = (m_head + 1) % m_numMultistep;
}



CUDADensityMatrix::CUDADensityMatrix(unsigned int numGridPoints,
				     unsigned int numLevels,
				     unsigned int numMultistep) :
    m_data(numGridPoints, numLevels, numMultistep)
    // m_numLevels(numLevels),
    //m_numMultistep(numMultistep)
{
    /* allocate GPU memory */
    chk_err(cudaMalloc(&m_gpuBuffer, m_data.getSize()));

    m_data = CUDADensityMatrixData(numGridPoints, numLevels, numMultistep,
				   m_gpuBuffer);
}

CUDADensityMatrix::~CUDADensityMatrix()
{
    cudaFree(m_gpuBuffer);
}

CUDADensityMatrixData&
CUDADensityMatrix::getData()
{
    return m_data;
}
/*
const CUDADensityMatrixData&
CUDADensityMatrix::getHostData() const
{
    return m_data;
}

void
CUDADensityMatrix::next()
{
    m_data.next();
    }*/

}
