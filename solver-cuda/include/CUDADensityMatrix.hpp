#ifndef CUDADENSITYMATRIX_H
#define CUDADENSITYMATRIX_H

#include <cuda.h>
#include <Types.hpp>

namespace mbsolve {

/* can be passed by value */
class CUDADensityMatrixData
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
    __host__ __device__ CUDADensityMatrixData(unsigned int numGridPoints,
					      unsigned int numLevels,
					      unsigned int numMultistep,
					      real *data = NULL);

    //    __host__ __device__ CUDADensityMatrixData(const CUDADensityMatrixData&
    //					      other);

    __host__ __device__ ~CUDADensityMatrixData();

    __host__ __device__ unsigned int getSize() const;

    __host__ __device__ unsigned int getNumLevels() const;

    __host__ __device__ unsigned int getNumMultistep() const;

    __host__ __device__ void next();

    __host__ __device__ real *oldDM(unsigned int row, unsigned int col) const;

    __host__ __device__ real *newDM(unsigned int row, unsigned int col) const;

    __host__ __device__ real *rhs(unsigned int row, unsigned int col,
				  unsigned int rhsIdx) const;
};

class CUDADensityMatrix
{
private:
    CUDADensityMatrixData m_data;

    /*   const unsigned int m_numLevels;
	 const unsigned int m_numMultistep;*/
    real *m_gpuBuffer;

public:
    CUDADensityMatrix(unsigned int numGridPoints, unsigned int numLevels,
		      unsigned int numMultistep);

    ~CUDADensityMatrix();

    CUDADensityMatrixData& getData();

    /*   void next();

    real *oldDM(unsigned int row, unsigned int col)
    {
	real *ret = m_data.oldDM(row, col);
	return ret;
	}*/
};

}

#endif
