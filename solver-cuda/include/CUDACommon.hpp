#ifndef CUDACOMMON_H
#define CUDACOMMON_H

#include <stdexcept>
#include <string>
#include <cuda.h>

namespace mbsolve {

static inline void chk_err(cudaError_t code)
{
    if (code != cudaSuccess) {
	throw std::runtime_error(std::string("CUDA: ") +
				 cudaGetErrorString(code));
    }
}

}

#endif
