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

/**
 * \defgroup MBSOLVE_SOLVER_CUDA solver-cuda
 * Different solvers that use NVDIA CUDA for parallelization.
 */

#ifndef MBSOLVE_SOLVER_CUDA_COMMON_H
#define MBSOLVE_SOLVER_CUDA_COMMON_H

#include <stdexcept>
#include <string>
#include <cuda.h>
//#include <internal/common_fdtd_2lvl.hpp>
#include <internal/copy_list_entry.hpp>

namespace mbsolve {

const unsigned int MB_CUDA_MAX_CLE = 32;

const unsigned int MB_CUDA_MAX_MATERIALS = 32;

const unsigned int MB_CUDA_MAX_SOURCES = 32;

/* TODO split into different headers for different level types? */
/*
extern __device__ __constant__
sim_constants_2lvl l_sim_consts[MB_CUDA_MAX_MATERIALS];
*/

extern __device__ __constant__
sim_source l_sim_sources[MB_CUDA_MAX_SOURCES];

/* copy list in constant GPU memory */
extern __device__ __constant__
copy_list_entry_dev l_copy_list[MB_CUDA_MAX_CLE];

extern __global__ void
init_memory(real *d, real *e, real *h, unsigned int *indices);

static inline void chk_err(cudaError_t code)
{
    if (code != cudaSuccess) {
	throw std::runtime_error(std::string("CUDA: ") +
				 cudaGetErrorString(code));
    }
}

}

#endif
