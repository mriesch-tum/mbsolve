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
 * \defgroup MBSOLVE_SOLVER_CPU solver-cpu
 * Different solvers that use OpenMP for parallelization.
 */

#include <stdexcept>
#include <stdlib.h>

#ifndef MBSOLVE_SOLVER_CPU_COMMON
#define MBSOLVE_SOLVER_CPU_COMMON

#define ALIGN 64

#ifdef __INTEL_COMPILER
__mb_on_device inline void *mb_aligned_alloc(size_t size)
{
    return _mm_malloc(size, ALIGN);
}

__mb_on_device inline void mb_aligned_free(void *ptr)
{
    _mm_free(ptr);
}

#define __mb_assume_aligned(ptr) __assume_aligned((ptr), ALIGN)

#else

inline void *mb_aligned_alloc(size_t size)
{
    void *addr;
    int ret;

    ret = posix_memalign(&addr, ALIGN, size);

    if (ret != 0) {
        throw std::invalid_argument("posix_memalign failed.");
    }

    return addr;
}

inline void mb_aligned_free(void *ptr)
{
    free(ptr);
}

#define __mb_assume_aligned(ptr) __builtin_assume_aligned(ptr, ALIGN)

#endif

#endif
