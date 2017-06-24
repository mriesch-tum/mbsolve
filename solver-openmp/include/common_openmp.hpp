#include <stdlib.h>

#ifndef MBSOLVE_SOLVER_OPENMP_COMMON
#define MBSOLVE_SOLVER_OPENMP_COMMON

#define ALIGN 64

#define __mb_phi_create alloc_if(1) free_if(0)
#define __mb_phi_use alloc_if(0) free_if(0)
#define __mb_phi_delete alloc_if(0) free_if(1)

#ifdef XEON_PHI_OFFLOAD

__attribute__((target(mic))) inline void *mb_aligned_alloc(unsigned int size)
{
    return _mm_malloc(size, ALIGN);
}

__attribute__((target(mic))) inline void mb_aligned_free(void *ptr)
{
    _mm_free(ptr);
}

#define __mb_assume_aligned(ptr) __assume_aligned((ptr), ALIGN)

#else

inline void *mb_aligned_alloc(unsigned int size)
{
    return _mm_malloc(size, ALIGN);
}

inline void mb_aligned_free(void *ptr)
{
    _mm_free(ptr);
}

#define __mb_assume_aligned(ptr) __assume_aligned((ptr), ALIGN)

#endif

#endif
