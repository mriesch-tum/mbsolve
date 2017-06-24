#include <stdlib.h>

#ifndef MBSOLVE_SOLVER_OPENMP_COMMON
#define MBSOLVE_SOLVER_OPENMP_COMMON

#define ALIGN 64

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
