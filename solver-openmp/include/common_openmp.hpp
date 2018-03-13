#include <stdexcept>
#include <stdlib.h>

#ifndef MBSOLVE_SOLVER_OPENMP_COMMON
#define MBSOLVE_SOLVER_OPENMP_COMMON

#define ALIGN 64

#define __mb_phi_create alloc_if(1) free_if(0)
#define __mb_phi_use alloc_if(0) free_if(0)
#define __mb_phi_delete alloc_if(0) free_if(1)

#ifdef XEON_PHI_OFFLOAD
#define __mb_on_device __attribute__((target(mic)))
#else
#define __mb_on_device
#endif

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
