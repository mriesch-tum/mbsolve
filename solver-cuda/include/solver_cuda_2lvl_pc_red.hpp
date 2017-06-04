#ifndef SOLVERCUDA2LVLRED_H
#define SOLVERCUDA2LVLRED_H

#include <cuda.h>
#include <Solver.hpp>
#include <CUDADensityMatrix.hpp>

namespace mbsolve {

/* TODO: make general? */
static const unsigned int MaxRegions = 8;

/* TODO: class with constructor(Device, Scenario) on CUDA ? */
struct sim_constants
{
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    real w12;
    real d12;
    real tau1;
    real gamma12;

    unsigned int idx_start;
    unsigned int idx_end;

    real d_x_inv;
    real d_t;

    real dm11_init;
    real dm22_init;
};

/* TODO: move generic helper classes for CUDA into separate file */
/* TODO: make inline where possible */
/* TODO: mark host or device where needed */


/* TODO: complex results support ? */
/* TODO: CopyListEntry base class */
/* subclasses for fields and dm? */
/* specialize functions getDst, getSrc? */
/* OR use generic programming/templates */

class CopyListEntry
{
private:
    real *m_dst;

    unsigned int m_interval;
    unsigned int m_position;
    unsigned int m_count;

    RecordType m_type;

public:
    __host__ __device__ CopyListEntry(real *dst,
				      unsigned int count,
				      unsigned int position,
				      unsigned int interval,
				      RecordType t) :
	m_dst(dst), m_count(count),
	m_position(position), m_interval(interval), m_type(t)
    {
    }

    //__host__ __device__ CopyListEntry(

    __host__ __device__ CopyListEntry() { }

    __host__ __device__ ~CopyListEntry()
    {
    }

    __host__ __device__ real *get_dst(unsigned int idx) const {
	return &m_dst[idx/m_interval * m_count];
    }

    __host__ __device__ unsigned int get_size() const {
	return m_count * sizeof(real);
    }

    __host__ __device__ unsigned int get_position() const {
	return m_position;
    }

    __host__ __device__ unsigned int get_count() const {
	return m_count;
    }

    __host__ __device__ RecordType get_type() const { return m_type; }

    __host__ __device__ bool record(unsigned int idx) const {
	return (idx % m_interval) == 0;
    }
};


class SolverCUDA2lvl_red : public ISolver
{
public:
    SolverCUDA2lvl_red(const Device& device, const Scenario& scenario);

    ~SolverCUDA2lvl_red();

    std::string getName() const;

    void run() const;

private:
    cudaStream_t comp_maxwell;
    cudaStream_t copy;

    real *m_h1;
    real *m_e1;
    real *m_d1;

    real *m_h2;
    real *m_e2;
    real *m_d2;


    real *m_src;

    unsigned int *m_indices;

    std::vector<real *> m_results_gpu;
};

}

#endif
