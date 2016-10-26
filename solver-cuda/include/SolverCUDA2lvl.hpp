#ifndef SOLVERCUDA2LVL_H
#define SOLVERCUDA2LVL_H

#include <Solver.hpp>

namespace mbsolve {

static const unsigned int NumEntries = 4;
static const unsigned int NumMultistep = 5;
static const unsigned int MaxRegions = 8;

struct sim_constants
{
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    real w12;
    real d12;
    real gamma1;
    real gamma2;

    unsigned int idx_start;
    unsigned int idx_end;

    real d_x;
    real d_t;
};


class DensityMatrix
{
public:
    DensityMatrix();

    ~DensityMatrix();

    void initialize(unsigned int numGridPoints);

    void next();

    __host__ __device__ __inline__ real *OldDM(unsigned int entry) const;

    __host__ __device__ __inline__ real *NewDM(unsigned int entry) const;

    __device__ __inline__ real *RHS(unsigned int entry, unsigned int row) const;

private:
    real *dm_a[NumEntries];
    real *dm_b[NumEntries];
    real *rhs[NumEntries][NumMultistep];

    bool a_is_old;
    unsigned int head;
};


/* functor classes */

class GetSrcField
{
private:
    real *m_address;

public:
    GetSrcField(real *address) :
	m_address(address)
    { }

    real *operator()()
    {
	return m_address;
    }

};

/* TODO: complex results support ? */

class GetSrcDensity
{
private:
    DensityMatrix *m_dm;
    unsigned int m_row;

public:
    GetSrcDensity(DensityMatrix *dm, unsigned int row) :
	m_dm(dm), m_row(row)
    {
    }

    real *operator()()
    {
	return m_dm->OldDM(m_row);
    }

};

/* TODO: CopyListEntry base class */
/* subclasses for fields and dm? */
/* specialize functions getDst, getSrc? */
/* OR use generic programming/templates */

class CopyListEntry
{
private:
    real *m_src;
    Result *m_res;
    unsigned int m_size;
    unsigned int m_interval;

    /* TODO: base address + offset (position */

public:
    CopyListEntry(real *src, Result *result, unsigned int count,
		  unsigned int interval) :
	m_src(src), m_res(result), m_size(sizeof(real) * count),
	m_interval(interval)
    {
    }

    real *getSrc() const { return m_src; }

    real *getDst(unsigned int idx) const {
	return m_res->data(idx / m_interval);
    }

    unsigned int getSize() const { return m_size; }

    bool record(unsigned int idx) const { return (idx % m_interval) == 0; }
};



class SolverCUDA2lvl : public ISolver
{
public:
    SolverCUDA2lvl(const Device& device, const Scenario& scenario);

    ~SolverCUDA2lvl();

    std::string getName() const;

    void run(const std::vector<Result *>& results) const;

private:
    cudaStream_t comp_maxwell;
    cudaStream_t comp_density;
    cudaStream_t copy;

    DensityMatrix dm;

    real *h;
    real *e;

    std::vector<CopyListEntry> m_copyListBlack;
    std::vector<CopyListEntry> m_copyListRed;
    std::vector<Result *> m_results;
};

}

#endif
