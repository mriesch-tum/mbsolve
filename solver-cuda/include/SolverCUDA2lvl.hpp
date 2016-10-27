#ifndef SOLVERCUDA2LVL_H
#define SOLVERCUDA2LVL_H

#include <Solver.hpp>

namespace mbsolve {

/* TODO: leave static members here? */
/* TODO: namespace separation? */
static const unsigned int NumLevels = 2;
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
    real tau1;
    real gamma12;

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

    __host__ __device__ __inline__ real *OldDM(unsigned int row,
					       unsigned int col) const;

    __host__ __device__ __inline__ real *NewDM(unsigned int row,
					       unsigned int col) const;

    __device__ __inline__ real *RHS(unsigned int row, unsigned int col,
				    unsigned int rhsIdx) const;

private:
    real *dm_a[NumLevels][NumLevels];
    real *dm_b[NumLevels][NumLevels];
    real *rhs[NumLevels][NumLevels][NumMultistep];

    bool a_is_old;
    unsigned int head;
};


/* functor classes */

class ISrcFunctor
{
public:
    ISrcFunctor() { }
    virtual ~ISrcFunctor() { }
    virtual real *operator()() const = 0;
};

class GetSrcField : public ISrcFunctor
{
private:
    real *m_address;

public:
    GetSrcField(real *address) :
	m_address(address)
    { }

    real *operator()() const
    {
	return m_address;
    }

};

/* TODO: complex results support ? */

class GetSrcDensity : public ISrcFunctor
{
private:
    DensityMatrix *m_dm;
    unsigned int m_row;
    unsigned int m_col;

public:
    GetSrcDensity(DensityMatrix *dm, unsigned int row, unsigned int col) :
	m_dm(dm), m_row(row), m_col(col)
    {
    }

    real *operator()() const
    {
	return m_dm->OldDM(m_row, m_col);
    }

};

/* TODO: CopyListEntry base class */
/* subclasses for fields and dm? */
/* specialize functions getDst, getSrc? */
/* OR use generic programming/templates */

class CopyListEntry
{
private:
    ISrcFunctor* m_srcFunctor;
    Result *m_res;
    unsigned int m_size;
    unsigned int m_interval;
    unsigned int m_position;

    /* TODO: base address + offset (position */

public:
    CopyListEntry(ISrcFunctor* srcFunctor, Result *result, unsigned int count,
		  unsigned int position, unsigned int interval) :
	m_srcFunctor(srcFunctor), m_res(result), m_size(sizeof(real) * count),
	m_position(position), m_interval(interval)
    {
    }

    ~CopyListEntry()
    {
	delete m_srcFunctor;
    }

    real *getSrc() const { return (*m_srcFunctor)() + m_position; }

    real *getDst(unsigned int idx) const {
	return m_res->data(idx/m_interval);
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
