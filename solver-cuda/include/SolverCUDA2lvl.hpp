#ifndef SOLVERCUDA2LVL_H
#define SOLVERCUDA2LVL_H

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

    real d_x;
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
protected:
    Result *m_res;
    unsigned int m_size;
    unsigned int m_interval;
    unsigned int m_position;

public:
    CopyListEntry(Result *result, unsigned int count,
		  unsigned int position, unsigned int interval) :
	m_res(result), m_size(sizeof(real) * count),
	m_position(position), m_interval(interval)
    {
    }

    ~CopyListEntry()
    {
    }

    virtual real *getSrc() const = 0;

    real *getDst(unsigned int idx) const {
	return m_res->data(idx/m_interval);
    }

    unsigned int getSize() const { return m_size; }

    bool record(unsigned int idx) const { return (idx % m_interval) == 0; }
};

class CLEField : public CopyListEntry
{
private:
    real *m_address;
public:
    CLEField(real *address, Result *result, unsigned int count,
	     unsigned int position, unsigned int interval) :
	CopyListEntry(result, count, position, interval), m_address(address)
    {
    }

    real *getSrc() const
    {
	return m_address + m_position;
    }
};

class CLEDensity : public CopyListEntry
{
private:
    CUDADensityMatrix *m_dm;
    unsigned int m_row;
    unsigned int m_col;

public:
    CLEDensity(CUDADensityMatrix *dm, unsigned int row, unsigned int col,
	     Result *result, unsigned int count,
	     unsigned int position, unsigned int interval) :
	CopyListEntry(result, count, position, interval), m_dm(dm), m_row(row),
	m_col(col)
    {
    }

    real *getSrc() const
    {
	return m_dm->getData().oldDM(m_row, m_col) + m_position;
    }
};


class SolverCUDA2lvl : public ISolver
{
public:
    SolverCUDA2lvl(const Device& device, const Scenario& scenario);

    ~SolverCUDA2lvl();

    std::string getName() const;

    void run() const;

private:
    cudaStream_t comp_maxwell;
    cudaStream_t comp_density;
    cudaStream_t copy;

    CUDADensityMatrix *m_dm;

    real *m_h;
    real *m_e;

    std::vector<CopyListEntry *> m_copyListBlack;
    std::vector<CopyListEntry *> m_copyListRed;
};

}

#endif
