#ifndef OPENMP_2LVL_PC_RED_H
#define OPENMP_2LVL_PC_RED_H

#include <iostream>
#include <omp.h>
#include <Solver.hpp>

namespace mbsolve {

const unsigned int MaxRegions = 8;

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

class CopyListEntry
{
protected:
    Result *m_res;
    unsigned int m_count;
    unsigned int m_interval;
    unsigned int m_position;

public:
    CopyListEntry(Result *result, unsigned int count,
		  unsigned int position, unsigned int interval) :
	m_res(result), m_count(count),
	m_position(position), m_interval(interval)
    {
    }

    ~CopyListEntry()
    {
    }

    virtual real *getSrc(unsigned int tid) const = 0;

    real *getDst(unsigned int idx) const {
	return m_res->data(idx/m_interval);
    }

    unsigned int getSize() const { return sizeof(real) * m_count; }

    unsigned int getCount() const { return m_count; }

    bool record(unsigned int idx) const { return (idx % m_interval) == 0; }
};

class CLEField_red : public CopyListEntry
{
private:
    real **m_address;
public:
    CLEField_red(real **address, Result *result, unsigned int count,
	     unsigned int position, unsigned int interval) :
	CopyListEntry(result, count, position, interval), m_address(address)
    {
    }

    real *getSrc(unsigned int tid) const
    {
	//std::cout << tid << std::endl;
	return m_address[tid];
    }
};
/*
class CLEDensity : public CopyListEntry
{
private:
    DensityMatrixData *m_dm;
    unsigned int m_row;
    unsigned int m_col;

public:
    CLEDensity(DensityMatrixData *dm, unsigned int row, unsigned int col,
	     Result *result, unsigned int count,
	     unsigned int position, unsigned int interval) :
	CopyListEntry(result, count, position, interval), m_dm(dm), m_row(row),
	m_col(col)
    {
    }

    real *getSrc() const
    {
	return m_dm->oldDM(m_row, m_col) + m_position;
    }
    };*/

class SolverOMP_2lvl_pc_red : public ISolver
{
public:
    SolverOMP_2lvl_pc_red(const Device& device, const Scenario& scenario);

    ~SolverOMP_2lvl_pc_red();

    std::string getName() const;

    void run() const;

private:

    real **m_dm11;
    real **m_dm12r;
    real **m_dm12i;
    real **m_dm22;

    real **m_h;
    real **m_e;

    unsigned int **region_indices;

    std::vector<CopyListEntry *> m_copyList;
};

}

#endif
