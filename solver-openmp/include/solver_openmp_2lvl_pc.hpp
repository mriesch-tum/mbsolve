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

#ifndef MBSOLVE_SOLVER_OPENMP_2LVL_PC_H
#define MBSOLVE_SOLVER_OPENMP_2LVL_PC_H

#include <solver.hpp>

namespace mbsolve {

/**
 * \defgroup MBSOLVE_SOLVER_OPENMP solver-openmp
 * Different solvers that use OpenMP for parallelization.
 */

/**
 *
 * \ingroup MBSOLVE_SOLVER_OPENMP
 */
struct sim_constants_2lvl
{
    /* electromagnetic constants */
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    /* quantum mechanical constants */
    real w12;
    real d12;
    real tau1;
    real gamma12;

    /* required? */
    unsigned int idx_start;
    unsigned int idx_end;

    /* simulation constants */
    real d_x_inv;
    real d_t;

    /* initialization constants */
    real dm11_init;
    real dm22_init;
};
/*
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

    virtual real *getSrc() const = 0;

    real *getDst(unsigned int idx) const {
	return m_res->data(idx/m_interval);
    }

    unsigned int getSize() const { return sizeof(real) * m_count; }

    unsigned int getCount() const { return m_count; }

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

class solver_openmp_2lvl_pc : public solver_int
{
public:
    solver_openmp_2lvl_pc(std::shared_ptr<const device> dev,
                          std::shared_ptr<scenario> scen);

    ~solver_openmp_2lvl_pc();

    const std::string& get_name() const;

    void run() const;

private:
    inline void estimate_step(int i, real src) const;

    //DensityMatrixData *m_dm;

    real *m_dm11;
    real *m_dm12r;
    real *m_dm12i;
    real *m_dm22;

    real *m_h;
    real *m_e;
    real *m_e_est;

    unsigned int *region_indices;

    //    std::vector<CopyListEntry *> m_copyListBlack;
    //std::vector<CopyListEntry *> m_copyListRed;
};

}

#endif
