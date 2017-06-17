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

#ifndef MBSOLVE_SOLVER_OPENMP_2LVL_PC_RED_H
#define MBSOLVE_SOLVER_OPENMP_2LVL_PC_RED_H

#include <iostream>
#include <omp.h>
#include <copy_list_entry_int.hpp>
#include <solver.hpp>
#include <internal/common_fdtd_2lvl.hpp>

namespace mbsolve {

class solver_openmp_2lvl_pc_red : public solver_int
{
public:
    solver_openmp_2lvl_pc_red(std::shared_ptr<const device> dev,
                              std::shared_ptr<scenario> scen);

    ~solver_openmp_2lvl_pc_red();

    const std::string& get_name() const;

    void run() const;

private:
    /* TODO: rule of three. make copy constructor etc. private?
     * or implement correctly
     */

    real **m_inv;
    real **m_dm12r;
    real **m_dm12i;

    real **m_h;
    real **m_e;

    real *m_result_scratch;

    real *m_source_data;

    unsigned int **m_mat_indices;

    std::vector<sim_constants_2lvl> m_sim_consts;

    std::vector<sim_source> m_sim_sources;

    std::vector<copy_list_entry> m_copy_list;
};

}

#endif
