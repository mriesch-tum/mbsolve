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

#ifndef MBSOLVE_SOLVER_CUDA_2LVL_PC_H
#define MBSOLVE_SOLVER_CUDA_2LVL_PC_H

#include <cuda.h>
#include <copy_list_entry_int.hpp>
#include <solver.hpp>

//#include <CUDADensityMatrix.hpp>

namespace mbsolve {

/* TODO: make general? */
static const unsigned int MAX_MATERIALS = 8;

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

class solver_cuda_2lvl_pc : public solver_int
{
public:
    solver_cuda_2lvl_pc(std::shared_ptr<const device> dev,
                        std::shared_ptr<scenario> scen);

    ~solver_cuda_2lvl_pc();

    const std::string& get_name() const;

    void run() const;

private:
    cudaStream_t comp_maxwell;
    cudaStream_t copy;

    real *m_h;
    real *m_e;
    real *m_d;

    unsigned int *m_mat_indices;



    std::vector<copy_list_entry> m_copy_list;
};

}

#endif
