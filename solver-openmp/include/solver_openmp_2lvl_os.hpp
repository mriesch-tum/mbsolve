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

#ifndef MBSOLVE_SOLVER_OPENMP_2LVL_OS_H
#define MBSOLVE_SOLVER_OPENMP_2LVL_OS_H

#include <Eigen/Core>
#include <solver.hpp>
#include <internal/common_fdtd_2lvl.hpp>
#include <internal/copy_list_entry.hpp>

namespace mbsolve {

class sim_constants_2lvl_os
{
public:
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
    real equi_inv;

    /* operator splitting */
    Eigen::Matrix3d prop_U02;
    Eigen::Matrix3d L_1E;

    /* simulation constants */
    real d_x_inv;
    real d_t;

    /* initialization constants */
    real inversion_init;
};


/**
 * OpenMP solver for 2-lvl systems using the FDTD scheme and the operator
 * splitting approach.
 * \ingroup MBSOLVE_SOLVER_OPENMP
 */
class solver_openmp_2lvl_os : public solver_int
{
public:
    solver_openmp_2lvl_os(std::shared_ptr<const device> dev,
                          std::shared_ptr<scenario> scen);

    ~solver_openmp_2lvl_os();

    const std::string& get_name() const;

    void run() const;

private:

    /* TODO: rule of three. make copy constructor etc. private?
     * or implement correctly
     */

    /*
     * Position-dependent density matrix in Liouville space. The population
     * inversion replaces the two populations. Order:
     * m_d[i][0] ... real part of rho12
     * m_d[i][1] ... imag part of rho12
     * m_d[i][2] ... population inversion rho11 - rho22
     */
    Eigen::Vector3d *m_d;


    real *m_h;
    real *m_e;

    real *m_result_scratch;

    real *m_source_data;

    unsigned int *m_mat_indices;

    //std::vector<sim_constants_2lvl> m_sim_consts;
    std::vector<sim_constants_2lvl_os> m_sim_consts;

    std::vector<sim_source> m_sim_sources;

    std::vector<copy_list_entry> m_copy_list;
};

}

#endif
