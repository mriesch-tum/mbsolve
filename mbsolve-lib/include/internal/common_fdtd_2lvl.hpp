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

#ifndef MBSOLVE_COMMON_FDTD_2LVL
#define MBSOLVE_COMMON_FDTD_2LVL

#include <memory>
#include <Eigen/Core>
#include <types.hpp>
#include <device.hpp>
#include <scenario.hpp>
#include <source.hpp>

namespace mbsolve {

/**
 * This class stores the material and simulation properties in concentrated
 * form.
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB
 */
class sim_constants_2lvl
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

    /* simulation constants */
    real d_x_inv;
    real d_t;

    /* initialization constants */
    real inversion_init;
};

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
    Eigen::Vector3d equi_corr;

    /* simulation constants */
    real d_x_inv;
    real d_t;

    /* initialization constants */
    real inversion_init;
};

/**
 * This class represents a source in concentrated form.
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB
 */
class sim_source
{
public:
    source::type type;
    unsigned int x_idx;
    unsigned int data_base_idx;
};

void
init_fdtd_simulation(std::shared_ptr<const device> dev,
                     std::shared_ptr<scenario> scen,
                     real courant = 1);

std::vector<sim_constants_2lvl>
init_sim_constants(std::shared_ptr<const device> dev,
                   std::shared_ptr<const scenario> scen,
                   std::map<std::string, unsigned int>& id_to_idx);


}

#endif
