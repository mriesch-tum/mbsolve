/*
 * mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
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

#ifndef MBSOLVE_COMMON_FDTD_H
#define MBSOLVE_COMMON_FDTD_H

#include <limits>
#include <memory>
#include <mbsolve/lib/device.hpp>
#include <mbsolve/lib/scenario.hpp>
#include <mbsolve/lib/source.hpp>
#include <mbsolve/lib/types.hpp>

namespace mbsolve {

/**
 * This class represents a source in concentrated form.
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB_INT
 */
class sim_source
{
public:
    source::type type;
    unsigned int x_idx;
    unsigned int data_base_idx;
};

/**
 * Initializes the FDTD simulation (grid point size and time step size
 * calculation).
 * \ingroup MBSOLVE_LIB_INT
 */
static void
init_fdtd_simulation(
    std::shared_ptr<const device> dev,
    std::shared_ptr<scenario> scen,
    real courant = 1)
{
    if (scen->get_num_gridpoints() > 0) {
        /* speed of light (use smallest value of relative permittivities) */
        real velocity =
            1.0 / sqrt(MU0 * EPS0 * dev->get_minimum_permittivity());

        /* get number of grid points */
        unsigned int n_x = scen->get_num_gridpoints();

        /* number of time steps */
        unsigned int n_t;

        /* special scenario for 0D simulations? */
        if (n_x == 1) {
            /*
             * if device length does not equal 0, this is probably not
             * intended!
             */
            if (dev->get_length() > 0.0) {
                std::cout << "Warning: Device with length = "
                          << dev->get_length()
                          << " simulated using only one grid point!"
                          << std::endl;
            }

            /*
             * the single grid point is infinitely large
             * note: It is not practical to use infinity, hence we use
             * approx. 1 light-year.
             */
            scen->set_gridpoint_size(1e13);

            /* use the number of time steps specified by user */
            n_t = scen->get_num_timesteps();
            if (n_t < 2) {
                throw std::invalid_argument("Invalid scenario.");
            }
        } else {
            /* divide device in equidistant grid */
            real d_x = dev->get_length() / (n_x - 1);
            scen->set_gridpoint_size(d_x);

            /* get time step size via Courant number and max. velocity */
            real d_t = courant * d_x / velocity;

            /* number of time steps */
            n_t = ceil(scen->get_endtime() / d_t) + 1;
            scen->set_num_timesteps(n_t);
        }

        /* (re)-adjust time step size in order to fit number of time steps */
        if (scen->get_endtime() > std::numeric_limits<real>::min()) {
            real d_t = scen->get_endtime() / (n_t - 1);
            scen->set_timestep_size(d_t);
        } else {
            scen->set_num_timesteps(1);
            /*
             * the single time step is infinitely large
             * note: It is not practical to use infinity, hence we use
             * approx. 1 century
             */
            scen->set_timestep_size(3.154e10);
        }
    } else {
        throw std::invalid_argument("Invalid scenario.");
    }
}

/**
 * This class stores the material and simulation properties for the FDTD
 * method in precomputed form.
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB_INT
 */
class sim_constants_fdtd
{
public:
    /* prefactor for current electric field in Ampere's law (update equation
     * for electric field).
     */
    real fac_a;

    /* prefactor for the \partial_x H_y and \partial_t P_z (update equation
     * for electric field)
     */
    real fac_b;

    /* prefactor for \partial_x E_z (update equation for magnetic field) */
    real fac_c;

    /* overlap factor (applied as pre-factor to polarization, which is
     * determined by the density matrix)
     */
    real gamma;
};

/**
 * Precomputes the factors for the electric and magnetic field update
 * equation, which are stored in \ref sim_constants_fdtd.
 * \ingroup MBSOLVE_LIB_INT
 */
static sim_constants_fdtd
get_fdtd_constants(
    std::shared_ptr<const device> dev,
    std::shared_ptr<const scenario> scen,
    std::shared_ptr<const material> mat)
{
    sim_constants_fdtd sc;

    /* convert loss term to conductivity */
    real sigma = sqrt(
                     EPS0 * mat->get_rel_permittivity() /
                     (MU0 * mat->get_rel_permeability())) *
        mat->get_losses() * 2.0;

    real dt = scen->get_timestep_size();
    real epsilon = EPS0 * mat->get_rel_permittivity();
    real x = (sigma * dt) / (2.0 * epsilon);

    /* prefactor a */
    sc.fac_a = (1.0 - x) / (1.0 + x);

    /* prefactor b */
    sc.fac_b = dt / epsilon / (1.0 + x);

    /* prefactor c */
    sc.fac_c =
        dt / (MU0 * mat->get_rel_permeability() * scen->get_gridpoint_size());

    /* overlap factor */
    sc.gamma = mat->get_overlap_factor();

    return sc;
}
}

#endif
