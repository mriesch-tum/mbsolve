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

#include <stdexcept>
#include <internal/common_fdtd_2lvl.hpp>

namespace mbsolve{

std::vector<sim_constants_2lvl>
init_sim_constants(std::shared_ptr<const device> dev,
                   std::shared_ptr<const scenario> scen,
                   std::map<std::string, unsigned int>& id_to_idx)
{
    std::vector<sim_constants_2lvl> sim_constants;

    unsigned int j = 0;

    for (const auto& mat_id : dev->get_used_materials()) {
        sim_constants_2lvl sc;

        auto mat = material::get_from_library(mat_id);

        /* factor for electric field update */
        sc.M_CE = scen->get_timestep_size()/
            (EPS0 * mat->get_rel_permittivity());

        /* factor for magnetic field update */
        sc.M_CH = scen->get_timestep_size()/
            (MU0 * mat->get_rel_permeability() * scen->get_gridpoint_size());

        /* convert loss term to conductivity */
        sc.sigma = sqrt(EPS0 * mat->get_rel_permittivity()/
                        (MU0 * mat->get_rel_permeability()))
            * mat->get_losses() * 2.0;

        /* active region in 2-lvl description? */
        /* TODO: remove ugly dynamic cast to qm_desc_2lvl */
        std::shared_ptr<qm_desc_2lvl> qm =
            std::dynamic_pointer_cast<qm_desc_2lvl>(mat->get_qm());
        if (qm) {
            /* factor for macroscopic polarization */
            sc.M_CP = -2.0 * HBAR * mat->get_overlap_factor() *
                qm->get_carrier_density();

            /* 2-lvl quantum mechanical system */
            sc.w12 = qm->get_transition_freq();
            sc.d12 = qm->get_dipole_moment() * E0/HBAR;
            sc.tau1 = qm->get_scattering_rate();
            sc.gamma12 = qm->get_dephasing_rate();
            sc.equi_inv = qm->get_equilibrium_inversion();

            if (scen->get_dm_init_type() == scenario::lower_full) {
                sc.inversion_init = -1.0;
            } else if (scen->get_dm_init_type() == scenario::upper_full) {
                sc.inversion_init = 1.0;
            } else {

            }
        } else {
            /* set all qm-related factors to zero */
            sc.M_CP = 0.0;
            sc.w12 = 0.0;
            sc.d12 = 0.0;
            sc.tau1 = 0.0;
            sc.gamma12 = 0.0;
            sc.equi_inv = 0.0;
            sc.inversion_init = 0.0;
        }

        /* simulation settings */
        sc.d_x_inv = 1.0/scen->get_gridpoint_size();
        sc.d_t = scen->get_timestep_size();

        sim_constants.push_back(sc);
        id_to_idx[mat->get_id()] = j;
        j++;
    }

    return sim_constants;
}

void init_fdtd_simulation(std::shared_ptr<const device> dev,
                          std::shared_ptr<scenario> scen,
                          real courant)
{
    if (scen->get_num_gridpoints() > 0) {
        /* speed of light (use smallest value of relative permittivities) */
        real velocity = 1.0/sqrt(MU0 * EPS0 * dev->get_minimum_permittivity());

        /* get number of grid points */
        unsigned int n_x = scen->get_num_gridpoints();

        /* grid point size */
        real d_x = dev->get_length()/(n_x - 1);
        scen->set_gridpoint_size(d_x);

        /* time step size */
        real d_t = courant * d_x/velocity;

        /* number of time steps */
        unsigned int n_t = ceil(scen->get_endtime()/d_t) + 1;
        scen->set_num_timesteps(n_t);

        /* re-adjust time step size in order to fit number of time steps */
        d_t = scen->get_endtime()/(n_t - 1);
        scen->set_timestep_size(d_t);

    } else {
        throw std::invalid_argument("Invalid scenario.");
    }
}

}
