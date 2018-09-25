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

#include <iostream>
#include <cmath>
#include <omp.h>
#include <solver_openmp_2lvl_pc.hpp>
#include <common_openmp.hpp>

namespace mbsolve {

static solver_factory<solver_openmp_2lvl_pc> factory("openmp-2lvl-pc");

solver_openmp_2lvl_pc::solver_openmp_2lvl_pc(std::shared_ptr<const device> dev,
                                             std::shared_ptr<scenario> scen) :
solver_int(dev, scen)
{
    /* TODO: scenario, device sanity check */
    /*
     * device.length > 0 (-> regions.size() > 0)
     * required materials found?
     * no gap in regions
     *
     */

    /* TODO: solver params
     * courant number
     */

    if (dev->get_regions().size() == 0) {
        throw std::invalid_argument("No regions in device!");
    }

    /* determine simulation settings */
    init_fdtd_simulation(dev, scen, 0.5);

    unsigned int num_gridpoints = m_scenario->get_num_gridpoints();
    unsigned int num_timesteps = m_scenario->get_num_timesteps();

    /* set up simulaton constants */
    std::map<std::string, unsigned int> id_to_idx;
    m_sim_consts = init_sim_constants(dev, scen, id_to_idx);

    /* allocate data arrays */
    m_inv = (real *) mb_aligned_alloc(sizeof(real) * num_gridpoints);
    m_dm12r = (real *) mb_aligned_alloc(sizeof(real) * num_gridpoints);
    m_dm12i = (real *) mb_aligned_alloc(sizeof(real) * num_gridpoints);
    m_h = (real *) mb_aligned_alloc(sizeof(real) * (num_gridpoints + 1));
    m_e = (real *) mb_aligned_alloc(sizeof(real) * num_gridpoints);
    m_mat_indices = (unsigned int *)
        mb_aligned_alloc(sizeof(unsigned int) * num_gridpoints);

    /* set up results and transfer data structures */
    unsigned int scratch_size = 0;
    for (const auto& rec : scen->get_records()) {
        /* create copy list entry */
        copy_list_entry entry(rec, scen, scratch_size);

        /* add result to solver */
        m_results.push_back(entry.get_result());

        /* calculate scratch size */
        scratch_size += entry.get_size();

        /* take imaginary part into account */
        if (rec->is_complex()) {
            scratch_size += entry.get_size();
        }

        /* TODO check if result is available */
        /*
           throw std::invalid_argument("Requested result is not available!");
        */

        m_copy_list.push_back(entry);
    }

    /* allocate scratchpad result memory */
    m_result_scratch = (real *) mb_aligned_alloc(sizeof(real) * scratch_size);
    m_scratch_size = scratch_size;

    /* create source data */
    m_source_data = new real[scen->get_num_timesteps() *
                             scen->get_sources().size()];
    unsigned int base_idx = 0;
    for (const auto& src : scen->get_sources()) {
        sim_source s;
        s.type = src->get_type();
        s.x_idx = src->get_position()/scen->get_gridpoint_size();
        s.data_base_idx = base_idx;
        m_sim_sources.push_back(s);

        /* calculate source values */
        for (unsigned int j = 0; j < scen->get_num_timesteps(); j++) {
            m_source_data[base_idx + j] =
                src->get_value(j * scen->get_timestep_size());
        }

        base_idx += scen->get_num_timesteps();
    }

    /* determine material indices of each grid point */
    unsigned int *l_mat_indices = new unsigned int[scen->get_num_gridpoints()];
    for (unsigned int i = 0; i < scen->get_num_gridpoints(); i++) {
        unsigned int mat_idx = 0;
        real x = i * scen->get_gridpoint_size();

        for (const auto& reg : dev->get_regions()) {
            if ((x >= reg->get_x_start()) && (x <= reg->get_x_end())) {
                mat_idx = id_to_idx[reg->get_material()->get_id()];
                break;
            }
        }

       /* TODO handle region not found */

        l_mat_indices[i] = mat_idx;
    }

    /* initialize arrays in parallel */
#ifndef XEON_PHI_OFFLOAD
    l_copy_list = m_copy_list.data();
    l_sim_consts = m_sim_consts.data();
    l_sim_sources = m_sim_sources.data();
#else
    /* prepare to offload sources */
    unsigned int num_sources = m_sim_sources.size();
    l_sim_sources = new sim_source[num_sources];
    for (int i = 0; i < num_sources; i++) {
        l_sim_sources[i] = m_sim_sources[i];
    }

    /* prepare to offload simulation constants */
    l_sim_consts = new sim_constants_2lvl[m_sim_consts.size()];
    for (unsigned int i = 0; i < m_sim_consts.size(); i++) {
        l_sim_consts[i] = m_sim_consts[i];
    }

    /* prepare to offload copy list entries */
    unsigned int num_copy = m_copy_list.size();
    l_copy_list = new copy_list_entry_dev[num_copy];
    for (int i = 0; i < m_copy_list.size(); i++) {
        l_copy_list[i] = m_copy_list[i].get_dev();
    }

#pragma offload target(mic:0) in(num_sources, num_copy, num_gridpoints) \
    in(l_mat_indices:length(num_gridpoints))                            \
    in(l_copy_list:length(num_copy) __mb_phi_create)                    \
    in(m_source_data:length(num_timesteps * num_sources) __mb_phi_create) \
    in(l_sim_sources:length(num_sources) __mb_phi_create)               \
    in(l_sim_consts:length(m_sim_consts.size()) __mb_phi_create)        \
    inout(m_e,m_h,m_inv:length(num_gridpoints) __mb_phi_create)         \
    inout(m_dm12i,m_dm12r:length(num_gridpoints) __mb_phi_create)       \
    inout(m_mat_indices:length(num_gridpoints) __mb_phi_create)
    {
#endif
#pragma omp parallel for schedule(static)
        for (unsigned int i = 0; i < num_gridpoints; i++) {
            unsigned int idx = l_mat_indices[i];
            m_mat_indices[i] = idx;

            /* TODO: evaluate flexible initialization in scenario */
            m_inv[i] = l_sim_consts[idx].inversion_init;
            m_dm12r[i] = 0.0;
            m_dm12i[i] = 0.0;
            m_e[i] = 0.0;
            m_h[i] = 0.0;
            if (i == num_gridpoints - 1) {
                m_h[i + 1] = 0.0;
            }
        }
#ifdef XEON_PHI_OFFLOAD
    }
#endif

    delete[] l_mat_indices;
}

solver_openmp_2lvl_pc::~solver_openmp_2lvl_pc()
{
#ifdef XEON_PHI_OFFLOAD
    unsigned int num_gridpoints = m_scenario->get_num_gridpoints();
    unsigned int num_timesteps = m_scenario->get_num_timesteps();
    unsigned int num_sources = m_sim_sources.size();

#pragma offload target(mic:0) \
    in(l_copy_list:length(m_copy_list.size()) __mb_phi_delete)          \
    in(m_source_data:length(num_timesteps * num_sources) __mb_phi_delete) \
    in(l_sim_sources:length(num_sources) __mb_phi_delete)               \
    in(l_sim_consts:length(m_sim_consts.size()) __mb_phi_delete)        \
    inout(m_e,m_h,m_inv:length(num_gridpoints) __mb_phi_delete)         \
    inout(m_dm12i,m_dm12r:length(num_gridpoints) __mb_phi_delete)       \
    inout(m_mat_indices:length(num_gridpoints) __mb_phi_delete)
    {
    }

    delete[] l_copy_list;
    delete[] l_sim_consts;
    delete[] l_sim_sources;
#endif

    mb_aligned_free(m_h);
    mb_aligned_free(m_e);
    mb_aligned_free(m_inv);
    mb_aligned_free(m_dm12r);
    mb_aligned_free(m_dm12i);
    mb_aligned_free(m_mat_indices);
    mb_aligned_free(m_result_scratch);

    delete[] m_source_data;
}

const std::string&
solver_openmp_2lvl_pc::get_name() const
{
    return factory.get_name();
}

void
solver_openmp_2lvl_pc::run() const
{
    unsigned int num_gridpoints = m_scenario->get_num_gridpoints();
    unsigned int num_timesteps = m_scenario->get_num_timesteps();
    unsigned int num_sources = m_sim_sources.size();
    unsigned int num_copy = m_copy_list.size();

#ifdef XEON_PHI_OFFLOAD
#pragma offload target(mic:0)                                           \
    in(num_gridpoints, num_timesteps, num_sources, num_copy)            \
    in(l_copy_list:length(num_copy) __mb_phi_use)                       \
    in(m_source_data:length(num_timesteps * num_sources) __mb_phi_use)  \
    in(l_sim_sources:length(num_sources) __mb_phi_use)                  \
    in(l_sim_consts:length(m_sim_consts.size()) __mb_phi_use)           \
    in(m_e,m_h,m_inv,m_dm12r:length(num_gridpoints) __mb_phi_use)       \
    in(m_dm12i,m_mat_indices:length(num_gridpoints) __mb_phi_use)       \
    inout(m_result_scratch:length(m_scratch_size))
    {
#endif
#pragma omp parallel
        {
            __mb_assume_aligned(m_e);
            __mb_assume_aligned(m_h);
            __mb_assume_aligned(m_inv);
            __mb_assume_aligned(m_dm12r);
            __mb_assume_aligned(m_dm12i);
            __mb_assume_aligned(m_mat_indices);
            __mb_assume_aligned(m_result_scratch);

            /* main loop */
            for (int n = 0; n < num_timesteps; n++) {
                /* update dm and e in parallel */
                //#pragma omp for simd schedule(static)
#pragma omp for schedule(static)
                for (int i = 0; i < num_gridpoints; i++) {
                    unsigned int mat_idx = m_mat_indices[i];

                    real inv_e = m_inv[i];
                    real rho12r_e = m_dm12r[i];
                    real rho12i_e = m_dm12i[i];
                    real field_e = m_e[i];

                    for (int pc_step = 0; pc_step < 4; pc_step++) {
                        /* execute prediction - correction steps */

                        real inv  = 0.5 * (m_inv[i] + inv_e);
                        real rho12r = 0.5 * (m_dm12r[i] + rho12r_e);
                        real rho12i = 0.5 * (m_dm12i[i] + rho12i_e);
                        real e = 0.5 * (m_e[i] + field_e);
                        real OmRabi = l_sim_consts[mat_idx].d12 * e;

                        inv_e = m_inv[i] + l_sim_consts[mat_idx].d_t *
                            (- 4.0 * OmRabi * rho12i
                             - l_sim_consts[mat_idx].tau1 *
                             (inv - l_sim_consts[mat_idx].equi_inv));

                        rho12i_e = m_dm12i[i] + l_sim_consts[mat_idx].d_t *
                            (- l_sim_consts[mat_idx].w12 * rho12r
                             + OmRabi * inv
                             - l_sim_consts[mat_idx].gamma12 * rho12i);

                        rho12r_e = m_dm12r[i] + l_sim_consts[mat_idx].d_t *
                            (+ l_sim_consts[mat_idx].w12 * rho12i
                             - l_sim_consts[mat_idx].gamma12 * rho12r);

                        real j = l_sim_consts[mat_idx].sigma * e;

                        real p_t = l_sim_consts[mat_idx].M_CP
                            * l_sim_consts[mat_idx].d12 *
                            (l_sim_consts[mat_idx].w12 * rho12i -
                             l_sim_consts[mat_idx].gamma12 * rho12r);

                        field_e = m_e[i] + l_sim_consts[mat_idx].M_CE *
                            (-j - p_t + (m_h[i + 1] - m_h[i])
                             * l_sim_consts[mat_idx].d_x_inv);
                    }

                    /* final update step */
                    m_inv[i] = inv_e;
                    m_dm12i[i] = rho12i_e;
                    m_dm12r[i] = rho12r_e;

                    m_e[i] = field_e;
                }

                /* apply sources */
                for (unsigned int k = 0; k < num_sources; k++) {
                    sim_source src = l_sim_sources[k];
                    if (src.type == source::type::hard_source) {
                        m_e[src.x_idx] = m_source_data[src.data_base_idx + n];
                    } else if (src.type == source::type::soft_source) {
                        m_e[src.x_idx] += m_source_data[src.data_base_idx + n];
                    } else {
                    }
                }

                /* update h in parallel */
                //#pragma omp for simd schedule(static)
#pragma omp for schedule(static)
                for (int i = 1; i < num_gridpoints; i++) {
                    unsigned int mat_idx = m_mat_indices[i - 1];

                    m_h[i] += l_sim_consts[mat_idx].M_CH *
                        (m_e[i] - m_e[i - 1]);
                }

                /* save results to scratchpad in parallel */
                for (int k = 0; k < num_copy; k++) {
                    if (l_copy_list[k].hasto_record(n)) {
                        unsigned int pos = l_copy_list[k].get_position();
                        unsigned int cols = l_copy_list[k].get_cols();
                        record::type t = l_copy_list[k].get_type();
                        int o_r = l_copy_list[k].get_offset_scratch_real(n, 0);

                        real *src_real;
                        if (t == record::type::electric) {
                            src_real = m_e;
                        } else if (t == record::type::inversion) {
                            src_real = m_inv;
                        } else {
                            /* TODO handle trouble */
                        }

                        //#pragma omp for simd schedule(static)
#pragma omp for schedule(static)
                        for (int i = pos; i < pos + cols; i++) {
                            m_result_scratch[o_r + i - pos] = src_real[i];
                        }

                        if (l_copy_list[k].is_complex()) {
                            int o_i = l_copy_list[k].get_offset_scratch_imag
                                (n, 0);

                            /* TODO handle imag part */
                        }
                    }
                }
            }
        }
#ifdef XEON_PHI_OFFLOAD
    }
#endif

    /* bulk copy results into result classes */
    for (const auto& cle : m_copy_list) {
        real *dr = m_result_scratch + cle.get_offset_scratch_real(0, 0);
        std::copy(dr, dr + cle.get_size(), cle.get_result_real(0, 0));
        if (cle.is_complex()) {
            real *di = m_result_scratch + cle.get_offset_scratch_imag(0, 0);
            std::copy(di, di + cle.get_size(), cle.get_result_imag(0, 0));
        }
    }
}

}
