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

#include <solver_openmp_2lvl_pc_red.hpp>

namespace mbsolve{

static solver_factory<solver_openmp_2lvl_pc_red> factory("openmp-2lvl-pc-red");

/** TODO necessary??
extern struct sim_constants gsc[MaxRegions];
extern unsigned int num_grid_points;
extern unsigned int num_time_steps;
extern real time_step_size;
*/

/* redundant calculation overlap */
//const unsigned int OL = 100;
unsigned int OL;

solver_openmp_2lvl_pc_red::solver_openmp_2lvl_pc_red
(std::shared_ptr<const device> dev, std::shared_ptr<scenario> scen) :
solver_int(dev, scen)
{
    /* TODO: scenario, device sanity check */

    /* TODO: solver params
     * courant number
     * overlap
     */

    OL = 1;

    /* determine simulation settings */
    init_fdtd_simulation(dev, scen, 0.5);

    /* set up simulaton constants */
    std::map<std::string, unsigned int> id_to_idx;
    m_sim_consts = init_sim_constants(dev, scen, id_to_idx);

#if 0
    num_grid_points = m_scenario.NumGridPoints;
    num_time_steps = m_scenario.NumTimeSteps;
    time_step_size = m_scenario.TimeStepSize;
#endif

    /* set up indices array and initialize data arrays */
    unsigned int P = omp_get_max_threads();
    unsigned int chunk = scen->get_num_gridpoints()/P;

    std::cout << "Number of threads: " << P << std::endl;
    m_inv = new real*[P];
    m_dm12r = new real*[P];
    m_dm12i = new real*[P];
    m_e = new real*[P];
    m_h = new real*[P];
    m_mat_indices = new unsigned int*[P];

    //#pragma offload target(mic) in(num_grid_points, chunk, P, OL)     \
    //inout(m_e,m_h,m_dm11,m_dm12i,m_dm12r,m_dm22,region_indices:length(P))
    //{
    for (int tid = 0; tid < P; tid++) {
        unsigned int size = chunk + 2 * OL;

        if (tid == P - 1) {
            //size += num_grid_points % P;
            size += scen->get_num_gridpoints() % P;
        }

        /* allocation */
        m_inv[tid] = new real[size];
        m_dm12r[tid] = new real[size];
        m_dm12i[tid] = new real[size];

        m_h[tid] = new real[size];
        m_e[tid] = new real[size];

        m_mat_indices[tid] = new unsigned int[size];
    }

    /* initialize memory in parallel */
#pragma omp parallel
    {
        unsigned int tid = omp_get_thread_num();
        //unsigned int chunk_base = num_grid_points/P;
        unsigned int chunk_base = scen->get_num_gridpoints()/P;
        unsigned int chunk = chunk_base;

        if (tid == P - 1) {
            //chunk += num_grid_points % P;
            chunk += scen->get_num_gridpoints() % P;
        }

        for (int i = 0; i < chunk + 2 * OL; i++) {
            unsigned int mat_idx = 0;
            unsigned int global_idx = tid * chunk_base + (i - OL);
            real x = global_idx * scen->get_gridpoint_size();

            for (const auto& reg : dev->get_regions()) {
                if ((x >= reg->get_start()) && (x <= reg->get_end())) {
                    mat_idx = id_to_idx[reg->get_material()->get_id()];
                    break;
                }
            }
            m_inv[tid][i] = m_sim_consts[mat_idx].inversion_init;
            m_mat_indices[tid][i] = mat_idx;
            m_dm12r[tid][i] = 0.0;
            m_dm12i[tid][i] = 0.0;
            m_e[tid][i] = 0.0;
            m_h[tid][i] = 0.0;
        }
#pragma omp barrier
    }
    //}

    /* set up results and transfer data structures */
    unsigned int scratch_size = 0;
    for (const auto& rec : scen->get_records()) {
        /* create copy list entry */
        copy_list_entry entry(rec, scen);

        /* add result to solver */
        m_results.push_back(entry.get_result());

        /* calculate scratch size */
        scratch_size += entry.get_size();

        /* TODO: make more generic? */
        /* TODO: move to parser in record class */
        /* add source address to copy list entry */
        if (rec->get_name() == "inv12") {
            entry.set_real(m_inv);
        } else if (rec->get_name() == "d12") {
            entry.set_real(m_dm12r);
            entry.set_imag(m_dm12i);

            /* take imaginary part into account */
            scratch_size += entry.get_size();
        } else if (rec->get_name() == "e") {
            entry.set_real(m_e);
        } else if (rec->get_name() == "h") {
            /* TODO: numGridPoints + 1 ? */
            entry.set_real(m_h);
        } else {
            throw std::invalid_argument("Requested result is not available!");
        }

        m_copy_list.push_back(entry);
    }

    /* allocate scratchpad result memory */
    m_result_scratch = new real[scratch_size];

    /* add scratchpad addresses to copy list entries */
    unsigned int scratch_offset = 0;
    for (auto& cle : m_copy_list) {

        std::cout << cle.get_position() << std::endl;

        cle.set_scratch_real(&m_result_scratch[scratch_offset]);
        scratch_offset += cle.get_size();

        if (cle.get_record()->get_name() == "d12") {
            /* complex result */
            cle.set_scratch_imag(&m_result_scratch[scratch_offset]);
            scratch_offset += cle.get_size();
        }
    }

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
}

solver_openmp_2lvl_pc_red::~solver_openmp_2lvl_pc_red()
{
    unsigned int P = omp_get_max_threads();

    //#pragma offload target(mic) in(P)                                 \
    //inout(m_e,m_h,m_dm11,m_dm12i,m_dm12r,m_dm22,region_indices:length(P))
    //{
    for (int tid = 0; tid < P; tid++) {
        delete[] m_h[tid];
        delete[] m_e[tid];

        delete[] m_inv[tid];
        delete[] m_dm12r[tid];
        delete[] m_dm12i[tid];

        delete[] m_mat_indices[tid];
    }

    delete[] m_h;
    delete[] m_e;
    delete[] m_inv;
    delete[] m_dm12r;
    delete[] m_dm12i;
    delete[] m_mat_indices;
    delete[] m_result_scratch;
    delete[] m_source_data;
}

const std::string&
solver_openmp_2lvl_pc_red::get_name() const
{
    return factory.get_name();
}

void
solver_openmp_2lvl_pc_red::run() const
{
    unsigned int P = omp_get_max_threads();
    //#pragma offload target(mic) in(gsc, num_grid_points, num_time_steps) \
    //in(time_step_size, P)                                             \
    //inout(m_e,m_h,m_dm11,m_dm12i,m_dm12r,m_dm22,region_indices:length(P))
    {
#pragma omp parallel
        {
            unsigned int tid = omp_get_thread_num();
            unsigned int chunk_base = m_scenario->get_num_gridpoints()/P;
            unsigned int chunk = chunk_base;

            if (tid == P - 1) {
                //chunk += num_grid_points % P;
                chunk += m_scenario->get_num_gridpoints() % P;
            }

            /* main loop */
            //for (unsigned int n = 0; n < num_time_steps/OL; n++) {
            for (unsigned int n = 0; n < m_scenario->get_num_timesteps()/OL;
                 n++) {
                /* exchange data */
                if (tid > 0) {
#pragma ivdep
                    for (unsigned int i = 0; i < OL; i++) {
                        m_inv[tid][i] = m_inv[tid - 1][chunk + i];
                        m_dm12r[tid][i] = m_dm12r[tid - 1][chunk + i];
                        m_dm12i[tid][i] = m_dm12i[tid - 1][chunk + i];
                        m_e[tid][i] = m_e[tid - 1][chunk + i];
                        m_h[tid][i] = m_h[tid - 1][chunk + i];
                    }
                }

                if (tid < P - 1) {
#pragma ivdep
                    for (unsigned int i = 0; i < OL; i++) {
                        m_inv[tid][OL + chunk + i] = m_inv[tid + 1][OL + i];
                        m_dm12r[tid][OL + chunk + i] = m_dm12r[tid + 1][OL + i];
                        m_dm12i[tid][OL + chunk + i] = m_dm12i[tid + 1][OL + i];
                        m_e[tid][OL + chunk + i] = m_e[tid + 1][OL + i];
                        m_h[tid][OL + chunk + i] = m_h[tid + 1][OL + i];
                    }
                }

                /* sync after communication */
#pragma omp barrier

                /* sub-loop */
                for (unsigned int m = 0; m < OL; m++) {

                    /* update dm and e */
#pragma simd
                    for (int i = m; i < chunk + 2 * OL - m - 1; i++) {
                        int mat_idx = m_mat_indices[tid][i];

                        real inv_e = m_inv[tid][i];
                        real rho12r_e = m_dm12r[tid][i];
                        real rho12i_e = m_dm12i[tid][i];
                        real field_e = m_e[tid][i];

                        for (int pc_step = 0; pc_step < 4; pc_step++) {
                            /* execute prediction - correction steps */

                            real inv  = 0.5 * (m_inv[tid][i] + inv_e);
                            real rho12r = 0.5 * (m_dm12r[tid][i] + rho12r_e);
                            real rho12i = 0.5 * (m_dm12i[tid][i] + rho12i_e);
                            real e = 0.5 * (m_e[tid][i] + field_e);
                            real OmRabi = 0.5 * m_sim_consts[mat_idx].d12 * e;

                            inv_e = m_inv[tid][i] + m_sim_consts[mat_idx].d_t *
                                (- 2.0 * OmRabi * rho12i
                                 - m_sim_consts[mat_idx].tau1 * inv);

                            rho12i_e = m_dm12i[tid][i]
                                + m_sim_consts[mat_idx].d_t *
                                (- m_sim_consts[mat_idx].w12 * rho12r
                                 + OmRabi * inv
                                 - m_sim_consts[mat_idx].gamma12 * rho12i);

                            rho12r_e = m_dm12r[tid][i]
                                + m_sim_consts[mat_idx].d_t *
                                (+ m_sim_consts[mat_idx].w12 * rho12i
                                 - m_sim_consts[mat_idx].gamma12 * rho12r);

                            real j = m_sim_consts[mat_idx].sigma * e;

                            real p_t = m_sim_consts[mat_idx].M_CP
                                * m_sim_consts[mat_idx].d12 *
                                (m_sim_consts[mat_idx].w12 * rho12i -
                                 m_sim_consts[mat_idx].gamma12 * rho12r);

                            field_e = m_e[tid][i]
                                + m_sim_consts[mat_idx].M_CE *
                                (-j - p_t + (m_h[tid][i + 1] - m_h[tid][i]) *
                                 m_sim_consts[mat_idx].d_x_inv);
                        }

                        /* final update step */
                        m_inv[tid][i] = inv_e;
                        m_dm12i[tid][i] = rho12i_e;
                        m_dm12r[tid][i] = rho12r_e;
                        m_e[tid][i] = field_e;

                        /* apply sources */
                        for (const auto& src : m_sim_sources) {
                            if (tid * chunk_base + (i - OL) == src.x_idx) {
                                if (src.type == source::type::hard_source) {
                                    m_e[tid][i] =
                                        m_source_data[src.data_base_idx
                                                      + (n * OL + m)];
                                } else if (src.type ==
                                           source::type::soft_source) {
                                    m_e[tid][i] +=
                                        m_source_data[src.data_base_idx
                                                      + (n * OL + m)];
                                } else {
                                }
                            }
                        }
                    }

                    /* update h */
#pragma ivdep
                    for (int i = m + 1; i < chunk + 2 * OL - m - 1; i++) {
                        int mat_idx = m_mat_indices[tid][i];

                        m_h[tid][i] += m_sim_consts[mat_idx].M_CH *
                            (m_e[tid][i] - m_e[tid][i - 1]);
                    }

                    /* apply boundary condition */
                    if (tid == 0) {
                        m_h[tid][OL] = 0;
                    }

                    /* save results to scratchpad in parallel */
                    for (const auto& cle : m_copy_list) {
                        if (cle.hasto_record(n * OL + m)) {
                            for (int i = OL; i < chunk + OL; i++) {
                                unsigned int pos = cle.get_position();
                                unsigned int idx = tid * chunk_base + (i - OL);
                                if ((idx >= pos) &&
                                    (idx < pos + cle.get_cols())) {
                                    *cle.get_scratch_real(n * OL + m,
                                                          i - pos) =
                                        *cle.get_real(i, tid);
                                    if (cle.is_complex()) {
                                        *cle.get_scratch_imag(n * OL + m,
                                                              i - pos) =
                                            *cle.get_imag(i, tid);
                                    }
                                }
                            }
                        }
                    }
                    /* sync after computation */
#pragma omp barrier
                }
            }
        }
    }

    /* bulk copy results into result classes */
    for (const auto& cle : m_copy_list) {
        std::copy(cle.get_scratch_real(0, 0), cle.get_scratch_real(0, 0) +
                  cle.get_size(), cle.get_result_real(0, 0));
        if (cle.is_complex()) {
            std::copy(cle.get_scratch_imag(0, 0), cle.get_scratch_imag(0, 0) +
                      cle.get_size(), cle.get_result_imag(0, 0));
        }
    }
}

}
