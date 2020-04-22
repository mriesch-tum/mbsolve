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

#ifndef MBSOLVE_SOLVER_CPU_FDTD_RED_TPP
#define MBSOLVE_SOLVER_CPU_FDTD_RED_TPP

#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_MALLOC
#define EIGEN_STRONG_INLINE inline

#include <iostream>
#include <random>
#include <omp.h>
#include <mbsolve/solver-cpu/internal/common_cpu.hpp>
#include <mbsolve/solver-cpu/solver_cpu_fdtd_red.hpp>

namespace mbsolve {

/* TODO make this flexible and useful */
const unsigned int VEC = 4;

template<unsigned int num_lvl, template<unsigned int> class density_algo>
solver_cpu_fdtd_red<num_lvl, density_algo>::solver_cpu_fdtd_red(
    std::shared_ptr<const device> dev,
    std::shared_ptr<scenario> scen)
  : solver(
        "cpu-fdtd-red-" + std::to_string(num_lvl) + "lvl-" +
            density_algo<num_lvl>::name(),
        dev,
        scen),
    OL((num_lvl == 6) ? 1 : 8)
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

    /* report overlap parameter */
    std::cout << "OpenMP overlap parameter: " << OL << std::endl;

    /* report number of used OpenMP threads */
    unsigned int P = omp_get_max_threads();
    std::cout << "Number of threads: " << P << std::endl;

    if (dev->get_regions().size() == 0) {
        throw std::invalid_argument("No regions in device!");
    }

    /* determine simulation settings */
    init_fdtd_simulation(dev, scen, 0.5);

    /* inverse grid point size for Ampere's law */
    m_dx_inv = 1.0 / scen->get_gridpoint_size();

    /* set up simulation constants for quantum mechanical description */
    std::vector<sim_constants_fdtd> m_sim_consts_fdtd;
    std::map<std::string, unsigned int> id_to_idx;
    unsigned int j = 0;
    for (const auto& mat_id : dev->get_used_materials()) {
        auto mat = material::get_from_library(mat_id);

        /* constants for electromagnetic (optical) field */
        m_sim_consts_fdtd.push_back(get_fdtd_constants(dev, scen, mat));

        /* constants for quantum mechanical system */
        m_sim_consts_qm.push_back(density_algo<num_lvl>::get_qm_constants(
            mat->get_qm(), scen->get_timestep_size()));

        /* enter material index */
        id_to_idx[mat->get_id()] = j;
        j++;
    }

    /* partition data in chunks per thread */
    uint64_t num_gridpoints = m_scenario->get_num_gridpoints();
    uint64_t chunk_base = m_scenario->get_num_gridpoints() / P;
    uint64_t chunk_rem = m_scenario->get_num_gridpoints() % P;
    uint64_t num_timesteps = m_scenario->get_num_timesteps();

    /* allocate data arrays */
    typedef typename density_algo<num_lvl>::density density_t;
    m_d = new density_t*[P];
    m_e = new real*[P];
    m_h = new real*[P];
    m_p = new real*[P];
    m_fac_a = new real*[P];
    m_fac_b = new real*[P];
    m_fac_c = new real*[P];
    m_gamma = new real*[P];
    m_mat_indices = new unsigned int*[P];

    for (unsigned int tid = 0; tid < P; tid++) {
        unsigned int chunk = chunk_base;

        if (tid == P - 1) {
            chunk += chunk_rem;
        }

        /* allocation */
        uint64_t size = chunk + 2 * OL;

        m_d[tid] = new density_t[size];
        m_h[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_e[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_p[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_fac_a[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_fac_b[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_fac_c[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_gamma[tid] = (real*) mb_aligned_alloc(size * sizeof(real));
        m_mat_indices[tid] =
            (unsigned int*) mb_aligned_alloc(size * sizeof(unsigned int));
    }

    /* set up indices array and initialize data arrays */
#pragma omp parallel
    {
        unsigned int tid = omp_get_thread_num();
        uint64_t chunk = chunk_base;

        if (tid == P - 1) {
            chunk += chunk_rem;
        }

        /* random number generation */
        std::random_device rand_dev;
        std::mt19937 rand_gen(rand_dev());
        std::normal_distribution<real> dis(0.0, 1.0);

        /* allocation */
        uint64_t size = chunk + 2 * OL;

        for (int i = 0; i < size; i++) {
            int64_t global_idx = tid * chunk_base + i - OL;
            real x = global_idx * scen->get_gridpoint_size();

            /* find region and corresponding material */
            int mat_idx = -1;
            bool has_qm = false;
            if ((global_idx >= 0) && (global_idx < num_gridpoints)) {
                for (const auto& reg : dev->get_regions()) {
                    if ((x >= reg->get_x_start()) &&
                        (x <= reg->get_x_end())) {
                        mat_idx = id_to_idx[reg->get_material()->get_id()];
                        has_qm =
                            (reg->get_material()->get_qm()) ? true : false;
                        break;
                    }
                }
            }

            /* initialize data arrays */
            if (mat_idx >= 0) {
                if (has_qm) {
                    auto ic_dm = scen->get_ic_density();
                    m_d[tid][i] = density_algo<num_lvl>::get_density(
                        ic_dm->initialize(x));
                } else {
                    m_d[tid][i] = density_algo<num_lvl>::get_density();
                }
                m_e[tid][i] = dis(rand_gen) * 1e-15;
                m_h[tid][i] = 0.0;
                m_p[tid][i] = 0.0;
                m_fac_a[tid][i] = m_sim_consts_fdtd[mat_idx].fac_a;
                m_fac_b[tid][i] = m_sim_consts_fdtd[mat_idx].fac_b;
                m_fac_c[tid][i] = m_sim_consts_fdtd[mat_idx].fac_c;
                m_gamma[tid][i] = m_sim_consts_fdtd[mat_idx].gamma;
                m_mat_indices[tid][i] = mat_idx;
            } else {
                m_d[tid][i] = density_algo<num_lvl>::get_density();
                m_e[tid][i] = 0.0;
                m_h[tid][i] = 0.0;
                m_p[tid][i] = 0.0;
                m_fac_a[tid][i] = 0.0;
                m_fac_b[tid][i] = 0.0;
                m_fac_c[tid][i] = 0.0;
                m_gamma[tid][i] = 0.0;
                m_mat_indices[tid][i] = 0;
            }
        }
#pragma omp barrier
    }

    /* set up results and transfer data structures */
    uint64_t scratch_size = 0;
    for (const auto& rec : scen->get_records()) {
        /* create copy list entry */
        copy_list_entry entry(rec, scen, scratch_size);

        std::cout << "Rows: " << entry.get_rows() << " "
                  << "Cols: " << entry.get_cols() << std::endl;

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
    m_result_scratch = (real*) mb_aligned_alloc(sizeof(real) * scratch_size);

    /* create source data */
    m_source_data =
        new real[scen->get_num_timesteps() * scen->get_sources().size()];
    unsigned int base_idx = 0;
    for (const auto& src : scen->get_sources()) {
        sim_source s;
        s.type = src->get_type();
        s.x_idx = src->get_position() / scen->get_gridpoint_size();
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

template<unsigned int num_lvl, template<unsigned int> class density_algo>
solver_cpu_fdtd_red<num_lvl, density_algo>::~solver_cpu_fdtd_red()
{
#pragma omp parallel
    {
        unsigned int tid = omp_get_thread_num();

        delete[] m_d[tid];
        mb_aligned_free(m_e[tid]);
        mb_aligned_free(m_h[tid]);
        mb_aligned_free(m_p[tid]);
        mb_aligned_free(m_fac_a[tid]);
        mb_aligned_free(m_fac_b[tid]);
        mb_aligned_free(m_fac_c[tid]);
        mb_aligned_free(m_gamma[tid]);
        mb_aligned_free(m_mat_indices[tid]);
    }

    mb_aligned_free(m_result_scratch);
    delete[] m_source_data;

    delete[] m_d;
    delete[] m_e;
    delete[] m_h;
    delete[] m_p;
    delete[] m_fac_a;
    delete[] m_fac_b;
    delete[] m_fac_c;
    delete[] m_gamma;
    delete[] m_mat_indices;
}

template<unsigned int num_lvl, template<unsigned int> class density_algo>
void
solver_cpu_fdtd_red<num_lvl, density_algo>::update_e(
    uint64_t size,
    unsigned int border,
    real* e,
    real* h,
    real* p,
    real* fac_a,
    real* fac_b,
    real* gamma) const
{
#if USE_OMP_SIMD
#pragma omp simd aligned(e, h, p, fac_a, fac_b, gamma : ALIGN)
#endif
    for (int i = border; i < size - border - 1; i++) {
        e[i] = fac_a[i] * e[i] +
            fac_b[i] * (-gamma[i] * p[i] + m_dx_inv * (h[i + 1] - h[i]));
    }
}

template<unsigned int num_lvl, template<unsigned int> class density_algo>
void
solver_cpu_fdtd_red<num_lvl, density_algo>::update_h(
    uint64_t size,
    unsigned int border,
    real* e,
    real* h,
    real* fac_c) const
{
#if USE_OMP_SIMD
#pragma omp simd aligned(e, h, fac_c : ALIGN)
#endif
    for (int i = border + 1; i < size - border; i++) {
        h[i] += fac_c[i] * (e[i] - e[i - 1]);
    }
}

template<unsigned int num_lvl, template<unsigned int> class density_algo>
void
solver_cpu_fdtd_red<num_lvl, density_algo>::apply_sources(
    real* t_e,
    real* source_data,
    unsigned int num_sources,
    uint64_t time,
    unsigned int base_pos,
    uint64_t chunk) const
{
    for (unsigned int k = 0; k < num_sources; k++) {
        int at = m_sim_sources[k].x_idx - base_pos + OL;
        if ((at > 0) && (at < chunk + 2 * OL)) {
            real src = source_data[m_sim_sources[k].data_base_idx + time];
            if (m_sim_sources[k].type == source::type::hard_source) {
                t_e[at] = src;
            } else if (m_sim_sources[k].type == source::type::soft_source) {
                /* TODO: fix source */
                t_e[at] += src;
            } else {
            }
        }
    }
}

template<unsigned int num_lvl, template<unsigned int> class density_algo>
void
solver_cpu_fdtd_red<num_lvl, density_algo>::run() const
{
    unsigned int P = omp_get_max_threads();
    uint64_t num_gridpoints = m_scenario->get_num_gridpoints();
    uint64_t chunk_base = m_scenario->get_num_gridpoints() / P;
    uint64_t chunk_rem = m_scenario->get_num_gridpoints() % P;
    uint64_t num_timesteps = m_scenario->get_num_timesteps();
    unsigned int num_sources = m_sim_sources.size();
    unsigned int num_copy = m_copy_list.size();

#pragma omp parallel
    {
        unsigned int tid = omp_get_thread_num();
        uint64_t chunk = chunk_base;
        if (tid == P - 1) {
            chunk += chunk_rem;
        }
        uint64_t size = chunk + 2 * OL;

        /* gather pointers */
        typedef typename density_algo<num_lvl>::density density_t;
        density_t* t_d;
        real *t_h, *t_e, *t_p;
        real *t_fac_a, *t_fac_b, *t_fac_c, *t_gamma;
        unsigned int* t_mat_indices;

        t_d = m_d[tid];
        t_e = m_e[tid];
        t_h = m_h[tid];
        t_p = m_p[tid];
        t_fac_a = m_fac_a[tid];
        t_fac_b = m_fac_b[tid];
        t_fac_c = m_fac_c[tid];
        t_gamma = m_gamma[tid];
        t_mat_indices = m_mat_indices[tid];

        __mb_assume_aligned(t_d);
        __mb_assume_aligned(t_e);
        __mb_assume_aligned(t_h);
        __mb_assume_aligned(t_p);
        __mb_assume_aligned(t_fac_a);
        __mb_assume_aligned(t_fac_b);
        __mb_assume_aligned(t_fac_c);
        __mb_assume_aligned(t_gamma);
        __mb_assume_aligned(t_mat_indices);
        __mb_assume_aligned(m_result_scratch);

        /* gather prev and next pointers from other threads */
        density_t *n_d, *p_d;
        real *n_h, *n_e;
        real *p_h, *p_e;

        __mb_assume_aligned(p_d);
        __mb_assume_aligned(p_e);
        __mb_assume_aligned(p_h);

        __mb_assume_aligned(n_d);
        __mb_assume_aligned(n_e);
        __mb_assume_aligned(n_h);

        if (tid > 0) {
            p_d = m_d[tid - 1];
            p_e = m_e[tid - 1];
            p_h = m_h[tid - 1];
        }

        if (tid < P - 1) {
            n_d = m_d[tid + 1];
            n_e = m_e[tid + 1];
            n_h = m_h[tid + 1];
        }

        /* main loop */
        for (uint64_t n = 0; n <= num_timesteps / OL; n++) {
            /* handle loop remainder */
            unsigned int subloop_ct =
                (n == num_timesteps / OL) ? num_timesteps % OL : OL;

            /* exchange data */
            if (tid > 0) {
#pragma ivdep
                for (unsigned int i = 0; i < OL; i++) {
                    t_d[i] = p_d[chunk_base + i];
                    t_e[i] = p_e[chunk_base + i];
                    t_h[i] = p_h[chunk_base + i];
                }
            }

            if (tid < P - 1) {
#pragma ivdep
                for (unsigned int i = 0; i < OL; i++) {
                    t_d[OL + chunk_base + i] = n_d[OL + i];
                    t_e[OL + chunk_base + i] = n_e[OL + i];
                    t_h[OL + chunk_base + i] = n_h[OL + i];
                }
            }

            /* sync after communication */
#pragma omp barrier

            /* sub-loop */
            for (unsigned int m = 0; m < subloop_ct; m++) {
                /* align border to vector length */
                unsigned int border = m - (m % VEC);

                for (int i = m; i < size - m - 1; i++) {
                    int mat_idx = t_mat_indices[i];

                    /* update density matrix */
                    density_algo<num_lvl>::update(
                        m_sim_consts_qm[mat_idx], t_d[i], t_e[i], &t_p[i]);
                }

                /* update electric field */
                update_e(
                    size, border, t_e, t_h, t_p, t_fac_a, t_fac_b, t_gamma);

                /* apply sources */
                apply_sources(
                    t_e,
                    m_source_data,
                    num_sources,
                    n * OL + m,
                    tid * chunk_base,
                    chunk);

                /* update magnetic field */
                update_h(size, border, t_e, t_h, t_fac_c);

                /* apply field boundary condition */
                if (tid == 0) {
                    t_h[OL] = 0;
                }
                if (tid == P - 1) {
                    t_h[OL + chunk] = 0;
                }

                /* save results to scratchpad in parallel */
                for (int k = 0; k < num_copy; k++) {
                    if (m_copy_list[k].hasto_record(n * OL + m)) {
                        uint64_t pos = m_copy_list[k].get_position();
                        uint64_t cols = m_copy_list[k].get_cols();
                        uint64_t ridx = m_copy_list[k].get_row_idx();
                        uint64_t cidx = m_copy_list[k].get_col_idx();
                        record::type t = m_copy_list[k].get_type();

                        int64_t base_idx = tid * chunk_base - OL;
                        int64_t off_r =
                            m_copy_list[k].get_offset_scratch_real(
                                n * OL + m, base_idx - pos);

                        /* TODO switch instead if else
                         * exchange type switch and for loop
                         * try again with vectorize
                         */

                        //#pragma omp simd
                        for (uint64_t i = OL; i < chunk + OL; i++) {
                            int64_t idx = base_idx + i;
                            if ((idx >= pos) && (idx < pos + cols)) {
                                if (t == record::type::electric) {
                                    m_result_scratch[off_r + i] = t_e[i];
                                } else if (t == record::type::magnetic) {
                                    m_result_scratch[off_r + i] = t_h[i];
                                } else if (t == record::type::inversion) {
                                    m_result_scratch[off_r + i] =
                                        density_algo<num_lvl>::calc_inversion(
                                            t_d[i]);

                                } else if (t == record::type::density) {

                                    /* right now only populations */
                                    if (ridx == cidx) {
                                        m_result_scratch[off_r + i] =
                                            density_algo<num_lvl>::
                                                calc_population(t_d[i], ridx);
                                    } else {
                                        /* coherence terms */

                                        /* TODO */
                                    }

                                    /* TODO: coherences
                                     * remove 1/3
                                     * consider only two/one corresponding
                                     * entry */

                                } else {
                                    /* TODO handle trouble */
                                    /* TODO calculate regular populations */
                                }
                            }
                        }
                    }
                }
            } /* end sub loop */

            /* sync after computation */
#pragma omp barrier
        } /* end main foor loop */

    } /* end openmp region */

    /* bulk copy results into result classes */
    for (const auto& cle : m_copy_list) {
        real* dr = m_result_scratch + cle.get_offset_scratch_real(0, 0);
        std::copy(dr, dr + cle.get_size(), cle.get_result_real(0, 0));
        if (cle.is_complex()) {
            real* di = m_result_scratch + cle.get_offset_scratch_imag(0, 0);
            std::copy(di, di + cle.get_size(), cle.get_result_imag(0, 0));
        }
    }
}
}

#endif
