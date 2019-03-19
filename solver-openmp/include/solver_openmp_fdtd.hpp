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

#ifndef MBSOLVE_SOLVER_OPENMP_FDTD_H
#define MBSOLVE_SOLVER_OPENMP_FDTD_H

#include <solver.hpp>
#include <internal/common_fdtd.hpp>
#include <internal/copy_list_entry.hpp>

namespace mbsolve {

/**
 * OpenMP solver for c-lvl systems using the FDTD method.
 * The number of levels c can be chosen arbitrarily, but must be known at
 * compile time.
 * \ingroup MBSOLVE_SOLVER_OPENMP
 */
template<unsigned int num_lvl, template<unsigned int> class density_algo>
class solver_openmp_fdtd : public solver_int
{
public:
    solver_openmp_fdtd(std::shared_ptr<const device> dev,
                       std::shared_ptr<scenario> scen);

    ~solver_openmp_fdtd();

    const std::string& get_name() const;

    void run() const;

private:
    const std::string m_name;

    /* TODO: rule of three. make copy constructor etc. private?
     * or implement correctly
     */

    typename density_algo<num_lvl>::density *m_d;

    real *m_h;
    real *m_e;
    real *m_p;

    real *m_fac_a;
    real *m_fac_b;
    real *m_fac_c;
    real *m_gamma;

    real m_dx_inv;

    real *m_result_scratch;

    real *m_source_data;

    unsigned int *m_mat_indices;

    typedef typename density_algo<num_lvl>::sim_constants qm_consts;
    typedef typename density_algo<num_lvl>::allocator qm_allocator;
    std::vector<qm_consts, qm_allocator> m_sim_consts_qm;

    std::vector<sim_source> m_sim_sources;

    std::vector<copy_list_entry> m_copy_list;
};

}

#endif
