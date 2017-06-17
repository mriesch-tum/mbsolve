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

namespace mbsolve {

static solver_factory<solver_openmp_2lvl_pc> factory("openmp-2lvl-pc");

/* TODO: necessary? */
unsigned int num_grid_points;
unsigned int num_time_steps;
real time_step_size;
#if 0
    num_grid_points = m_scenario.NumGridPoints;
    num_time_steps = m_scenario.NumTimeSteps;
    time_step_size = m_scenario.TimeStepSize;
#endif


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

    /* set up simulaton constants */
    std::map<std::string, unsigned int> id_to_idx;
    m_sim_consts = init_sim_constants(dev, scen, id_to_idx);

    /* allocate data arrays */
    m_inv = new real[scen->get_num_gridpoints()];
    m_dm12r = new real[scen->get_num_gridpoints()];
    m_dm12i = new real[scen->get_num_gridpoints()];
    m_h = new real[scen->get_num_gridpoints() + 1];
    m_e = new real[scen->get_num_gridpoints()];
    m_mat_indices = new unsigned int[scen->get_num_gridpoints()];

    /* set up indices array and initialize data arrays */
#pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < scen->get_num_gridpoints(); i++) {
        /* determine index of material */
        int idx = -1;
        real x = i * scen->get_gridpoint_size();
        for (const auto& reg : dev->get_regions()) {
            if ((x >= reg->get_start()) && (x <= reg->get_end())) {
                idx = id_to_idx[reg->get_material()->get_id()];
                break;
            }
        }
        /* TODO: assert/bug if idx == -1 */
        if ((idx < 0) || (idx >= dev->get_used_materials().size())) {
            std::cout << "At index " << i << std::endl;
            throw std::invalid_argument("region not found");
        }
        m_mat_indices[i] = idx;

        /* TODO: evaluate flexible initialization in scenario */
        m_inv[i] = m_sim_consts[idx].inversion_init;
        m_dm12r[i] = 0.0;
        m_dm12i[i] = 0.0;
        m_e[i] = 0.0;
        m_h[i] = 0.0;
        if (i == scen->get_num_gridpoints() - 1) {
            m_h[i + 1] = 0.0;
        }
    }

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

solver_openmp_2lvl_pc::~solver_openmp_2lvl_pc()
{
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
solver_openmp_2lvl_pc::get_name() const
{
    return factory.get_name();
}

void
solver_openmp_2lvl_pc::run() const
{
    /*
#pragma offload target(mic) in(gsc, num_grid_points, num_time_steps) \
  in(time_step_size) \
  inout(m_e:length(m_scenario.NumGridPoints)) \
  inout(m_h:length(m_scenario.NumGridPoints + 1)) \
  inout(m_dm11:length(m_scenario.NumGridPoints)) \
  inout(m_dm12r:length(m_scenario.NumGridPoints)) \
  inout(m_dm12i:length(m_scenario.NumGridPoints)) \
  inout(m_dm22:length(m_scenario.NumGridPoints)) \
  in(region_indices:length(m_scenario.NumGridPoints))
  {
    */
#pragma omp parallel
    {
          /* main loop */
          for (int n = 0; n < m_scenario->get_num_timesteps(); n++) {
              /* update dm and e in parallel */
              //#pragma omp for simd schedule(static)
#pragma omp for schedule(static)
              for (int i = 0; i < m_scenario->get_num_gridpoints(); i++) {
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
                      real OmRabi = m_sim_consts[mat_idx].d12 * e;

                      inv_e = m_inv[i] + m_sim_consts[mat_idx].d_t *
                          (- 4.0 * OmRabi * rho12i
                           - m_sim_consts[mat_idx].tau1 *
                           (inv - m_sim_consts[mat_idx].equi_inv));

                      rho12i_e = m_dm12i[i] + m_sim_consts[mat_idx].d_t *
                          (- m_sim_consts[mat_idx].w12 * rho12r
                           + OmRabi * inv
                           - m_sim_consts[mat_idx].gamma12 * rho12i);

                      rho12r_e = m_dm12r[i] + m_sim_consts[mat_idx].d_t *
                          (+ m_sim_consts[mat_idx].w12 * rho12i
                           - m_sim_consts[mat_idx].gamma12 * rho12r);

                      real j = m_sim_consts[mat_idx].sigma * e;

                      real p_t = m_sim_consts[mat_idx].M_CP
                          * m_sim_consts[mat_idx].d12 *
                          (m_sim_consts[mat_idx].w12 * rho12i -
                           m_sim_consts[mat_idx].gamma12 * rho12r);

                      field_e = m_e[i] + m_sim_consts[mat_idx].M_CE *
                          (-j - p_t + (m_h[i + 1] - m_h[i])
                           * m_sim_consts[mat_idx].d_x_inv);
                  }

                  /* final update step */
                  m_inv[i] = inv_e;
                  m_dm12i[i] = rho12i_e;
                  m_dm12r[i] = rho12r_e;

                  m_e[i] = field_e;
              }

              /* apply sources */
              for (const auto& src : m_sim_sources) {
                  /* TODO: support other source types than hard sources */
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
              for (int i = 1; i < m_scenario->get_num_gridpoints(); i++) {
                  unsigned int mat_idx = m_mat_indices[i - 1];

                  m_h[i] += m_sim_consts[mat_idx].M_CH * (m_e[i] - m_e[i - 1]);
              }

              /* save results to scratchpad in parallel */
              for (const auto& cle : m_copy_list) {
                  if (cle.hasto_record(n)) {
#pragma omp for schedule(static)
                      for (int i = 0; i < m_scenario->get_num_gridpoints();
                           i++) {
                          unsigned int pos = cle.get_position();
                          if ((i >= pos) && (i < pos + cle.get_cols())) {
                              *cle.get_scratch_real(n, i - pos) =
                                  *cle.get_real(i);
                              if (cle.is_complex()) {
                                  *cle.get_scratch_imag(n, i - pos) =
                                      *cle.get_imag(i);
                              }
                          }
                      }
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
