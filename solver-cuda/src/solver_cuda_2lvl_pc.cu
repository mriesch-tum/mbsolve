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

#include <curand.h>
#include <curand_kernel.h>
#include <common_cuda.hpp>
#include <solver_cuda_2lvl_pc.hpp>

namespace mbsolve {

static solver_factory<solver_cuda_2lvl_pc> factory("cuda-2lvl-pc");

/* initialize memory kernel */
__global__ void init_memory(real *d, real *e, real *h, unsigned int *indices)
{
    unsigned int gsize = blockDim.x * gridDim.x;
    unsigned int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int mat_idx = indices[gidx];

    /* TODO: alternative initializations */
    curandState_t rand_state;

    /* initialize random number generator */
    curand_init(clock64(), gidx, 0, &rand_state);

    if (gidx == blockDim.x * gridDim.x - 1) {
	h[gidx + 1] = 0.0;
    }
    h[gidx] = 0.0;
    e[gidx] = 0.0;
    //   e[gidx] = curand_uniform(&rand_state) * 1e-15;

    d[gsize * 0 + gidx] = l_sim_consts[mat_idx].inversion_init;
    d[gsize * 1 + gidx] = 0.0;
    d[gsize * 2 + gidx] = 0.0;
}

__global__ void makestep_h(const real *ge, real *gh, unsigned int *indices)
{
    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int mat_idx = indices[gidx];

    extern __shared__ real e[];

    if ((idx == 0) && (gidx != 0)) {
	e[idx] = ge[gidx - 1];
    }
    e[idx + 1] = ge[gidx];

    __syncthreads();

    /* TODO: alternative boundary conditions? */
    /* TODO: different kernel or templated version?? */
    /* open circuit boundary conditions already set */
    /* gh_ghz[0] = 0; */
    /* gh_ghz[N_x] = 0; */

    if (gidx != 0) {
	gh[gidx] += l_sim_consts[mat_idx].M_CH * (e[idx + 1] - e[idx]);
    }
}


__global__ void makestep_e_dm(real *d, const real *gh, real *ge,
			      unsigned int *indices, real *scratch, real *sd,
                              unsigned int source_ct, unsigned int copy_ct,
                              unsigned int n)
{
    int gsize = blockDim.x * gridDim.x;
    int size = blockDim.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int idx = threadIdx.x;
    int mat_idx = indices[gidx];

    extern __shared__ real h[];
    real *inv = &h[size + 1];
    real *dm12i = &h[2 * size + 1];
    real *dm12r = &h[3 * size + 1];
    real *e = &h[5 * size + 1];

    h[idx] = gh[gidx];
    if (idx == blockDim.x - 1) {
	h[idx + 1] = gh[gidx + 1];
    }
    inv[idx] = d[gsize * 0 + gidx];
    dm12i[idx] = d[gsize * 1 + gidx];
    dm12r[idx] = d[gsize * 2 + gidx];
    e[idx] = ge[gidx];

    __syncthreads();

    real inv_e = inv[idx];
    real dm12i_e = dm12i[idx];
    real dm12r_e = dm12r[idx];
    real e_e = e[idx];

    /* execute prediction - correction steps */
    for (int pc_step = 0; pc_step < 4; pc_step++) {

	real inversion = 0.5 * (inv[idx] + inv_e);
	real rho12i = 0.5 * (dm12i[idx] + dm12i_e);
	real rho12r = 0.5 * (dm12r[idx] + dm12r_e);
        real field = 0.5 * (e[idx] + e_e);
	real OmRabi = l_sim_consts[mat_idx].d12 * field;

	/* dm11 */
	inv_e = inv[idx] + l_sim_consts[mat_idx].d_t *
	    (- 4.0 * OmRabi * rho12i
             - l_sim_consts[mat_idx].tau1 *
             (inversion - l_sim_consts[mat_idx].equi_inv));

	/* imag dm12 */
	dm12i_e = dm12i[idx] + l_sim_consts[mat_idx].d_t *
	    (- l_sim_consts[mat_idx].w12 * rho12r
             + OmRabi * inversion
	     - l_sim_consts[mat_idx].gamma12 * rho12i);

	/* real dm12 */
	dm12r_e = dm12r[idx] + l_sim_consts[mat_idx].d_t *
	    (+ l_sim_consts[mat_idx].w12 * rho12i
             - l_sim_consts[mat_idx].gamma12 * rho12r);

	real j = l_sim_consts[mat_idx].sigma * field;

	real p_t = l_sim_consts[mat_idx].M_CP * l_sim_consts[mat_idx].d12 *
	    (+ l_sim_consts[mat_idx].w12 * rho12i
             - l_sim_consts[mat_idx].gamma12 * rho12r);

	e_e = e[idx] + l_sim_consts[mat_idx].M_CE *
	    (-j - p_t + (h[idx + 1] - h[idx]) * l_sim_consts[mat_idx].d_x_inv);
    }

    /* apply sources */
    for (unsigned int k = 0; k < source_ct; k++) {
        int at = l_sim_sources[k].x_idx;
	if (gidx == at) {
            if (l_sim_sources[k].type == source::type::hard_source) {
                e_e = sd[l_sim_sources[k].data_base_idx + n];
            } else if (l_sim_sources[k].type == source::type::soft_source) {
                e_e += sd[l_sim_sources[k].data_base_idx + n];
            }
	}
    }

    ge[gidx] = e_e;
    d[gsize * 0 + gidx] = inv_e;
    d[gsize * 1 + gidx] = dm12i_e;
    d[gsize * 2 + gidx] = dm12r_e;

    /* copy results into scratchpad memory */
    for (unsigned int k = 0; k < copy_ct; k++) {
        if (l_copy_list[k].hasto_record(n)) {
            unsigned int pos = l_copy_list[k].get_position();
            record::type t = l_copy_list[k].get_type();
            real src_real;

            if (t == record::type::electric) {
                src_real = e_e;
            } else if (t == record::type::inversion) {
                src_real = inv_e;
            } else {
                /* TODO handle trouble, handle complex quantities */
            }

            if ((gidx >= pos) && (gidx < pos + l_copy_list[k].get_cols())) {
                int off_r = l_copy_list[k].get_offset_scratch_real(n,
                                                                   gidx - pos);
                scratch[off_r] = src_real;
            }
        }
    }
}

/* host members */
solver_cuda_2lvl_pc::solver_cuda_2lvl_pc(std::shared_ptr<const device> dev,
                                         std::shared_ptr<scenario> scen) :
solver_int(dev, scen)
{
    /* determine simulation settings */
    init_fdtd_simulation(dev, scen, 0.5);

    /* set up simulaton constants */
    std::map<std::string, unsigned int> id_to_idx;
    m_sim_consts = init_sim_constants(dev, scen, id_to_idx);

    /* TODO: sanity check scenario? */

    /* TODO: handle gridpoint number with %32 != 0 ? */
    if (scen->get_num_gridpoints() % 32 != 0) {
	throw std::invalid_argument("Number of grid points must be multiple"
				    " of 32");
    }

    if (dev->get_regions().size() == 0) {
        throw std::invalid_argument("No regions in device!");
    }

    if (m_sim_consts.size() >= MB_CUDA_MAX_MATERIALS) {
        throw std::invalid_argument("Too many materials in device!");
    }

    if (scen->get_sources().size() >= MB_CUDA_MAX_SOURCES) {
        throw std::invalid_argument("Too many sources in scenario!");
    }

    if (scen->get_records().size() >= MB_CUDA_MAX_CLE) {
        throw std::invalid_argument("Too many records requested in scenario!");
    }

    /* determine indices */
    unsigned int *mat_indices = new unsigned int[scen->get_num_gridpoints()];
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
            throw std::invalid_argument("region not found");
        }
        mat_indices[i] = idx;
    }

    /* copy settings to CUDA constant memory */
    chk_err(cudaMemcpyToSymbol(l_sim_consts, m_sim_consts.data(),
                               m_sim_consts.size() *
                               sizeof(sim_constants_2lvl)));

    /* allocate buffers on GPU */
    chk_err(cudaMalloc(&m_e, sizeof(real) * scen->get_num_gridpoints()));
    chk_err(cudaMalloc(&m_h, sizeof(real) * (scen->get_num_gridpoints() + 1)));
    chk_err(cudaMalloc(&m_d, sizeof(real) * scen->get_num_gridpoints() * 3));
    chk_err(cudaMalloc(&m_mat_indices, sizeof(unsigned int) *
                       scen->get_num_gridpoints()));

    /* copy indices to GPU */
    chk_err(cudaMemcpy(m_mat_indices, mat_indices, sizeof(unsigned int) *
                       scen->get_num_gridpoints(), cudaMemcpyHostToDevice));

    delete[] mat_indices;

    /* initialize memory */
    /* TODO move to class member */
    unsigned int threads = 128;
    unsigned int blocks = scen->get_num_gridpoints()/threads;

    init_memory<<<blocks, threads>>>(m_d, m_e, m_h, m_mat_indices);

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
    chk_err(cudaMalloc(&m_result_scratch, sizeof(real) * scratch_size));

    /* set up sources */
    unsigned int source_data_size = scen->get_num_timesteps() *
        scen->get_sources().size();
    chk_err(cudaMalloc(&m_source_data, sizeof(real) * source_data_size));
    real *source_data = new real[source_data_size];
    unsigned int base_idx = 0;
    for (const auto& src : scen->get_sources()) {
        sim_source s;
        s.type = src->get_type();
        s.x_idx = src->get_position()/scen->get_gridpoint_size();
        s.data_base_idx = base_idx;
        m_sim_sources.push_back(s);

        /* calculate source values */
        for (unsigned int j = 0; j < scen->get_num_timesteps(); j++) {
            source_data[base_idx + j] =
                src->get_value(j * scen->get_timestep_size());
        }

        base_idx += scen->get_num_timesteps();
    }

    /* copy indices to GPU */
    chk_err(cudaMemcpy(m_source_data, source_data, sizeof(unsigned int) *
                       source_data_size, cudaMemcpyHostToDevice));
    delete[] source_data;

    /* copy source properties to GPU constant memory */
    chk_err(cudaMemcpyToSymbol(l_sim_sources, m_sim_sources.data(),
                               m_sim_sources.size() *
                               sizeof(sim_source)));

    /* copy copy list entries to GPU constant memory */
    for (unsigned int i = 0; i < m_copy_list.size(); i++) {
        chk_err(cudaMemcpyToSymbol(l_copy_list, &m_copy_list[i].get_dev(),
                                   sizeof(copy_list_entry_dev),
                                   sizeof(copy_list_entry_dev) * i));

    }

}

solver_cuda_2lvl_pc::~solver_cuda_2lvl_pc()
{
    /* free CUDA memory */
    cudaFree(m_h);
    cudaFree(m_e);
    cudaFree(m_d);
    cudaFree(m_mat_indices);
    cudaFree(m_result_scratch);
    cudaFree(m_source_data);

    /* reset device */
    cudaDeviceReset();
}

const std::string&
solver_cuda_2lvl_pc::get_name() const
{
    return factory.get_name();
}

void
solver_cuda_2lvl_pc::run() const
{
    unsigned int threads = 128;
    unsigned int blocks = m_scenario->get_num_gridpoints()/threads;
    /* TODO handle roundoff errors in thread/block partition */

    /* main loop */
    for (unsigned int i = 0; i < m_scenario->get_num_timesteps(); i++) {

	/* makestep e and density matrix */
	makestep_e_dm<<<blocks, threads,
	    (6 * threads + 1) * sizeof(real)>>>(m_d, m_h, m_e, m_mat_indices,
                                                m_result_scratch,
                                                m_source_data,
                                                m_sim_sources.size(),
                                                m_copy_list.size(),
                                                i);

	/* makestep h */
	makestep_h<<<blocks, threads, (threads + 1) * sizeof(real)>>>
	    (m_e, m_h, m_mat_indices);

    }

    /* bulk copy results into result classes */
    for (const auto& cle : m_copy_list) {
        chk_err(cudaMemcpy(cle.get_result()->get_data_real_raw(),
                           m_result_scratch +
                           cle.get_offset_scratch_real(0, 0),
                           cle.get_size() * sizeof(real),
                           cudaMemcpyDeviceToHost));
        if (cle.is_complex()) {
            /*
            chk_err(cudaMemcpy(cle.get_result_imag(0, 0).data(),
                               cle.get_scratch_imag(0, 0), cle.get_size(),
                               cudaMemcpyDeviceToHost));
            */
        }
    }
}

}
