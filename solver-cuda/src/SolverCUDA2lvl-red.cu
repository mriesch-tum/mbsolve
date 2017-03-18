#include <boost/foreach.hpp>

#include <curand.h>
#include <curand_kernel.h>
#include <CUDACommon.hpp>
#include <SolverCUDA2lvl-red.hpp>

namespace mbsolve {

static SolverFactory<SolverCUDA2lvl_red> factory("cuda-2lvl-red");

/* CUDA memory and kernels */
__device__ __constant__ struct sim_constants gsc_red[MaxRegions];

__device__ __constant__ CopyListEntry copy_list[10];

static const unsigned int threads = 256;

/* TODO: initialize may be reused by other CUDA solvers, make general */
/* TODO: region-wise? initialization */
__global__ void init_memory_red(real *d, real *e, real *h,
				unsigned int *indices)
{
    unsigned int gsize = blockDim.x * gridDim.x;
    unsigned int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = 0;

    /* determine region index */
    for (unsigned int i = 0; i < MaxRegions; i++) {
	if (gidx <= gsc_red[i].idx_end) {
	    region = i;
	    break;
	}
    }

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

    d[gsize * 0 + gidx] = gsc_red[region].dm11_init;
    d[gsize * 1 + gidx] = 0.0;
    d[gsize * 2 + gidx] = 0.0;
    d[gsize * 3 + gidx] = gsc_red[region].dm22_init;

    indices[gidx] = region;
}

const unsigned int OL = 32;

__global__ void makestep(real *d, real *gh, real *ge, unsigned int n,
			     unsigned int *indices, unsigned int copy_list_ct,
			     real *src, real *d_new, real *gh_new,
			     real *ge_new)
{
    int gsize = threads * gridDim.x;
    int gidx = threads * blockIdx.x + threadIdx.x - OL;
    int idx = threadIdx.x;
    int region = 0;

    __shared__ real dm11[threads + 2 * OL];
    __shared__ real dm12i[threads + 2 * OL];
    __shared__ real dm12r[threads + 2 * OL];
    __shared__ real dm22[threads + 2 * OL];
    __shared__ real h[threads + 2 * OL];
    __shared__ real e[threads + 2 * OL];

    if ((gidx >= 0) && (gidx < gsize)) {
	h[idx] = gh[gidx];
	e[idx] = ge[gidx];
	dm11[idx] = d[gsize * 0 + gidx];
	dm12i[idx] = d[gsize * 1 + gidx];
	dm12r[idx] = d[gsize * 2 + gidx];
	dm22[idx] = d[gsize * 3 + gidx];

	region = indices[gidx];
    }
    if (gidx == gsize) {
	h[idx] = 0;
    }

    for (int i = 0; i < OL; i++) {
	__syncthreads();

	real dm11_e = dm11[idx];
	real dm12i_e = dm12i[idx];
	real dm12r_e = dm12r[idx];
	real dm22_e = dm22[idx];
	real e_e = e[idx];

	if ((idx >= i) && (idx < blockDim.x - i - 1)) {
	    /* execute prediction - correction steps */
#pragma unroll
	    for (int pc_step = 0; pc_step < 4; pc_step++) {
		real rho11 = 0.5 * (dm11[idx] + dm11_e);
		real rho12i = 0.5 * (dm12i[idx] + dm12i_e);
		real rho12r = 0.5 * (dm12r[idx] + dm12r_e);
		real rho22 = 0.5 * (dm22[idx] + dm22_e);
		real OmRabi = gsc_red[region].d12 * 0.5 * (e[idx] + e_e);

		/* dm11 */
		dm11_e = dm11[idx] + gsc_red[region].d_t *
		    (-2.0 * OmRabi * rho12i - gsc_red[region].tau1 * rho11);

		/* imag dm12 */
		dm12i_e = dm12i[idx] + gsc_red[region].d_t *
		    (- gsc_red[region].w12 * rho12r + OmRabi * (rho11 - rho22)
		     - gsc_red[region].gamma12 * rho12i);

		/* real dm12 */
		dm12r_e = dm12r[idx] + gsc_red[region].d_t *
		    (+ gsc_red[region].w12 * rho12i
		     - gsc_red[region].gamma12 * rho12r);

		/* dm22 */
		dm22_e = dm22[idx] + gsc_red[region].d_t *
		    (2.0 * OmRabi * rho12i + gsc_red[region].tau1 * rho11);

		real j = 0; /* ge[gidx] * gsc[region].sigma;*/

		real p_t = gsc_red[region].M_CP * gsc_red[region].d12 *
		    (+ gsc_red[region].w12 * rho12i
		     - gsc_red[region].gamma12 * rho12r);

		e_e = e[idx] + gsc_red[region].M_CE *
		    (-j - p_t +
		     (h[idx + 1] - h[idx]) * gsc_red[region].d_x_inv);

		if (gidx == 0) {
		    e_e = src[OL * n + i]; /* hard source */
		}
	    }
	    e[idx] = e_e;
	    dm11[idx] = dm11_e;
	    dm12i[idx] = dm12i_e;
	    dm12r[idx] = dm12r_e;
	    dm22[idx] = dm22_e;
	}

	__syncthreads();

	if ((idx > i) && (idx < blockDim.x - i - 1)) {
	    h[idx] += gsc_red[region].M_CH * (e[idx] - e[idx - 1]);
	}
	if (gidx == 0) {
	    h[idx] = 0;
	}
	if (gidx == gsize) {
	    h[idx] = 0;
	}

	/* copy result data */
	for (int k = 0; k < copy_list_ct; k++) {
	    if (copy_list[k].record(n * OL + i)) {
		if ((idx >= OL) && (idx < OL + threads)) {
		    //if ((gidx >= copy_list[k].get_position()) &&
		    //	(gidx < copy_list[k].get_position()
		    //	 + copy_list[k].get_count())) {
			switch (copy_list[k].get_type()) {
			case EField:
			    copy_list[k].get_dst(n * OL + i)[gidx] = e[idx];
			    break;
			case D11:
			    copy_list[k].get_dst(n * OL + i)[gidx] = dm11[idx];
			    break;
			case D22:
			    copy_list[k].get_dst(n * OL + i)[gidx] = dm22[idx];
			    break;
			default:
			    break;
			}
			//}
		}
	    }
	}
    }

    if ((idx >= OL) && (idx < OL + threads)) {
	gh_new[gidx] = h[idx];
	ge_new[gidx] = e[idx];
	d_new[gsize * 0 + gidx] = dm11[idx];
	d_new[gsize * 1 + gidx] = dm12i[idx];
	d_new[gsize * 2 + gidx] = dm12r[idx];
	d_new[gsize * 3 + gidx] = dm22[idx];
	/*
	gh[gidx] = h[idx];
	ge[gidx] = e[idx];
	d[gsize * 0 + gidx] = dm11[idx];
	d[gsize * 1 + gidx] = dm12i[idx];
	d[gsize * 2 + gidx] = dm12r[idx];
	d[gsize * 3 + gidx] = dm22[idx];
	*/
    }
}

/* host members */
SolverCUDA2lvl_red::SolverCUDA2lvl_red(const Device& device,
				       const Scenario& scenario) :
    ISolver(device, scenario), comp_maxwell(0), copy(0)
{
    /* total device length */
    Quantity length = device.XDim();

    /* minimum relative permittivity */
    Quantity minRelPermittivity = device.MinRelPermittivity();

    /* TODO: sanity check scenario? */
    /*    if (m_scenario.NumGridPoints % 32 != 0) {
	throw std::invalid_argument("Number of grid points must be multiple"
				    " of 32");
				    }*/

    /* determine grid point and time step size */
    real C = 0.5; /* courant number */
    real velocity_inv = sqrt(MU0() * EPS0() * minRelPermittivity());
    m_scenario.GridPointSize = length()/(m_scenario.NumGridPoints - 1);
    real timestep  = C * m_scenario.GridPointSize * velocity_inv;
    m_scenario.NumTimeSteps = ceil(m_scenario.SimEndTime/timestep) + 1;
    m_scenario.TimeStepSize = m_scenario.SimEndTime /
	(m_scenario.NumTimeSteps - 1);


    /* determine border indices and initialize region settings */
    if (device.Regions.size() > MaxRegions) {
	throw std::invalid_argument("Too many regions requested");
    }
    struct sim_constants sc[MaxRegions];

    unsigned int i = 0;
    BOOST_FOREACH(Region reg, device.Regions) {
	if (i > 0) {
	    sc[i - 1].idx_end = round(reg.X0()/m_scenario.GridPointSize) - 1;
	}
	sc[i].idx_start = round(reg.X0()/m_scenario.GridPointSize);
	sc[i].M_CE = m_scenario.TimeStepSize/(EPS0() * reg.RelPermittivity());
	sc[i].M_CH = m_scenario.TimeStepSize/(MU0() *
					      m_scenario.GridPointSize);
	/* TODO: overlap factor? */
	sc[i].M_CP = -2.0 * reg.DopingDensity * HBAR();
	sc[i].sigma = 2.0 * sqrt(EPS0 * reg.RelPermittivity/MU0) * reg.Losses;

	sc[i].w12 = (reg.TransitionFrequencies.size() < 1) ? 0.0 :
	    reg.TransitionFrequencies[0]();
	/* TODO rename to rabi freqs or something */
	sc[i].d12 = (reg.DipoleMoments.size() < 1) ? 0.0 :
	    reg.DipoleMoments[0]() * E0()/HBAR();
	sc[i].tau1 = (reg.ScatteringRates.size() < 1) ? 0.0 :
	    reg.ScatteringRates[0]();
	sc[i].gamma12 = (reg.DephasingRates.size() < 1) ? 0.0 :
	    reg.DephasingRates[0]();

	sc[i].d_x_inv = 1.0/m_scenario.GridPointSize;
	sc[i].d_t = m_scenario.TimeStepSize;

	if (reg.DopingDensity() < 1.0) {
	    sc[i].dm11_init = 0.0;
	    sc[i].dm22_init = 0.0;
	} else {
	    sc[i].dm11_init = 0.0;
	    sc[i].dm22_init = 1.0;
	}

	i++;
    }
    if (i > 0) {
	sc[i - 1].idx_end = m_scenario.NumGridPoints - 1;
    }

    /* copy settings to CUDA constant memory */
    chk_err(cudaMemcpyToSymbol(gsc_red, &sc, MaxRegions *
			       sizeof(struct sim_constants)));

    /* initialize streams */
    chk_err(cudaStreamCreate(&comp_maxwell));
    chk_err(cudaStreamCreate(&copy));

    /* allocate space */
    chk_err(cudaMalloc(&m_e1, sizeof(real) * m_scenario.NumGridPoints));
    chk_err(cudaMalloc(&m_h1, sizeof(real) * (m_scenario.NumGridPoints + 1)));
    chk_err(cudaMalloc(&m_d1, sizeof(real) * m_scenario.NumGridPoints * 4));
    chk_err(cudaMalloc(&m_e2, sizeof(real) * m_scenario.NumGridPoints));
    chk_err(cudaMalloc(&m_h2, sizeof(real) * (m_scenario.NumGridPoints + 1)));
    chk_err(cudaMalloc(&m_d2, sizeof(real) * m_scenario.NumGridPoints * 4));
    chk_err(cudaMalloc(&m_indices, sizeof(unsigned int) *
		       m_scenario.NumGridPoints));
    chk_err(cudaMalloc(&m_src, sizeof(real) * m_scenario.NumTimeSteps));

    real *src_host = new real[m_scenario.NumTimeSteps];
    for (unsigned int k = 0; k < m_scenario.NumTimeSteps; k++) {
	/* calculate source value */
	real t = k * m_scenario.TimeStepSize;

	/* source 0 */
	real f_0 = 2e14;
	real T_p = 20/f_0;
	real E_0 = 4.2186e9;
	//E_0 /= 2; /* pi pulse */
	real gamma = 2 * t/T_p - 1;
	src_host[k] = E_0 * 1/std::cosh(10 * gamma) *
	    std::sin(2 * M_PI * f_0 * t);
    }

    chk_err(cudaMemcpy(m_src, src_host, m_scenario.NumTimeSteps * sizeof(real),
		       cudaMemcpyHostToDevice));
    delete[] src_host;

    /* set up results transfer data structures */
    CopyListEntry *list = new CopyListEntry[m_scenario.Records.size()];

    for (int k = 0; k < m_scenario.Records.size(); k++) {
	Record rec = m_scenario.Records[k];

	unsigned int row_ct = m_scenario.SimEndTime/rec.Interval;
	unsigned int interval = ceil(1.0 * m_scenario.NumTimeSteps/row_ct);
	unsigned int position_idx;
	unsigned int col_ct;

	if (rec.Position() < 0.0) {
	    /* copy complete grid */
	    position_idx = 0;
	    col_ct = m_scenario.NumGridPoints;
	} else {
	    position_idx = round(rec.Position()/m_scenario.GridPointSize);
	    col_ct = 1;
	}

	/* allocate result memory */
	Result *res = new Result(rec.Name, col_ct, row_ct);
	m_results.push_back(res);

	/* allocate results memory on GPU */
	real *data;
	chk_err(cudaMalloc(&data, sizeof(real) * col_ct * row_ct));
	m_results_gpu.push_back(data);

	/* create copy list entry */

	if (rec.Type == EField) {
	    list[k] = CopyListEntry(data, col_ct, position_idx, interval,
				    EField);
	} else if (rec.Type == HField) {
	    /* TODO: numGridPoints + 1 */
	    list[k] = CopyListEntry(data, col_ct, position_idx, interval,
				    HField);
	} else {
	    if ((rec.I - 1 < 2) && (rec.J - 1 < 2)) {
		if ((rec.I == 1) && (rec.J == 1)) {
		    list[k] = CopyListEntry(data, col_ct, position_idx,
					    interval, D11);
		} else if ((rec.I == 2) && (rec.J == 2)) {
		    list[k] = CopyListEntry(data, col_ct, position_idx,
					    interval, D22);
		} else {
		    /* off-diagonal entry */
		    /* TODO */
		    /* if complex */
		    /* create two list entries */
		    /* create two Results, or one complex Result */

		    /* real part: GetSrcDensity(&dm, rec.I, rec.J); */
		    /* imag part: GetSrcDensity(&dm, rec.J, rec.I); */
		}
	    } else {
		// throw exc
	    }
	}
    }

    /* copy copy_list entries to CUDA constant memory */
    chk_err(cudaMemcpyToSymbol(copy_list, list, m_scenario.Records.size() *
			       sizeof(CopyListEntry)));
    delete[] list;

    /* initialize memory */
    unsigned int blocks = m_scenario.NumGridPoints/threads;
    init_memory_red<<<blocks, threads>>>(m_d1, m_e1, m_h1, m_indices);

    /* sync */
    cudaDeviceSynchronize();
}

SolverCUDA2lvl_red::~SolverCUDA2lvl_red()
{
    /* sync */
    cudaDeviceSynchronize();

    /* delete result data on GPU */
    BOOST_FOREACH(real *data, m_results_gpu) {
	cudaFree(data);
    }

    /* free CUDA memory */
    cudaFree(m_h1);
    cudaFree(m_e1);
    cudaFree(m_d1);
    cudaFree(m_h2);
    cudaFree(m_e2);
    cudaFree(m_d2);

    cudaFree(m_indices);
    cudaFree(m_src);

    /* clean up streams */
    if (comp_maxwell) {
	cudaStreamDestroy(comp_maxwell);
    }
    if (copy) {
	cudaStreamDestroy(copy);
    }

    /* reset device */
    cudaDeviceReset();
}

std::string
SolverCUDA2lvl_red::getName() const
{
    return factory.getName();
}

void
SolverCUDA2lvl_red::run() const
{
    unsigned int blocks = m_scenario.NumGridPoints/threads;
    /* TODO handle roundoff errors in thread/block partition */

    /* main loop */
    for (unsigned int n = 0; n < m_scenario.NumTimeSteps/OL; n++) {

	/* makestep */
	if (n % 2 == 0) {
	    makestep<<<blocks, threads + 2 * OL>>>
		(m_d1, m_h1, m_e1, n, m_indices, m_scenario.Records.size(),
		 m_src, m_d2, m_h2, m_e2);
	} else {
	    makestep<<<blocks, threads + 2 * OL>>>
		(m_d2, m_h2, m_e2, n, m_indices, m_scenario.Records.size(),
		 m_src, m_d1, m_h1, m_e1);
	}
    }

    /* copy result data to host */
    for (int k = 0; k < m_results.size(); k++) {
	chk_err(cudaMemcpy(m_results[k]->data(), m_results_gpu[k],
			   m_results[k]->count() * sizeof(real),
			   cudaMemcpyDeviceToHost));
    }
}

}
