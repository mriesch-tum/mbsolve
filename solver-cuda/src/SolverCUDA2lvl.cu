#include <boost/foreach.hpp>

#include <curand.h>
#include <curand_kernel.h>
#include <CUDACommon.hpp>
#include <SolverCUDA2lvl.hpp>

namespace mbsolve {

static SolverFactory<SolverCUDA2lvl> factory("cuda-2lvl");

/* CUDA memory and kernels */
__device__ __constant__ struct sim_constants gsc[MaxRegions];

/* TODO: hash function or something? */
/* TODO: otherwise divergence within warp */
__device__ __inline__ unsigned int get_region(unsigned int idx)
{
    for (unsigned int i = 0; i < MaxRegions; i++) {
	if (idx <= gsc[i].idx_end) {
	    return i;
	}
    }
    return 0;
}

/* TODO: initialize may be reused by other CUDA solvers, make general */
/* TODO: region-wise? initialization */
__global__ void init_memory(real *d, real *e, real *h, unsigned int *indices)
{
    unsigned int gsize = blockDim.x * gridDim.x;
    unsigned int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);

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

    d[gsize * 0 + gidx] = gsc[region].dm11_init;
    d[gsize * 1 + gidx] = 0.0;
    d[gsize * 2 + gidx] = 0.0;
    d[gsize * 3 + gidx] = gsc[region].dm22_init;

    indices[gidx] = region;
}

__global__ void makestep_h(const real *ge, real *gh, unsigned int *indices)
{
    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = indices[gidx];

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
	gh[gidx] += gsc[region].M_CH * (e[idx + 1] - e[idx]);
    }
}

__global__ void makestep_e_dm(real *d, const real *gh, real *ge, real src,
			      unsigned int *indices)
{
    int gsize = blockDim.x * gridDim.x;
    int size = blockDim.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int idx = threadIdx.x;
    int region = indices[gidx];

    extern __shared__ real h[];
    real *dm11 = &h[size + 1];
    real *dm12i = &h[2 * size + 1];
    real *dm12r = &h[3 * size + 1];
    real *dm22 = &h[4 * size + 1];
    real *e = &h[5 * size + 1];

    h[idx] = gh[gidx];
    if (idx == blockDim.x - 1) {
	h[idx + 1] = gh[gidx + 1];
    }
    dm11[idx] = d[gsize * 0 + gidx];
    dm12i[idx] = d[gsize * 1 + gidx];
    dm12r[idx] = d[gsize * 2 + gidx];
    dm22[idx] = d[gsize * 3 + gidx];
    e[idx] = ge[gidx];

    __syncthreads();

    real dm11_e = dm11[idx];
    real dm12i_e = dm12i[idx];
    real dm12r_e = dm12r[idx];
    real dm22_e = dm22[idx];
    real e_e = e[idx];

    /* execute prediction - correction steps */
    for (int pc_step = 0; pc_step < 4; pc_step++) {

	real rho11 = 0.5 * (dm11[idx] + dm11_e);
	real rho12i = 0.5 * (dm12i[idx] + dm12i_e);
	real rho12r = 0.5 * (dm12r[idx] + dm12r_e);
	real rho22 = 0.5 * (dm22[idx] + dm22_e);
	real OmRabi = gsc[region].d12 * 0.5 * (e[idx] + e_e);

	/* dm11 */
	dm11_e = dm11[idx] + gsc[region].d_t *
	    (-2.0 * OmRabi * rho12i - gsc[region].tau1 * rho11);

	/* imag dm12 */
	dm12i_e = dm12i[idx] + gsc[region].d_t *
	    (- gsc[region].w12 * rho12r + OmRabi * (rho11 - rho22)
	     - gsc[region].gamma12 * rho12i);

	/* real dm12 */
	dm12r_e = dm12r[idx] + gsc[region].d_t *
	    (gsc[region].w12 * rho12i - gsc[region].gamma12 * rho12r);

	/* dm22 */
	dm22_e = dm22[idx] + gsc[region].d_t *
	    (2.0 * OmRabi * rho12i + gsc[region].tau1 * rho11);


	real j = 0; /*TODO revise e *//* ge[gidx] * gsc[region].sigma;*/

	real p_t = gsc[region].M_CP * gsc[region].d12 *
	    (gsc[region].w12 * rho12i - gsc[region].gamma12 * rho12r);

	e_e = e[idx] + gsc[region].M_CE *
	    (-j - p_t + (h[idx + 1] - h[idx])/gsc[region].d_x);

	if (gidx == 0) {
	    /* soft source must be re-implemented, if required */
	    //ge[gidx] += src; /* soft source */
	    e_e = src; /* hard source */
	}
    }

    ge[gidx] = e_e;
    d[gsize * 0 + gidx] = dm11_e;
    d[gsize * 1 + gidx] = dm12i_e;
    d[gsize * 2 + gidx] = dm12r_e;
    d[gsize * 3 + gidx] = dm22_e;
}

/* host members */
SolverCUDA2lvl::SolverCUDA2lvl(const Device& device,
			       const Scenario& scenario) :
    ISolver(device, scenario), comp_maxwell(0), copy(0)
{
    /* total device length */
    Quantity length = device.XDim();

    /* minimum relative permittivity */
    Quantity minRelPermittivity = device.MinRelPermittivity();

    /* TODO: sanity check scenario? */
    if (m_scenario.NumGridPoints % 32 != 0) {
	throw std::invalid_argument("Number of grid points must be multiple"
				    " of 32");
    }

    /* determine grid point and time step size */
    real C = 0.5; /* courant number */
    real velocity = sqrt(MU0() * EPS0() * minRelPermittivity());
    m_scenario.GridPointSize = length()/(m_scenario.NumGridPoints - 1);
    real timestep  = C * m_scenario.GridPointSize * velocity;
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

	sc[i].d_x = m_scenario.GridPointSize;
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
    chk_err(cudaMemcpyToSymbol(gsc, &sc, MaxRegions *
			       sizeof(struct sim_constants)));

    /* initialize streams */
    chk_err(cudaStreamCreate(&comp_maxwell));
    chk_err(cudaStreamCreate(&copy));

    /* allocate space */
    chk_err(cudaMalloc(&m_e, sizeof(real) * m_scenario.NumGridPoints));
    chk_err(cudaMalloc(&m_h, sizeof(real) * (m_scenario.NumGridPoints + 1)));
    chk_err(cudaMalloc(&m_d, sizeof(real) * m_scenario.NumGridPoints * 4));
    chk_err(cudaMalloc(&m_indices, sizeof(unsigned int) *
		       m_scenario.NumGridPoints));

    /* initialize memory */
    unsigned int threads = 128;
    unsigned int blocks = m_scenario.NumGridPoints/threads;

	    init_memory<<<blocks, threads>>>(m_d, m_e, m_h, m_indices);

    /* set up results transfer data structures */
    BOOST_FOREACH(Record rec, m_scenario.Records) {
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

	/* create copy list entry */
	CopyListEntry *entry;
	if (rec.Type == EField) {
	    entry = new CLEField(m_e, res, col_ct, position_idx,
				 interval);
	    m_copyListRed.push_back(entry);
	} else if (rec.Type == HField) {
	    /* TODO: numGridPoints + 1 */
	    entry = new CLEField(m_h, res, col_ct, position_idx,
				 interval);
	    m_copyListBlack.push_back(entry);
	} else if (rec.Type == Density) {
	    if ((rec.I - 1 < 2) && (rec.J - 1 < 2)) {
		if (rec.I == rec.J) {
		    /* main diagonal entry */
		    /*
		    entry = new CLEDensity(m_dm, rec.I - 1, rec.J - 1, res,
					   col_ct, position_idx,
					   interval);*/
		    //m_copyListBlack.push_back(entry);
		    if (rec.I == 2) {
			position_idx += m_scenario.NumGridPoints * 3;
		    }
		    entry = new CLEField(m_d, res, col_ct, position_idx,
					 interval);
		    m_copyListRed.push_back(entry);

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
	} else {
	    // throw exc
	}
    }

    /* sync */
    cudaDeviceSynchronize();
}

SolverCUDA2lvl::~SolverCUDA2lvl()
{
    /* delete copy lists */
    BOOST_FOREACH(CopyListEntry *entry, m_copyListRed) {
	delete entry;
    }
    BOOST_FOREACH(CopyListEntry *entry, m_copyListBlack) {
	delete entry;
    }

    /* free CUDA memory */
    cudaFree(m_h);
    cudaFree(m_e);
    cudaFree(m_d);
    cudaFree(m_indices);

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
SolverCUDA2lvl::getName() const
{
    return factory.getName(); //std::string("CUDA two-level solver");
}

void
SolverCUDA2lvl::run() const
{
    unsigned int threads = 128;
    unsigned int blocks = m_scenario.NumGridPoints/threads;
    /* TODO handle roundoff errors in thread/block partition */

    dim3 block_maxwell(blocks);
    //dim3 block_density(blocks, 2, 2);

    /* main loop */
    for (unsigned int i = 0; i < m_scenario.NumTimeSteps; i++) {

	/* gather h field in copy stream */
	BOOST_FOREACH(CopyListEntry *entry, m_copyListBlack) {
            if (entry->record(i)) {
		/*
		chk_err(cudaMemcpyAsync(entry->getDst(i), entry->getSrc(),
					entry->getSize(),
					cudaMemcpyDeviceToHost, copy));
		*/
	    }
	}

	/* calculate source value */
	/* TODO */
	real f_0 = 2e14;
	real t = i * m_scenario.TimeStepSize;
	real T_p = 20/f_0;
	real gamma = 2 * t/T_p - 1;
	real E_0 = 4.2186e9;
	real src = E_0 * 1/std::cosh(10 * gamma) * sin(2 * M_PI * f_0 * t);
	src = src/2; // pi pulse
	//src = src*2; // 4*pi pulse

	/*
	makestep_e_dm<<<block_maxwell, threads,
	    (6 * threads + 1) * sizeof(real),
	    comp_maxwell>>>(m_d, m_h, m_e, src);
	*/
	makestep_e_dm<<<block_maxwell, threads,
	    (6 * threads + 1) * sizeof(real)>>>(m_d, m_h, m_e, src, m_indices);


	/* sync */
	/*
	chk_err(cudaStreamSynchronize(copy));
	chk_err(cudaStreamSynchronize(comp_maxwell));
	*/

	/* gather e field and dm entries in copy stream */
	BOOST_FOREACH(CopyListEntry *entry, m_copyListRed) {
	    if (entry->record(i)) {
		/*chk_err(cudaMemcpyAsync(entry->getDst(i), entry->getSrc(),
					entry->getSize(),
					cudaMemcpyDeviceToHost, copy));*/
		chk_err(cudaMemcpy(entry->getDst(i), entry->getSrc(),
				   entry->getSize(), cudaMemcpyDeviceToHost));
	    }
	}

	/* makestep_h in maxwell stream */
	/*	makestep_h<<<block_maxwell, threads, (threads + 1) * sizeof(real),
	    comp_maxwell>>>(m_e, m_h);
	*/
	makestep_h<<<block_maxwell, threads, (threads + 1) * sizeof(real)>>>
	    (m_e, m_h, m_indices);

	/* sync */
	/*
	chk_err(cudaStreamSynchronize(copy));
	chk_err(cudaStreamSynchronize(comp_maxwell));*/
    }

    /* sync */
    cudaDeviceSynchronize();
}

}
