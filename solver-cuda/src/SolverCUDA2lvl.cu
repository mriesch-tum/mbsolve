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
__device__ __inline__ unsigned int get_region(unsigned int idx)
{
    for (unsigned int i = 0; i < MaxRegions; i++) {
	if (idx < gsc[i].idx_end) {
	    return i;
	}
    }
    return 0;
}

/* TODO: initialize may be reused by other CUDA solvers, make general */
/* TODO: region-wise? initialization */
__global__ void init_memory(CUDADensityMatrixData dm, real *e, real *h)
{
    unsigned int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int max = blockDim.x * gridDim.x - 1;
    int region = get_region(gidx);

    /* TODO: alternative initializations */
    curandState_t rand_state;

    /* initialize random number generator */
    curand_init(clock64(), gidx, 0, &rand_state);

    if (gidx == max - 1) {
	h[max] = 0.0;
    }
    h[gidx] = 0.0;
    e[gidx] = 0.0;
    //   e[gidx] = curand_uniform(&rand_state) * 1e-15;

    real trace = 0.0;

    for (unsigned int row = 0; row < dm.getNumLevels(); row++) {
	for (unsigned int col = 0; col < dm.getNumLevels(); col++) {
	    if (row == col) {
		/*
		real entry = curand_uniform(&rand_state);
		dm.oldDM(row, col)[gidx] = entry;
		trace += entry;
		*/
	    } else {
		dm.oldDM(row, col)[gidx] = 0.0;
	    }
	    for (int i = 0; i < dm.getNumMultistep(); i++) {
		dm.rhs(row, col, i)[gidx] = 0.0;
	    }
	}
    }

    dm.oldDM(0, 0)[gidx] = gsc[region].dm11_init;
    dm.oldDM(1, 1)[gidx] = gsc[region].dm22_init;

    /*
    for (unsigned int i = 0; i < dm.getNumLevels(); i++) {
	dm.oldDM(i, i)[gidx] /= trace;
    }
    */
}

__global__ void makestep_h(const real *ge, real *gh)
{
    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);

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

__global__ void makestep_e(CUDADensityMatrixData dm, const real *gh,
			   real *ge, real src)
{
    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);

    extern __shared__ real h[];

    h[idx] = gh[gidx];
    if (idx == blockDim.x - 1) {
	h[idx + 1] = gh[gidx + 1];
    }

    __syncthreads();

    real j = ge[gidx] * gsc[region].sigma;
    /*
      real p_t = gsc[region].M_CP * gsc[region].d12 * dm.rhs(0, 1, 1)[gidx];


    */
    // real p_t = gsc[region].M_CP * gsc[region].d12
    real p_t = gsc[region].M_CP * gsc[region].d12 * dm.rhs(1, 0, 1)[gidx];

    ge[gidx] += gsc[region].M_CE *
	(-j - p_t + (h[idx + 1] - h[idx])/gsc[region].d_x);

    if (gidx == 0) {
	ge[gidx] += src;
    }
}

__global__ void makestep_dm(CUDADensityMatrixData dm, const real *ge)
{
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);
    int row = blockIdx.y;
    int col = blockIdx.z;

    real rhs = 0.0;

    if ((row == 0) && (col == 0)) {
	/* dm11 */
	rhs = - dm.oldDM(0, 1)[gidx] * 2.0 * gsc[region].d12 * ge[gidx]
	    - dm.oldDM(0, 0)[gidx] * gsc[region].tau1;
    } else if ((row == 0) && (col == 1)) {
	/* imag dm12 */
	rhs = - dm.oldDM(1, 0)[gidx] * gsc[region].w12
	    + (dm.oldDM(0, 0)[gidx] - dm.oldDM(1, 1)[gidx]) * gsc[region].d12 * ge[gidx]
	    - dm.oldDM(0, 1)[gidx] * gsc[region].gamma12;
    } else if ((row == 1) && (col == 0)) {
	/* real dm12 */
	rhs = + dm.oldDM(0, 1)[gidx] * gsc[region].w12
	    - dm.oldDM(1, 0)[gidx] * gsc[region].gamma12;
    } else if ((row == 1) && (col == 1)) {
	/* dm22 */
	rhs = dm.oldDM(0, 1)[gidx] * 2.0 * gsc[region].d12 * ge[gidx]
	    + dm.oldDM(0, 0)[gidx] * gsc[region].tau1;
    } else {
	/* do nothing */
    }

    dm.rhs(row, col, 0)[gidx] = rhs;
    dm.newDM(row, col)[gidx] = dm.oldDM(row, col)[gidx] + gsc[region].d_t *
	(+ dm.rhs(row, col, 0)[gidx] * 1901.0/720.0
	 - dm.rhs(row, col, 1)[gidx] * 1387.0/360.0
	 + dm.rhs(row, col, 2)[gidx] * 109.0/30.0
	 - dm.rhs(row, col, 3)[gidx] * 637.0/360.0
	 + dm.rhs(row, col, 4)[gidx] * 251.0/720.0);
}

/* host members */
SolverCUDA2lvl::SolverCUDA2lvl(const Device& device,
			       const Scenario& scenario) :
    ISolver(device, scenario), comp_maxwell(0), comp_density(0), copy(0)
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
    real C = 0.9; /* courant number */
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
    chk_err(cudaStreamCreate(&comp_density));
    chk_err(cudaStreamCreate(&copy));

    /* allocate space */
    chk_err(cudaMalloc(&m_e, sizeof(real) * m_scenario.NumGridPoints));
    chk_err(cudaMalloc(&m_h, sizeof(real) * (m_scenario.NumGridPoints + 1)));
    m_dm = new CUDADensityMatrix(m_scenario.NumGridPoints, 2, 5);

    /* initialize memory */
    unsigned int threads = 128;
    unsigned int blocks = m_scenario.NumGridPoints/threads;

    init_memory<<<blocks, threads>>>(m_dm->getData(), m_e, m_h);

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
		    entry = new CLEDensity(m_dm, rec.I - 1, rec.J - 1, res,
					   col_ct, position_idx,
					   interval);
		    m_copyListBlack.push_back(entry);
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

    /* delete density matrix */
    delete(m_dm);

    /* free CUDA memory */
    cudaFree(m_h);
    cudaFree(m_e);

    /* clean up streams */
    if (comp_maxwell) {
	cudaStreamDestroy(comp_maxwell);
    }
    if (comp_density) {
	cudaStreamDestroy(comp_density);
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
    dim3 block_density(blocks, 2, 2);

    /* main loop */
    for (unsigned int i = 0; i < m_scenario.NumTimeSteps; i++) {
	/* makestep_dm in density stream */
	makestep_dm<<<block_density, threads, 0, comp_density>>>
	    (m_dm->getData(), m_e);

	/* makestep_h in maxwell stream */
	makestep_h<<<block_maxwell, threads + 1, (threads + 1) * sizeof(real),
	    comp_maxwell>>>(m_e, m_h);

	/* gather e field in copy stream */
	BOOST_FOREACH(CopyListEntry *entry, m_copyListRed) {
	    if (entry->record(i)) {
		chk_err(cudaMemcpyAsync(entry->getDst(i), entry->getSrc(),
					entry->getSize(),
					cudaMemcpyDeviceToHost, copy));
	    }
	}

	/* sync */
	chk_err(cudaStreamSynchronize(copy));
	chk_err(cudaStreamSynchronize(comp_maxwell));
	chk_err(cudaStreamSynchronize(comp_density));

	/* switch density matrix double buffer */
	m_dm->getData().next();

	/* calculate source value */
	/* TODO */
	real f_0 = 2e14;
	real t = i * m_scenario.TimeStepSize;
	real T_p = 20/f_0;
	real gamma = 2 * t/T_p - 1;
	real E_0 = 4.2186e9;
	real src = E_0 * 1/std::cosh(10 * gamma) * sin(2 * M_PI * f_0 * t);

	/* makestep_e kernel */
	makestep_e<<<block_maxwell, threads + 1, (threads + 1) * sizeof(real),
	    comp_maxwell>>>(m_dm->getData(), m_h, m_e, src);

	/* gather h field and dm entries in copy stream */
	BOOST_FOREACH(CopyListEntry *entry, m_copyListBlack) {
            if (entry->record(i)) {
		chk_err(cudaMemcpyAsync(entry->getDst(i), entry->getSrc(),
					entry->getSize(),
					cudaMemcpyDeviceToHost, copy));
	    }
	}

	/* sync */
	chk_err(cudaStreamSynchronize(copy));
	chk_err(cudaStreamSynchronize(comp_maxwell));
	chk_err(cudaStreamSynchronize(comp_density));
    }

    /* sync */
    cudaDeviceSynchronize();
}

}
