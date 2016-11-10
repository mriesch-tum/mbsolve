#include <boost/foreach.hpp>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <SolverCUDA2lvl.hpp>

namespace mbsolve {

static SolverFactory<SolverCUDA2lvl> factory("cuda-2lvl");

static inline void chk_err(cudaError_t code)
{
    if (code != cudaSuccess) {
	throw std::runtime_error(std::string("CUDA: ") +
				 cudaGetErrorString(code));
    }
}

/* CUDA memory and kernels */
__device__ __constant__ struct sim_constants gsc[MaxRegions];

__device__ __inline__ unsigned int get_region(unsigned int idx)
{
    for (unsigned int i = 0; i < MaxRegions; i++) {
	if (idx < gsc[i].idx_end) {
	    return i;
	}
    }
    return 0;
}

__global__ void init_memory(const DensityMatrix& dm, real *e, real *h)
{
    unsigned int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int max = blockDim.x * gridDim.x - 1;

    /* TODO: alternative initializations */
    curandState_t rand_state;

    /* initialize random number generator */
    curand_init(clock64(), gidx, 0, &rand_state);

    if (gidx == max - 1) {
	h[max] = 0.0;
    }
    h[gidx] = 0.0;
    e[gidx] = curand_uniform(&rand_state) * 1e-15;

    real populations[NumLevels];
    real trace = 0.0;

    for (unsigned int i = 0; i < NumLevels; i++) {
	populations[i] = curand_uniform(&rand_state);
	trace += populations[i];
    }

    for (unsigned int row = 0; row < NumLevels; row++) {
	for (unsigned int col = 0; col < NumLevels; col++) {
	    if (row == col) {
		if (row == NumLevels - 1) {
		    dm.OldDM(row, col)[gidx] = populations[row]/trace;
		} else {
		    dm.OldDM(row, col)[gidx] = 0.0;
		}
		for (int i = 0; i < NumMultistep; i++) {
		    dm.RHS(row, col, i)[gidx] = 0.0;
		}
	    }
	}
    }
}

__global__ void makestep_h(const real *ge, real *gh)
{
    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);

    extern __shared__ real e[];

    if ((idx == 0) && (gidx != 0)) {
	e[0] = ge[gidx - 1];
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

__global__ void makestep_e(const DensityMatrix& dm, const real *gh, real *ge)
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
    real p_t = gsc[region].M_CP * gsc[region].d12 * dm.OldDM(0, 1)[gidx];

    ge[gidx] += gsc[region].M_CE *
	(-j - p_t + (h[idx + 1] - h[idx])/gsc[region].d_x);
}

__global__ void makestep_dm(const DensityMatrix& dm, const real *ge)
{
    //    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);
    int row = blockIdx.y;
    int col = blockIdx.z;

    real rhs = 0.0;

    if ((row == 0) && (col == 0)) {
	/* dm11 */
	rhs = - dm.OldDM(0, 1)[gidx] * 2.0 * gsc[region].d12 * ge[gidx]
	    - dm.OldDM(0, 0)[gidx] * gsc[region].tau1;
    } else if ((row == 0) && (col == 1)) {
	/* imag dm12 */
	rhs = dm.OldDM(1, 0)[gidx] * gsc[region].w12
	    - dm.OldDM(0, 1)[gidx] * gsc[region].gamma12;
    } else if ((row == 1) && (col == 0)) {
	/* real dm12 */
	rhs = - dm.OldDM(0, 1)[gidx] * gsc[region].w12
	    - dm.OldDM(1, 0)[gidx] * gsc[region].gamma12;
    } else if ((row == 1) && (col == 1)) {
	/* dm22 */
	rhs = dm.OldDM(0, 1)[gidx] * 2.0 * gsc[region].d12 * ge[gidx]
	    + dm.OldDM(0, 0)[gidx] * gsc[region].tau1;
    } else {
	/* do nothing */
    }

    dm.RHS(row, col, 0)[gidx] = rhs;
    dm.NewDM(row, col)[gidx] = dm.OldDM(row, col)[gidx] + gsc[region].d_t *
	(+ dm.RHS(row, col, 0)[gidx] * 1901.0/720.0
	 - dm.RHS(row, col, 1)[gidx] * 1387.0/360.0
	 + dm.RHS(row, col, 2)[gidx] * 109.0/30.0
	 - dm.RHS(row, col, 3)[gidx] * 637.0/360.0
	 + dm.RHS(row, col, 4)[gidx] * 251.0/720.0);
}


DensityMatrix::DensityMatrix() : a_is_old(true), head(0)
{
}

DensityMatrix::~DensityMatrix()
{
    for (unsigned int i = 0; i < NumLevels; i++) {
	for (unsigned int j = 0; j < NumLevels; j++) {
	    cudaFree(dm_a[i][j]);
	    cudaFree(dm_b[i][j]);
	    for (unsigned int k = 0; k < NumMultistep; k++) {
		cudaFree(rhs[i][j][k]);
	    }
	}
    }
}

__device__ __inline__ real *
DensityMatrix::OldDM(unsigned int row, unsigned int col) const
{
    return a_is_old ? dm_a[row][col] : dm_b[row][col];
}

__device__ __inline__ real *
DensityMatrix::NewDM(unsigned int row, unsigned int col) const
{
    return a_is_old ? dm_b[row][col] : dm_a[row][col];
}

__device__ __inline__ real *
DensityMatrix::RHS(unsigned int row, unsigned int col,
		   unsigned int rhsIdx) const
{
    return rhs[row][col][(rhsIdx + head) % NumMultistep];
}

void
DensityMatrix::next()
{
    a_is_old = !a_is_old;
    head = (head + 1) % NumMultistep;
}

void
DensityMatrix::initialize(unsigned int numGridPoints)
{
    for (unsigned int i = 0; i < NumLevels; i++) {
	for (unsigned int j = 0; j < NumLevels; j++) {
	    chk_err(cudaMalloc(&dm_a[i][j], sizeof(real) * numGridPoints));
	    chk_err(cudaMalloc(&dm_b[i][j], sizeof(real) * numGridPoints));
	    for (unsigned int k = 0; k < NumMultistep; k++) {
		chk_err(cudaMalloc(&rhs[i][j][k],
				   sizeof(real) * numGridPoints));
	    }
	}
    }
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

    unsigned int i;
    BOOST_FOREACH(Region reg, device.Regions) {
	if (i > 0) {
	    sc[i - 1].idx_end = round(reg.X0()/m_scenario.GridPointSize) - 1;
	}
	sc[i].idx_start = round(reg.X0()/m_scenario.GridPointSize);
	sc[i].M_CE = m_scenario.TimeStepSize/(EPS0() * reg.RelPermittivity());
	sc[i].M_CH = m_scenario.TimeStepSize/(MU0() *
					      m_scenario.GridPointSize);
	sc[i].M_CP = -2.0 * reg.DopingDensity * E0;
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
	i++;
    }
    if (i > 0) {
	sc[i - 1].idx_end = m_scenario.NumGridPoints - 1;
    }

    /* initialize streams */
    chk_err(cudaStreamCreate(&comp_maxwell));
    chk_err(cudaStreamCreate(&comp_density));
    chk_err(cudaStreamCreate(&copy));

    /* allocate space */
    chk_err(cudaMalloc(&e, sizeof(real) * m_scenario.NumGridPoints));
    chk_err(cudaMalloc(&h, sizeof(real) * (m_scenario.NumGridPoints + 1)));
    dm.initialize(m_scenario.NumGridPoints);

    /* initalize memory */
    /* TODO: kernel call */

    /* copy settings to CUDA constant memory */
    chk_err(cudaMemcpyToSymbol(gsc, &sc, MaxRegions *
			       sizeof(struct sim_constants)));

    /* set up results transfer data structures */
    BOOST_FOREACH(Record rec, m_scenario.Records) {
	unsigned int interval = ceil(rec.Interval()/m_scenario.TimeStepSize);
	unsigned int row_ct = m_scenario.NumTimeSteps/interval;
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
	    entry = new CLEField(e, res, col_ct * row_ct, position_idx,
				 interval);
	    m_copyListRed.push_back(entry);
	} else if (rec.Type == HField) {
	    entry = new CLEField(h, res, col_ct * row_ct, position_idx,
				 interval);
	    m_copyListBlack.push_back(entry);
	} else if (rec.Type == Density) {
	    if ((rec.I - 1 < 2) && (rec.J - 1 < 2)) {
		if (rec.I == rec.J) {
		    /* main diagonal entry */
		    entry = new CLEDensity(&dm, rec.I, rec.J, res,
					   col_ct * row_ct, position_idx,
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

    /* delete results */
    BOOST_FOREACH(Result *result, m_results) {
	delete result;
    }

    /* free CUDA memory */
    cudaFree(h);
    cudaFree(e);

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
    return std::string("CUDA two-level solver");
}

void
SolverCUDA2lvl::run(const std::vector<Result *>& results) const
{
    unsigned int threads = 128;
    unsigned int blocks = 10; // NumGridPoint / threads
    /* TODO handle roundoff errors in thread/block partition */

    dim3 block(blocks);
    dim3 thread(threads);

    /* main loop */
    for (unsigned int i = 1; i < m_scenario.NumTimeSteps; i++) {
	/* makestep_h in maxwell stream */
	/* makestep_dm in density stream */
	/* gather e field in copy stream */
	BOOST_FOREACH(CopyListEntry *entry, m_copyListRed) {
	    if (entry->record(i)) {
		cudaMemcpyAsync(entry->getDst(i), entry->getSrc(),
				entry->getSize(), cudaMemcpyDeviceToHost,
				copy);
	    }
	}


	/* sync */

	/* call toggle */

	/* calculate source value -> makestep_e kernel */

	/* gather h field and dm entries in copy stream */
	BOOST_FOREACH(CopyListEntry *entry, m_copyListBlack) {
	    if (entry->record(i)) {
		cudaMemcpyAsync(entry->getDst(i), entry->getSrc(),
				entry->getSize(), cudaMemcpyDeviceToHost,
				copy);
	    }
	}


	/* makestep_e */
	/* sync */

    }

    //    makestep_h<<<block, thread, sizeof(real) * (threads + 1)>>>(e, h);

}

}
