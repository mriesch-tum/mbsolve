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
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    int type = blockIdx.y;
    int max = blockDim.x * gridDim.x - 1;

    /* TODO: alternative initializations */

    if (type < NumEntries) {
	dm.OldDM(type)[idx] = 0.0;
	for (int i = 0; i < NumMultistep; i++) {
	    dm.RHS(type, i)[idx] = 0.0;
	}
    } else if (type == NumEntries) {
	e[idx] = 0.0;
    } else if (type == NumEntries + 1) {
	if (idx == max - 1) {
	    h[max] = 0.0;
	}
	h[idx] = 0.0;
    } else {
	// handle error
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
    real p_t = 0.0; /* TODO: ? */

    real p_t = gsc[region].C_P * gsc[region].d12 * dm.OldDM(2)[gidx];

    ge[gidx] += gsc[region].M_CE *
	(-j - p_t + (h[idx + 1] - h[idx])/gsc[region].d_x);
}

__global__ void makestep_dm(const DensityMatrix& dm, const real *ge)
{
    //    int idx = threadIdx.x;
    int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    int region = get_region(gidx);
    int type = blockIdx.y;

    real rhs = 0.0;

    /* if blah
       depending on dm entry
     */

    dm.RHS(type, 0)[gidx] = rhs;
    dm.NewDM(type)[gidx] = dm.OldDM(type)[gidx] + gsc[region].d_t *
	(+ dm.RHS(type, 0)[gidx] * 1901.0/720.0
	 - dm.RHS(type, 1)[gidx] * 1387.0/360.0
	 + dm.RHS(type, 2)[gidx] * 109.0/30.0
	 - dm.RHS(type, 3)[gidx] * 637.0/360.0
	 + dm.RHS(type, 4)[gidx] * 251.0/720.0);
}


DensityMatrix::DensityMatrix() : a_is_old(true), head(0)
{
}

DensityMatrix::~DensityMatrix()
{
    for (unsigned int i = 0; i < NumEntries; i++) {
	cudaFree(dm_a[i]);
	cudaFree(dm_b[i]);
	for (unsigned int j = 0; j < NumMultistep; j++) {
	    cudaFree(rhs[i][j]);
	}
    }
}

__device__ __inline__ real *
DensityMatrix::OldDM(unsigned int entry) const
{
    return a_is_old ? dm_a[entry] : dm_b[entry];
}

__device__ __inline__ real *
DensityMatrix::NewDM(unsigned int entry) const
{
    return a_is_old ? dm_b[entry] : dm_a[entry];
}

__device__ __inline__ real *
DensityMatrix::RHS(unsigned int entry, unsigned int row) const
{
    return rhs[entry][(row + head) % NumMultistep];
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
    for (unsigned int i = 0; i < NumEntries; i++) {
	chk_err(cudaMalloc(&dm_a[i], sizeof(real) * numGridPoints));
	chk_err(cudaMalloc(&dm_b[i], sizeof(real) * numGridPoints));
	for (unsigned int j = 0; j < NumMultistep; j++) {
	    chk_err(cudaMalloc(&rhs[i][j], sizeof(real) * numGridPoints));
	}
    }
}

/* host members */
SolverCUDA2lvl::SolverCUDA2lvl() : Solver("CUDA 2-Level Solver")
{
}

SolverCUDA2lvl::~SolverCUDA2lvl()
{
}

/* TODO: RAII with CUDA pointers/memory? */
void SolverCUDA2lvl::do_setup(const Device& device, Scenario& scenario)
{
    /* total device length */
    Quantity length = device.XDim();

    /* minimum relative permittivity */
    Quantity minRelPermittivity = device.MinRelPermittivity();

    /* determine grid point and time step size */
    real C = 0.9; /* courant number */
    real velocity = sqrt(MU0() * EPS0() * minRelPermittivity());
    scenario.GridPointSize = length()/(scenario.NumGridPoints - 1);
    real timestep  = C * scenario.GridPointSize * velocity;
    scenario.NumTimeSteps = ceil(scenario.SimEndTime/timestep) + 1;
    scenario.TimeStepSize = scenario.SimEndTime/(scenario.NumTimeSteps - 1);

    /* determine border indices and initialize region settings */
    if (device.Regions.size() > MaxRegions) {
	throw std::invalid_argument("Too many regions requested");
    }
    struct sim_constants sc[MaxRegions];

    std::vector<Region>::const_iterator it;
    unsigned int i;
    for (it = device.Regions.begin(), i = 0; it != device.Regions.end();
	 it++, i++) {
	if (i > 0) {
	    sc[i - 1].idx_end = round(it->X0()/scenario.GridPointSize) - 1;
	}
	sc[i].idx_start = round(it->X0()/scenario.GridPointSize);
	sc[i].M_CE = scenario.TimeStepSize/(EPS0() * it->RelPermittivity());
	sc[i].M_CH = scenario.TimeStepSize/(MU0() * scenario.GridPointSize);
	sc[i].M_CP = -2.0 * it->DopingDensity * E0;
	sc[i].sigma = 2.0 * sqrt(EPS0 * it->RelPermittivity/MU0) * it->Losses;

	sc[i].w12 = (it->TransitionFrequencies.size() < 1) ? 0.0 :
	    it->TransitionFrequencies[i]();
	sc[i].d12 = (it->DipoleMoments.size() < 1) ? 0.0 :
	    it->DipoleMoments[i]();
	sc[i].gamma1 = (it->ScatteringRates.size() < 1) ? 0.0 :
	    it->ScatteringRates[i]();
	sc[i].gamma2 = (it->DephasingRates.size() < 1) ? 0.0 :
	    it->DephasingRates[i]();

	sc[i].d_x = scenario.GridPointSize;
	sc[i].d_t = scenario.TimeStepSize;
    }
    if (i > 0) {
	sc[i - 1].idx_end = scenario.NumGridPoints - 1;
    }

    /* initialize streams */
    chk_err(cudaStreamCreate(&comp_maxwell));
    chk_err(cudaStreamCreate(&comp_density));
    chk_err(cudaStreamCreate(&copy));

    /* allocate space */
    chk_err(cudaMalloc(&e, sizeof(real) * scenario.NumGridPoints));
    chk_err(cudaMalloc(&h, sizeof(real) * (scenario.NumGridPoints + 1)));
    dm.initialize(scenario.NumGridPoints);

    /* initalize memory */
    /* TODO: kernel call */

    /* copy settings to CUDA constant memory */
    chk_err(cudaMemcpyToSymbol(gsc, &sc, MaxRegions *
			       sizeof(struct sim_constants)));
}

void SolverCUDA2lvl::do_cleanup()
{
    /* free CUDA memory */
    cudaFree(h);
    cudaFree(e);

    /* clean up streams */
    cudaStreamDestroy(comp_maxwell);
    cudaStreamDestroy(comp_density);
    cudaStreamDestroy(copy);

    cudaDeviceReset();
}

void SolverCUDA2lvl::do_run(std::vector<Result *>& results)
{


    int threads = 1024;
    int blocks = 10; // NumGridPoint / threads
    /* TODO handle roundoff errors */

    dim3 block(blocks);
    dim3 thread(threads);

    /* for i = 2:N_t */
    /* makestep_h in maxwell stream */
    /* makestep_dm in density stream */
    /* gather e field in copy stream */
    /* call toggle */
    /* sync */
    /* calculate source value -> makestep_e kernel */
    /* gather h field and dm entries in copy stream */
    /* makestep_e */
    /* sync */

    /* end for */
    makestep_h<<<block, thread, sizeof(real) * (threads + 1)>>>(e, h);

}

}
