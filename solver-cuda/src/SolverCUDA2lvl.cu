#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <SolverCUDA2lvl.hpp>

namespace mbsolve {

static SolverFactory<SolverCUDA2lvl> factory("cuda-2lvl");

/* vector of different settings struct */

__global__ void makestep_black()
{

}

SolverCUDA2lvl::SolverCUDA2lvl() : Solver("CUDA 2-Level Solver")
{
}

SolverCUDA2lvl::~SolverCUDA2lvl()
{
}

void SolverCUDA2lvl::do_setup(const Device& device, const Scenario& scenario)
{

    /* overall length */

    /* determine grid point size */
    /* determine time step size */


    /* determine border indices */

    /* allocate space */

}

void SolverCUDA2lvl::do_cleanup()
{
}

void SolverCUDA2lvl::do_run(std::vector<Result *>& results)
{

    /* invoke kernels with different settings index */

    makestep_black<<<12, 1>>>();


}

}
