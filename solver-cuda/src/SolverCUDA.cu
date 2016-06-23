#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <SolverCUDA.hpp>

namespace mbsolve {

__global__ void makestep_black()
{

}

SolverCUDA::SolverCUDA() : Solver("CUDA Solver")
{
}

SolverCUDA::~SolverCUDA()
{
}

void SolverCUDA::do_setup(const Device& device, const Scenario& scenario)
{
}

void SolverCUDA::do_cleanup()
{
}

void SolverCUDA::do_run()
{
    makestep_black<<<12, 1>>>();
}

}
