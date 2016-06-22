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

void SolverCUDA::setup()
{
}

void SolverCUDA::cleanup()
{
}

void SolverCUDA::run()
{
    makestep_black<<<12, 1>>>();
}

}
