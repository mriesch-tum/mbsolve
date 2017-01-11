#include <solver_2lvl_pc.hpp>

namespace mbsolve{

static SolverFactory<SolverOMP_2lvl_pc> factory("openmp-2lvl-pc");

SolverOMP_2lvl_pc::SolverOMP_2lvl_pc(const Device& device,
				     const Scenario& scenario) :
    ISolver(device, scenario)
{

}

SolverOMP_2lvl_pc::~SolverOMP_2lvl_pc()
{

}

std::string
SolverOMP_2lvl_pc::getName() const
{
    return factory.getName();
}

void
SolverOMP_2lvl_pc::run() const
{


}

}
