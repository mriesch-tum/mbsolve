#include <SolverGeneric.hpp>

namespace mbsolve {

static SolverFactory<SolverGeneric> factory("Generic");

SolverGeneric::SolverGeneric() : Solver("Generic")
{
}

SolverGeneric::~SolverGeneric()
{
}

void SolverGeneric::do_setup(const Device& device, const Scenario& scenario)
{
}

void SolverGeneric::do_cleanup()
{
}

void SolverGeneric::do_run(std::vector<Result *>& results)
{

}

}
