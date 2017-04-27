#include <solver_generic.hpp>

namespace mbsolve {

static SolverFactory<solver_generic> factory("generic");

solver_generic::solver_generic() : solver("generic")
{
}

solver_generic::~solver_generic()
{
}

void
solver_generic::do_setup(const device& dev, const scenario& scen)
{
}

void solver_generic::do_cleanup()
{
}

void solver_generic::do_run(std::vector<Result *>& results)
{

}

}
