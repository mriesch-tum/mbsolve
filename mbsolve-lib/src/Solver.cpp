#include <Solver.hpp>

namespace mbsolve {

std::map<std::string, ISolverFactory *> Solver::m_solvers;

void
Solver::register_new(const std::string& name, ISolverFactory *factory)
{
    // TODO: check doubles
    m_solvers[name] = factory;
}

Solver *
Solver::create(const std::string& name)
{
    std::map<std::string, ISolverFactory *>::iterator it;
    it = m_solvers.find(name);
    if (it == m_solvers.end()) {
	throw std::invalid_argument("Unknown solver " + name);
    }
    return it->second->create();
}

Solver::Solver(std::string name) : m_initialized(false), m_name(name)
{
}

Solver::~Solver()
{
    if (m_initialized) {
	do_cleanup();
    }
}

const std::string&
Solver::name()
{
    return m_name;
}

void
Solver::setup(const Device& device, Scenario& scenario)
{
    do_setup(device, scenario);
    m_initialized = true;
}

void
Solver::run(std::vector<Result *>& results)
{
    return do_run(results);
}

}
