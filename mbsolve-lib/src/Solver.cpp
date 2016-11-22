#include <Solver.hpp>

#include <boost/foreach.hpp>

namespace mbsolve {

std::map<std::string, ISolverFactory *>
Solver::m_factories;

ISolver::~ISolver()
{
    /* clean up results */
    BOOST_FOREACH(Result *result, m_results) {
	delete result;
    }
}

void
Solver::registerFactory(const std::string& name, ISolverFactory *factory)
{
    if (m_factories[name]) {
	throw std::invalid_argument("Solver already registered.");
    }
    m_factories[name] = factory;
}

Solver::Solver(const std::string& name, const Device& device,
	       const Scenario& scenario)
{
    /* create solver */
    std::map<std::string, ISolverFactory *>::iterator it;
    it = m_factories.find(name);
    if (it == m_factories.end()) {
	throw std::invalid_argument("Unknown solver " + name);
    }
    m_solver = it->second->createInstance(device, scenario);
}

Solver::~Solver()
{
    /* clean up solver */
    delete m_solver;
}

std::string
Solver::getName() const
{
    return m_solver->getName();
}

void
Solver::run() const
{
    m_solver->run();
}

const std::vector<Result *>&
Solver::getResults() const
{
    return m_solver->getResults();
}

}
