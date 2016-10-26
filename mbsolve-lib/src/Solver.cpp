#include <Solver.hpp>

#include <boost/foreach.hpp>

namespace mbsolve {

std::map<std::string, ISolverFactory *>
Solver::m_factories;

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

    /* allocate space for results */
    unsigned int N_t = m_solver->getScenario().NumTimeSteps;
    /* m_solver->getScenario == scenario? */
    unsigned int N_x = m_solver->getScenario().NumGridPoints;

/*    std::vector<Record>::iterator rec_it;
    for (rec_it = scenario.Records.begin(); rec_it != scenario.Records.end();
	 rec_it++) {

	/* calculate size */

	/* allocate */

	/* link Record <-> Result (implicetely by index?) */

	/* push_back */

//}

}

Solver::~Solver()
{
    /* clean up results */
    BOOST_FOREACH(Result *result, m_results) {
	delete result;
    }
    /* clean up solver */
    delete m_solver;
}

std::string
Solver::getName() const
{
    return m_solver->getName();
}

void
Solver::run()
{
    m_solver->run(m_results);
}

const std::vector<Result *>&
Solver::getResults() const
{
    return m_results;
}

}
