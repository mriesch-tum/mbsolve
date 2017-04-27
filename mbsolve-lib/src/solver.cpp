/*
 * mbsolve: Framework for solving the Maxwell-Bloch/-Lioville equations
 *
 * Copyright (c) 2016, Computational Photonics Group, Technical University of
 * Munich.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#include <solver.hpp>

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
