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

namespace mbsolve {

std::map<std::string, solver_factory_int *>
solver_int::m_factories;

solver_int::solver_int(std::shared_ptr<const device> dev,
                       std::shared_ptr<scenario> scen) :
    m_device(dev), m_scenario(scen)
{
}

solver_int::~solver_int()
{
}

void
solver_int::register_factory(const std::string& name,
                             solver_factory_int *factory)
{
    if (m_factories[name]) {
	throw std::invalid_argument("Solver already registered.");
    }
    m_factories[name] = factory;
}

solver_factory_int *
solver_int::find_factory(const std::string& name)
{
    auto it = m_factories.find(name);
    if (it == m_factories.end()) {
        throw std::invalid_argument("Unknown solver " + name);
    }
    return it->second;
}

solver::solver(const std::string& name, std::shared_ptr<const device> dev,
	       std::shared_ptr<scenario> scen)
{
    /* create solver */
    solver_factory_int *factory = solver_int::find_factory(name);
    m_solver = factory->create_instance(dev, scen);
}

solver::~solver()
{
}

const std::string&
solver::get_name() const
{
    return m_solver->get_name();
}

void
solver::run() const
{
    m_solver->run();
}

const std::vector<std::shared_ptr<result> >&
solver::get_results() const
{
    return m_solver->get_results();
}

}
