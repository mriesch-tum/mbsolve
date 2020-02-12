/*
 * mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
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

#include <mbsolve/lib/solver.hpp>

namespace mbsolve {

std::map<std::string, solver::bootstrap_t> solver::m_bootstraps;

std::shared_ptr<solver>
solver::create_instance(
    const std::string& name,
    std::shared_ptr<const device> dev,
    std::shared_ptr<scenario> scen)
{
    auto it = m_bootstraps.find(name);
    if (it == m_bootstraps.end()) {
        throw std::invalid_argument("Unknown solver " + name);
    }
    return it->second(dev, scen);
}

void
solver::register_bootstrap(const std::string& name, bootstrap_t b)
{
    if (m_bootstraps[name]) {
        throw std::invalid_argument("Solver already registered.");
    }
    m_bootstraps[name] = b;
}

solver::solver(
    const std::string& name,
    std::shared_ptr<const device> dev,
    std::shared_ptr<scenario> scen)
  : m_name(name), m_device(dev), m_scenario(scen)
{}

solver::~solver() {}

const std::string&
solver::get_name() const
{
    return m_name;
}

const scenario&
solver::get_scenario() const
{
    return *m_scenario;
}

const device&
solver::get_device() const
{
    return *m_device;
}

const std::vector<std::shared_ptr<result> >&
solver::get_results() const
{
    return m_results;
}

std::vector<std::string>
solver::get_avail_solvers()
{
    std::vector<std::string> solvers;

    for (const auto& s : m_bootstraps) {
        solvers.push_back(s.first);
    }

    return solvers;
}
}
