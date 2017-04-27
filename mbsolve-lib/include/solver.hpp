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

#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <device.hpp>
#include <scenario.hpp>
#include <types.hpp>

namespace mbsolve {

class ISolver
{
protected:
    Scenario m_scenario;
    Device m_device;

    std::vector<Result *> m_results;
public:
    ISolver(const Device& device, const Scenario& scenario) :
	m_device(device), m_scenario(scenario) { }

    virtual ~ISolver();

    const Scenario& getScenario() const { return m_scenario; }

    const Device& getDevice() const { return m_device; }

    virtual std::string getName() const = 0;

    virtual void run() const = 0;

    const std::vector<Result *>& getResults() const
    {
	return m_results;
    }
};

class ISolverFactory;

class Solver
{
private:
    static std::map<std::string, ISolverFactory *> m_factories;
    ISolver *m_solver;

public:
    Solver(const std::string& name, const Device& device,
	   const Scenario& scenario);

    ~Solver();

    std::string getName() const;

    const Scenario& getScenario() const { return m_solver->getScenario(); }

    const Device& getDevice() const { return m_solver->getDevice(); }

    void run() const;

    const std::vector<Result *>& getResults() const;

    static void registerFactory(const std::string& name,
				ISolverFactory *factory);
};

class ISolverFactory
{
public:
    ISolverFactory() { }
    virtual ~ISolverFactory() { }
    virtual ISolver* createInstance(const Device& device,
				    const Scenario& scenario) const = 0;
};

template<typename T>
class SolverFactory : ISolverFactory
{
private:
    std::string m_name;
public:
    explicit SolverFactory(const std::string& name) : m_name(name) {
        Solver::registerFactory(name, this);
    }

    ISolver* createInstance(const Device& device,
			    const Scenario& scenario) const {
	return new T(device, scenario);
    }

    const std::string& getName() const { return m_name; }
};

}

#endif
