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

#ifndef MBSOLVE_SOLVER_H
#define MBSOLVE_SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <device.hpp>
#include <scenario.hpp>
#include <solver_int.hpp>
#include <types.hpp>

namespace mbsolve {

/**
 * This class provides the interface to create an instance of a solver
 * implementation. Each implementation is a subclass of \ref solver_int and
 * is created internally.
 * \ingroup MBSOLVE_LIB
 */
class solver
{
private:
    solver_int *m_solver;

public:
    solver(const std::string& name, device * const dev,
	   const Scenario& scenario);

    ~solver();

    std::string get_name() const;

    const Scenario& get_scenario() const { return m_solver->get_scenario(); }

    //const Device& getDevice() const { return m_solver->getDevice(); }

    void run() const;

    const std::vector<Result *>& get_results() const;

};

}

#endif
