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

#include <solver_generic.hpp>

namespace mbsolve {

static SolverFactory<solver_generic> factory("generic");

solver_generic::solver_generic(const Device& dev, const Scenario& scen) :
    ISolver(dev, scen)
{
}

solver_generic::~solver_generic()
{
}

std::string
solver_generic::getName() const
{
    return factory.getName();
}

void
solver_generic::run() const
{

}

}
