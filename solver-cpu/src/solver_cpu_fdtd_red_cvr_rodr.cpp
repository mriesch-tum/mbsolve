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

#include <iostream>
#include <mbsolve/lib/internal/algo_lindblad_cvr_rodr.hpp>
#include <mbsolve/solver-cpu/solver_cpu_fdtd_red.hpp>

#include <mbsolve/solver-cpu/solver_cpu_fdtd_red.tpp>

namespace mbsolve {

template class solver_cpu_fdtd_red<2, lindblad_cvr_rodr>;
template class solver_cpu_fdtd_red<3, lindblad_cvr_rodr>;
template class solver_cpu_fdtd_red<6, lindblad_cvr_rodr>;
}
