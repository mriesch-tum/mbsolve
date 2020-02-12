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

#include <iostream>
#include <solver_cpu_fdtd.hpp>
#include <solver_cpu_fdtd_red.hpp>
#include <solver_cpu_loader.hpp>
#include <mbsolve/lib/internal/algo_lindblad_cvr_rodr.hpp>
#include <mbsolve/lib/internal/algo_lindblad_noop.hpp>

namespace mbsolve {

solver_cpu_loader::solver_cpu_loader()
{
    /* fdtd noop */
    solver::register_solver<solver_cpu_fdtd<0, lindblad_noop> >(
        "cpu-fdtd-noop");

    /* fdtd-red noop */
    solver::register_solver<solver_cpu_fdtd_red<0, lindblad_noop> >(
        "cpu-fdtd-red-noop");

    /* fdtd coherence vector representation/rodrigues formula */
    solver::register_solver<solver_cpu_fdtd<2, lindblad_cvr_rodr> >(
        "cpu-fdtd-2lvl-cvr-rodr");
    solver::register_solver<solver_cpu_fdtd<3, lindblad_cvr_rodr> >(
        "cpu-fdtd-3lvl-cvr-rodr");
    solver::register_solver<solver_cpu_fdtd<6, lindblad_cvr_rodr> >(
        "cpu-fdtd-6lvl-cvr-rodr");

    /* fdtd-red coherence vector representation/rodrigues formula */
    solver::register_solver<solver_cpu_fdtd_red<2, lindblad_cvr_rodr> >(
        "cpu-fdtd-red-2lvl-cvr-rodr");
    solver::register_solver<solver_cpu_fdtd_red<3, lindblad_cvr_rodr> >(
        "cpu-fdtd-red-3lvl-cvr-rodr");
    solver::register_solver<solver_cpu_fdtd_red<6, lindblad_cvr_rodr> >(
        "cpu-fdtd-red-6lvl-cvr-rodr");
}
}
