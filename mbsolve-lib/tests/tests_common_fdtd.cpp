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

#include <catch2/catch.hpp>

#include <mbsolve/lib/device.hpp>
#include <mbsolve/lib/internal/common_fdtd.hpp>
#include <mbsolve/lib/scenario.hpp>

TEST_CASE("FDTD params -- Zero device length, zero end time", "[FDTD]")
{
    auto dev = std::make_shared<mbsolve::device>("dut");

    REQUIRE(dev->get_length() == 0.0);

    auto scen = std::make_shared<mbsolve::scenario>(
        "test", 1, 0, mbsolve::qm_operator({ 1, 0 }));

    REQUIRE_NOTHROW(init_fdtd_simulation(dev, scen));

    REQUIRE(scen->get_endtime() == 0.0);
    REQUIRE(scen->get_num_timesteps() == 1);
    REQUIRE(scen->get_timestep_size() == 3.154e10);
    REQUIRE(scen->get_num_gridpoints() == 1);
    REQUIRE(scen->get_gridpoint_size() == 1e13);
}

TEST_CASE("FDTD params -- Zero device length, nonzero end time", "[FDTD]")
{
    auto dev = std::make_shared<mbsolve::device>("dut");

    REQUIRE(dev->get_length() == 0.0);

    mbsolve::real t_e = 1e-9;
    unsigned int n_t = 10000;

    auto scen = std::make_shared<mbsolve::scenario>(
        "test", 1, t_e, mbsolve::qm_operator({ 1, 0 }), n_t);

    REQUIRE_NOTHROW(init_fdtd_simulation(dev, scen));

    REQUIRE(scen->get_endtime() == t_e);
    REQUIRE(scen->get_num_timesteps() == n_t);
    REQUIRE(scen->get_timestep_size() == Approx(t_e / (n_t - 1)));
    REQUIRE(scen->get_num_gridpoints() == 1);
    REQUIRE(scen->get_gridpoint_size() == 1e13);
}

TEST_CASE("FDTD params -- Nonzero device length, zero end time", "[FDTD]")
{
    mbsolve::real L = 1e-3;

    auto mat = std::make_shared<mbsolve::material>("Vacuum");

    auto dev = std::make_shared<mbsolve::device>("dut");
    dev->add_region(std::make_shared<mbsolve::region>("rut", mat, 0, L));

    REQUIRE(dev->get_length() == L);

    unsigned int n_x = 32768;

    auto scen = std::make_shared<mbsolve::scenario>(
        "test", n_x, 0.0, mbsolve::qm_operator({ 1, 0 }), 12345);

    REQUIRE_NOTHROW(init_fdtd_simulation(dev, scen));

    REQUIRE(scen->get_endtime() == 0.0);
    REQUIRE(scen->get_num_timesteps() == 1);
    REQUIRE(scen->get_timestep_size() == 3.154e10);
    REQUIRE(scen->get_num_gridpoints() == n_x);
    REQUIRE(scen->get_gridpoint_size() == Approx(L / (n_x - 1)));
}

TEST_CASE("FDTD params -- Nonzero device length, nonzero end time", "[FDTD]")
{
    mbsolve::real L = 1e-3;

    auto mat = std::make_shared<mbsolve::material>("Vacuum");

    auto dev = std::make_shared<mbsolve::device>("dut");
    dev->add_region(std::make_shared<mbsolve::region>("rut", mat, 0, L));

    REQUIRE(dev->get_length() == L);

    mbsolve::real t_e = 1e-9;
    unsigned int n_x = 32768;

    auto scen = std::make_shared<mbsolve::scenario>(
        "test", n_x, t_e, mbsolve::qm_operator({ 1, 0 }), 12345);

    REQUIRE_NOTHROW(init_fdtd_simulation(dev, scen));

    REQUIRE(scen->get_endtime() == t_e);
    REQUIRE(scen->get_num_gridpoints() == n_x);

    mbsolve::real d_x = L / (n_x - 1);
    REQUIRE(scen->get_gridpoint_size() == Approx(d_x));

    mbsolve::real velocity = 1.0 / sqrt(mbsolve::MU0 * mbsolve::EPS0);
    mbsolve::real d_t = d_x / velocity;
    unsigned int n_t = ceil(t_e / d_t) + 1;
    d_t = t_e / (n_t - 1);

    REQUIRE(scen->get_num_timesteps() == n_t);
    REQUIRE(scen->get_timestep_size() == Approx(d_t));
}
