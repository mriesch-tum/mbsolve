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

#include <mbsolve/lib/scenario.hpp>

namespace mbsolve {

scenario::scenario(
    const std::string& name,
    unsigned int num_gridpoints,
    real endtime,
    const qm_operator& rho_init,
    unsigned int num_timesteps)
  : scenario(
        name,
        num_gridpoints,
        endtime,
        std::make_shared<ic_density_const>(rho_init),
        std::make_shared<ic_field_random>(),
        std::make_shared<ic_field_const>(0),
        num_timesteps)
{}

scenario::scenario(
    const std::string& name,
    unsigned int num_gridpoints,
    real endtime,
    std::shared_ptr<ic_density> density_init,
    std::shared_ptr<ic_field> electric_init,
    std::shared_ptr<ic_field> magnetic_init,
    unsigned int num_timesteps)
  : m_name(name), m_num_gridpoints(num_gridpoints), m_endtime(endtime),
    m_dens_init(density_init), m_e_init(electric_init),
    m_h_init(magnetic_init), m_num_timesteps(num_timesteps)
{}

void
scenario::add_record(std::shared_ptr<record> rec)
{
    m_records.push_back(rec);
}

const std::vector<std::shared_ptr<record> >&
scenario::get_records() const
{
    return m_records;
}

void
scenario::add_source(std::shared_ptr<source> src)
{
    m_sources.push_back(src);
}

const std::vector<std::shared_ptr<source> >&
scenario::get_sources() const
{
    return m_sources;
}

const std::string&
scenario::get_name() const
{
    return m_name;
}

unsigned int
scenario::get_num_timesteps() const
{
    return m_num_timesteps;
}

void
scenario::set_num_timesteps(unsigned int value)
{
    m_num_timesteps = value;
}

unsigned int
scenario::get_num_gridpoints() const
{
    return m_num_gridpoints;
}

void
scenario::set_num_gridpoints(unsigned int value)
{
    m_num_gridpoints = value;
}

real
scenario::get_timestep_size() const
{
    return m_timestep_size;
}

void
scenario::set_timestep_size(real value)
{
    m_timestep_size = value;
}

real
scenario::get_gridpoint_size() const
{
    return m_gridpoint_size;
}

void
scenario::set_gridpoint_size(real value)
{
    m_gridpoint_size = value;
}

real
scenario::get_endtime() const
{
    return m_endtime;
}

void
scenario::set_endtime(real value)
{
    m_endtime = value;
}

std::shared_ptr<ic_density>
scenario::get_ic_density() const
{
    return m_dens_init;
}

void
scenario::set_ic_density(std::shared_ptr<ic_density> density_init)
{
    m_dens_init = density_init;
}

std::shared_ptr<ic_field>
scenario::get_ic_electric() const
{
    return m_e_init;
}

void
scenario::set_ic_electric(std::shared_ptr<ic_field> electric_init)
{
    m_e_init = electric_init;
}

std::shared_ptr<ic_field>
scenario::get_ic_magnetic() const
{
    return m_h_init;
}

void
scenario::set_ic_magnetic(std::shared_ptr<ic_field> magnetic_init)
{
    m_h_init = magnetic_init;
}
}
