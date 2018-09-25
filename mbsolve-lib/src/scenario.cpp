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

#include <scenario.hpp>

namespace mbsolve {

scenario::scenario(const std::string& name, unsigned int num_gridpoints,
                   real endtime, const qm_operator& rho_init) :
    m_name(name), m_num_gridpoints(num_gridpoints), m_endtime(endtime),
    m_rho_init(rho_init)
{
}

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


qm_operator
scenario::get_rho_init() const
{
    return m_rho_init;
}

void
scenario::set_rho_init(const qm_operator& rho_init)
{
    /* TODO check whether rho_init is a valid matrix */
    /* possibly exceptions for zero matrix, indicating that random init
     * is desired */

    m_rho_init = rho_init;
}

}
