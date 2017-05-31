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

#ifndef MBSOLVE_SCENARIO_H
#define MBSOLVE_SCENARIO_H

#include <memory>
#include <string>
#include <vector>
#include <record.hpp>
#include <source.hpp>

namespace mbsolve {

/**
 * Stores simulation scenario (simulation settings as well as \ref source
 * objects and a collection of \ref record.
 * \ingroup MBSOLVE_LIB
 */
class scenario
{
private:
    std::string m_name;

    unsigned int m_num_timesteps;

    unsigned int m_num_gridpoints;

    real m_timestep_size;

    real m_gridpoint_size;

    real m_endtime;

    std::vector<std::shared_ptr<record> > m_records;

    std::vector<std::shared_ptr<source> > m_sources;

    /* TODO: initial conditions fields */
    /* choice: random, zero */

    /* TODO: initial conditions density matrix */
    /* choice: equilibrium (all in lowest level), random population */
    /* off-diagonal elements always zero? */

public:

    scenario(const std::string& name, unsigned int num_gridpoints,
             real endtime);

    void add_record(std::shared_ptr<record> rec);

    const std::vector<std::shared_ptr<record> >& get_records() const;

    void add_source(std::shared_ptr<source> src);

    const std::vector<std::shared_ptr<source> >& get_sources() const;

    const std::string& get_name() const;

    unsigned int get_num_timesteps() const;

    void set_num_timesteps(unsigned int value);

    unsigned int get_num_gridpoints() const;

    void set_num_gridpoints(unsigned int value);

    real get_timestep_size() const;

    void set_timestep_size(real value);

    real get_gridpoint_size() const;

    void set_gridpoint_size(real value);

    real get_endtime() const;

    void set_endtime(real value);

};

}

#endif
