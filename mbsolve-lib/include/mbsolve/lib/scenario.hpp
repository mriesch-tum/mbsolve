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

#ifndef MBSOLVE_LIB_SCENARIO_H
#define MBSOLVE_LIB_SCENARIO_H

#include <memory>
#include <string>
#include <vector>
#include <mbsolve/lib/qm_description.hpp>
#include <mbsolve/lib/record.hpp>
#include <mbsolve/lib/source.hpp>

namespace mbsolve {

/**
 * Abstract base class that represents the initial conditions of the
 * density matrix, which potentially depend on space.
 */
class ic_density
{
public:
    /**
     * Returns initial density matrix for given position x.
     */
    virtual const qm_operator& initialize(real x) = 0;
};

/**
 * Represents initial conditions of the density matrix which are constant
 * with respect to space.
 */
class ic_density_const : public ic_density
{
private:
    qm_operator m_rho;

public:
    /**
     * Constructs constant initial conditions from given initial density
     * matrix \p rho.
     *
     * \param [in] rho   Initial density matrix.
     */
    explicit ic_density_const(const qm_operator& rho) : m_rho(rho)
    {
        /* TODO check whether rho_init is a valid matrix
         * possibly exceptions for zero matrix, indicating that random init
         * is desired
         */
    }

    /**
     * Returns initial density matrix (independent of position).
     */
    const qm_operator& initialize(real /* x */) { return m_rho; }
};

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

    /* initial conditions for density matrix */
    std::shared_ptr<ic_density> m_dens_init;

public:
    /**
     * Constructs scenario.
     *
     * \param [in] name           Name of the scenario.
     * \param [in] num_gridpoints Number of spatial grid points.
     * \param [in] endtime        Simulation end time.
     * \param [in] rho_init       Initial density matrix.
     */
    scenario(
        const std::string& name,
        unsigned int num_gridpoints,
        real endtime,
        const qm_operator& rho_init);

    /**
     * Constructs scenario.
     *
     * \param [in] name           Name of the scenario.
     * \param [in] num_gridpoints Number of spatial grid points.
     * \param [in] endtime        Simulation end time.
     * \param [in] density_init   Initial conditions for density matrix.
     */
    scenario(
        const std::string& name,
        unsigned int num_gridpoints,
        real endtime,
        std::shared_ptr<ic_density> density_init);

    /**
     * Adds a record that specifies which data trace is collected.
     */
    void add_record(std::shared_ptr<record> rec);

    /**
     * Gets all records of this scenario.
     */
    const std::vector<std::shared_ptr<record> >& get_records() const;

    /**
     * Adds a source to the scenario.
     */
    void add_source(std::shared_ptr<source> src);

    /**
     * Gets all sources of this scenario.
     */
    const std::vector<std::shared_ptr<source> >& get_sources() const;

    /**
     * Returns name of the scenario.
     */
    const std::string& get_name() const;

    /**
     * Gets number of time steps.
     */
    unsigned int get_num_timesteps() const;

    /**
     * Sets the number of time steps manually. Use with care!
     */
    void set_num_timesteps(unsigned int value);

    /**
     * Gets number of spatial grid points.
     */
    unsigned int get_num_gridpoints() const;

    /**
     * Sets the number of grid points manually. Use with care!
     */
    void set_num_gridpoints(unsigned int value);

    /**
     * Gets time step size.
     */
    real get_timestep_size() const;

    /**
     * Sets the time step size manually. Use with care!
     */
    void set_timestep_size(real value);

    /**
     * Gets grid point size.
     */
    real get_gridpoint_size() const;

    /**
     * Sets the grid point size manually. Use with care!
     */
    void set_gridpoint_size(real value);

    /**
     * Gets the simulation end time.
     */
    real get_endtime() const;

    /**
     * Sets the simulation end time manually. Use with care!
     */
    void set_endtime(real value);

    /**
     * Gets initial conditions for density matrix.
     */
    std::shared_ptr<ic_density> get_ic_density() const;

    /**
     * Sets initial conditions for density matrix.
     */
    void set_ic_density(std::shared_ptr<ic_density> density_init);
};
}

#endif
