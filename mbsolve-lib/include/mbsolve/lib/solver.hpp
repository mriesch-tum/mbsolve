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

#ifndef MBSOLVE_LIB_SOLVER_H
#define MBSOLVE_LIB_SOLVER_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <mbsolve/lib/device.hpp>
#include <mbsolve/lib/result.hpp>
#include <mbsolve/lib/scenario.hpp>
#include <mbsolve/lib/types.hpp>

namespace mbsolve {

/**
 * This class provides the static interface to create an instance of a solver
 * implementation and the base class fo each solver implementation.
 * \ingroup MBSOLVE_LIB
 */
class solver
{
public:
    typedef std::function<std::shared_ptr<solver>(
        std::shared_ptr<const device> dev,
        std::shared_ptr<scenario> scen)>
        bootstrap_t;

private:
    static std::map<std::string, bootstrap_t> m_bootstraps;

protected:
    /**
     * Constructs solver (only available for derived classes).
     *
     * \param [in] name      Name of the writer method.
     */
    solver(
        const std::string& name,
        std::shared_ptr<const device> dev,
        std::shared_ptr<scenario> scen);

    /* solver name, set during registration */
    std::string m_name;

    /* device to be simulated */
    std::shared_ptr<const device> m_device;

    /* simulation scenario */
    std::shared_ptr<scenario> m_scenario;

    /* simulation results */
    std::vector<std::shared_ptr<result> > m_results;

public:
    /**
     * Destructs solver.
     */
    virtual ~solver();

    /**
     * Gets available solvers.
     */
    static std::vector<std::string> get_avail_solvers();

    /**
     * Registers solver bootstrap.
     */
    static void register_bootstrap(const std::string& name, bootstrap_t b);

    /**
     * Provides a shortcut function to register a solver of given type T.
     */
    template<typename T>
    static void register_solver(const std::string& name)
    {
        register_bootstrap(
            name,
            [](std::shared_ptr<const device> dev,
               std::shared_ptr<scenario> scen) {
                return std::make_shared<T>(dev, scen);
            });
    }

    /**
     * Constructs solver with a given \p name.
     *
     * \param [in] name Name of the solver method.
     * \param [in] dev  Specify the \ref device to be simulated.
     * \param [in] scen Specify the \ref scenario.
     */
    static std::shared_ptr<solver> create_instance(
        const std::string& name,
        std::shared_ptr<const device> dev,
        std::shared_ptr<scenario> scen);

    /**
     * Gets solver name.
     */
    const std::string& get_name() const;

    /**
     * Gets scenario.
     */
    const scenario& get_scenario() const;

    /**
     * Gets device.
     */
    const device& get_device() const;

    /**
     * Executes solver (empty).
     */
    virtual void run() const = 0;

    /**
     * Gets results.
     */
    const std::vector<std::shared_ptr<result> >& get_results() const;
};
}

#endif
