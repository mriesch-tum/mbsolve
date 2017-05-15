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

#ifndef MBSOLVE_SOLVER_INT_H
#define MBSOLVE_SOLVER_INT_H

#include <string>
#include <memory>
#include <vector>
#include <device.hpp>
#include <scenario.hpp>
#include <result.hpp>
#include <types.hpp>

namespace mbsolve {

class solver_int;

/**
 * Abstract solver factory class.
 * \ingroup MBSOLVE_LIB
 */
class solver_factory_int
{
public:
    solver_factory_int() { }

    virtual ~solver_factory_int() { }

    virtual std::shared_ptr<solver_int>
    create_instance(std::shared_ptr<const device> dev,
                    std::shared_ptr<scenario> scen) const = 0;

    virtual const std::string& get_name() const = 0;

};

/**
 * This internal class provides the base class for the different solver
 * implementations and collects the corresponding factories in a static array.
 * \ingroup MBSOLVE_LIB
 */
class solver_int
{
private:
    static std::map<std::string, solver_factory_int *> m_factories;

protected:
    std::shared_ptr<const device> m_device;

    std::shared_ptr<scenario> m_scenario;

    std::vector<std::shared_ptr<result> > m_results;

public:
    solver_int(std::shared_ptr<const device> dev,
               std::shared_ptr<scenario> scen);

    virtual ~solver_int();

    const scenario& get_scenario() const { return *m_scenario; }

    const device& get_device() const { return *m_device; }

    virtual std::string get_name() const = 0;

    virtual void run() const = 0;

    const std::vector<std::shared_ptr<result> >& get_results() const
    {
	return m_results;
    }

    static void register_factory(const std::string& name,
                                 solver_factory_int *factory);

    static solver_factory_int *find_factory(const std::string& name);

};

/*
 * Solver factory template. Every solver implementation T has to provide
 * a \ref solver_factory<T> to create an instance of the solver class
 * \ingroup MBSOLVE_LIB
 */
template<typename T>
class solver_factory : public solver_factory_int
{
private:
    std::string m_name;
public:
    explicit solver_factory(const std::string& name) : m_name(name) {
        solver_int::register_factory(name, this);
    }

    std::shared_ptr<solver_int>
    create_instance(std::shared_ptr<const device> dev,
                    std::shared_ptr<scenario> scen) const {
	return std::make_shared<T>(dev, scen);
    }

    const std::string& get_name() const { return m_name; }
};

}

#endif
