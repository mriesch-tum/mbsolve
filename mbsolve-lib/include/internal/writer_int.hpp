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

#ifndef MBSOLVE_WRITER_INT_H
#define MBSOLVE_WRITER_INT_H

#include <map>
#include <string>
#include <memory>
#include <vector>
#include <device.hpp>
#include <scenario.hpp>
#include <result.hpp>
#include <types.hpp>

namespace mbsolve {

class writer_int;

/**
 * Abstract writer factory class.
 * \ingroup MBSOLVE_LIB_INT
 */
class writer_factory_int
{
public:
    writer_factory_int() { }

    virtual ~writer_factory_int() { }

    virtual std::shared_ptr<writer_int> create_instance() const = 0;

    virtual const std::string& get_name() const = 0;
};

/**
 * This internal class provides the base class for the different writer
 * implementations and collects the corresponding factories in a static array.
 * \ingroup MBSOLVE_LIB_INT
 */
class writer_int
{
private:
    static std::map<std::string, writer_factory_int *> m_factories;

public:
    writer_int();

    virtual ~writer_int();

    virtual const std::string& get_name() const = 0;

    virtual void write(const std::string& file,
                       const std::vector<std::shared_ptr<result> >& results,
                       std::shared_ptr<const device> dev,
                       std::shared_ptr<const scenario> scen) const = 0;

    virtual const std::string& get_extension() const = 0;

    static void register_factory(const std::string& name,
                                 writer_factory_int *factory);

    static writer_factory_int *find_factory(const std::string& name);
};

/*
 * Writer factory template. Every writer implementation T has to provide
 * a \ref writer_factory<T> to create an instance of the writer class.
 * \ingroup MBSOLVE_LIB_INT
 */
template<typename T>
class writer_factory : public writer_factory_int
{
private:
    std::string m_name;
public:
    explicit writer_factory(const std::string& name) : m_name(name) {
        writer_int::register_factory(name, this);
    }

    std::shared_ptr<writer_int> create_instance() const {
	return std::make_shared<T>();
    }

    const std::string& get_name() const { return m_name; }
};

}

#endif
