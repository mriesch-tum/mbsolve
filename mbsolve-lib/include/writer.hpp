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

#ifndef MBSOLVE_WRITER_H
#define MBSOLVE_WRITER_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <device.hpp>
#include <result.hpp>
#include <scenario.hpp>
#include <types.hpp>

namespace mbsolve {

/**
 * This class provides the static interface to create an instance of a writer
 * implementation and the base class fo each writer implementation.
 * \ingroup MBSOLVE_LIB
 */
class writer
{
public:
    typedef std::function<std::shared_ptr<writer>()> bootstrap_t;

private:
    static std::map<std::string, bootstrap_t> m_bootstraps;

protected:
    /**
     * Constructs writer (only available for derived classes).
     *
     * \param [in] name      Name of the writer method.
     * \param [in] extension File extension used by writer.
     */
    writer(const std::string& name, const std::string& extension);

    /* writer name, set during registration */
    std::string m_name;

    /* file extension, set by derived class */
    std::string m_ext;

public:
    /**
     * Destructs writer.
     */
    virtual ~writer();

    /**
     * Gets available writers.
     */
    static std::vector<std::string> get_avail_writers();

    /**
     * Registers writer bootstrap.
     */
    static void register_bootstrap(const std::string& name, bootstrap_t b);

    /**
     * Provides a shortcut function to register a writer of given type T.
     */
    template<typename T>
    static void register_writer(const std::string& name)
    {
        register_bootstrap(name, []() { return std::make_shared<T>(); });
    }

    /**
     * Constructs writer with a given \p name.
     *
     * \param [in] name Name of the writer method.
     */
    static std::shared_ptr<writer> create_instance(const std::string& name);

    /**
     * Gets writer name.
     */
    const std::string& get_name() const;

    /**
     * Gets file extension.
     */
    const std::string& get_extension() const;

    /**
     * Writes results to a \p file.
     *
     * \param [in] file     Filename.
     * \param [in] results  Results to be written.
     * \param [in] dev      Device that was simulated.
     * \param [in] scenario Scenario that was used.
     */
    virtual void write(
        const std::string& file,
        const std::vector<std::shared_ptr<result> >& results,
        std::shared_ptr<const device> dev,
        std::shared_ptr<const scenario> scen) const = 0;
};
}

#endif
