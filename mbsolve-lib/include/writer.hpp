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

#include <string>
#include <vector>
#include <device.hpp>
#include <internal/writer_int.hpp>
#include <scenario.hpp>
#include <types.hpp>

namespace mbsolve {

/**
 * This class provides the interface to create an instance of a writer
 * implementation. Each implementation is a subclass of \ref writer_int and
 * is created internally.
 * \ingroup MBSOLVE_LIB
 */
class writer
{
private:
    std::shared_ptr<writer_int> m_writer;

public:
    /**
     * Constructs writer of a given \p name.
     */
    writer(const std::string& name);

    ~writer();

    /**
     * Writes results to a \p file.
     *
     * \param [in] file     Filename.
     * \param [in] results  Results to be written.
     * \param [in] dev      Device that was simulated.
     * \param [in] scenario Scenario that was used.
     */
    void write(const std::string& file,
               const std::vector<std::shared_ptr<result> >& results,
               std::shared_ptr<const device> dev,
               std::shared_ptr<const scenario> scen) const;

    /**
     * Gets file extension.
     */
    const std::string& get_extension() const;
};


}

#endif
