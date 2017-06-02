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

#ifndef MBSOLVE_RECORD_H
#define MBSOLVE_RECORD_H

#include <string>
#include <types.hpp>

namespace mbsolve {

enum record_type { electric, magnetic, density };

/**
 * Represents a request to store certain simulation results (optionally at a
 * given index or given interval).
 * \ingroup MBSOLVE_LIB
 */
class record
{
private:
    std::string m_name;

    record_type m_type;

    unsigned int m_col;
    unsigned int m_row;

    real m_position;

    real m_interval;

public:
    record(const std::string& name, real interval = 0.0,
           real position = -1.0) :
        m_name(name), m_interval(interval), m_position(position)
    {

        /* TODO: parse name to type and indices */
    }

    const std::string& get_name() const { return m_name; }

    record_type get_type() const { return m_type; }

    unsigned int get_col() const { return m_col; }

    unsigned int get_row() const { return m_row; }

    real get_position() const { return m_position; }

    real get_interval() const { return m_interval; }
};

}

#endif
