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

#ifndef MBSOLVE_RESULT_H
#define MBSOLVE_RESULT_H

#include <stdexcept>
#include <valarray>
#include <types.hpp>

namespace mbsolve {

/**
 * Represents a simulation result trace requested with a corresponding
 * \ref Record.
 *
 * \ingroup MBSOLVE_LIB
 */
class result
{
private:
    std::string m_name;
    unsigned int m_cols;
    unsigned int m_rows;
    unsigned int m_count;
    std::valarray<complex> m_values;

public:
    explicit result(const std::string& name, unsigned int cols,
		    unsigned int rows) :
	m_name(name), m_cols(cols), m_rows(rows), m_count(cols * rows),
        m_values(cols * rows)
    {
    }

    ~result() {
    }

    std::string get_name() const { return m_name; }

    unsigned int count() const { return m_count; }

    unsigned int cols() const { return m_cols; }

    unsigned int rows() const { return m_rows; }

    complex *data(unsigned int row = 0) { return &m_values[row * m_cols]; }

};

}

#endif
