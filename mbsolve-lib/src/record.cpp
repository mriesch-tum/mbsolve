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

#include <stdexcept>
#include <mbsolve/lib/record.hpp>

namespace mbsolve {

record::record(const std::string& name, real interval, real position)
  : m_name(name), m_interval(interval), m_position(position)
{
    /* type of requested result? */
    switch (m_name[0]) {
    case 'e':
        m_type = type::electric;
        break;
    case 'h':
        m_type = type::magnetic;
        break;
    case 'd':
        m_type = type::density;
        break;
    case 'i':
        m_type = type::inversion;
        break;
    default:
        throw std::invalid_argument("Unknown record type");
        break;
    }

    /* TODO parse numbers to identify e.g. density matrix entries */

    /* complex quantity? only off-diagonal density matrix entries */
    m_is_complex = (m_type == type::density) && (m_row != m_col);
}

record::record(
    const std::string& name,
    type record_type,
    unsigned int row_idx,
    unsigned int col_idx,
    real interval,
    real position)
  : m_name(name), m_type(record_type), m_row(row_idx - 1), m_col(col_idx - 1),
    m_interval(interval), m_position(position)
{
    m_is_complex = (m_type == type::density) && (m_row != m_col);
}
}
