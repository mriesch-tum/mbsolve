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

#ifndef MBSOLVE_LIB_RECORD_H
#define MBSOLVE_LIB_RECORD_H

#include <string>
#include <mbsolve/lib/types.hpp>

namespace mbsolve {

/**
 * Represents a request to store certain simulation results (optionally at a
 * given index or given interval).
 * \ingroup MBSOLVE_LIB
 */
class record
{
public:
    enum type
    {
        electric,
        polar_dt,
        magnetic,
        density,
        inversion
    };

private:
    std::string m_name;

    type m_type;

    unsigned int m_col;
    unsigned int m_row;

    real m_position;

    real m_interval;

    bool m_is_complex;

public:
    /**
     * Constructs record by deducing the type from the given \p name.
     *
     * \param [in] name     Name of the record, e.g., "e" for electric field.
     * \param [in] interval Sampling interval in seconds. If set to 0.0,
                            the value at every time step is stored.
     * \param [in] position Sampling position in meter. If set to -1.0,
                            the complete grid is stored.
     */
    record(
        const std::string& name,
        real interval = 0.0,
        real position = -1.0);

    /**
     * Constructs record.
     *
     * \param [in] name        Name of the record.
     * \param [in] record_type Type of the record, e.g., "density".
     * \param [in] row_idx     Only for density matrix entries. One-based.
     * \param [in] col_idx     Only for density matrix entries. One-based.
     * \param [in] interval    Sampling interval in seconds. If set to 0.0,
                               the value at every time step is stored.
     * \param [in] position    Sampling position in meter. If set to -1.0,
                               the complete grid is stored.
    */
    record(
        const std::string& name,
        type record_type,
        unsigned int row_idx,
        unsigned int col_idx,
        real interval = 0.0,
        real position = -1.0);

    /**
     * Gets name of the record.
     */
    const std::string& get_name() const { return m_name; }

    /**
     * Gets type of the record.
     */
    type get_type() const { return m_type; }

    /**
     * Gets column index (zero-based).
     */
    unsigned int get_col() const { return m_col; }

    /**
     * Gets row index (zero-based).
     */
    unsigned int get_row() const { return m_row; }

    /**
     * Gets sampling position (in meter).
     */
    real get_position() const { return m_position; }

    /**
     * Gets sampling interval (in seconds).
     */
    real get_interval() const { return m_interval; }

    /**
     * Is result complex valued?
     */
    bool is_complex() const { return m_is_complex; }
};
}

#endif
