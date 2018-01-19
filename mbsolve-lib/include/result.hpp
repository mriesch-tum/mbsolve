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
#include <vector>
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

    std::vector<real> m_real;
    std::vector<real> m_imag;

    bool m_is_complex;

public:
    explicit result(const std::string& name, unsigned int cols,
		    unsigned int rows, bool is_complex = false) :
	m_name(name), m_cols(cols), m_rows(rows), m_count(cols * rows),
        m_real(cols * rows), m_imag(cols * rows), m_is_complex(is_complex)
    {
    }

    ~result() {
    }

    std::string get_name() const { return m_name; }

    unsigned int get_count() const { return m_count; }

    unsigned int get_cols() const { return m_cols; }

    unsigned int get_rows() const { return m_rows; }

    bool is_complex() const { return m_is_complex; }

    //   complex *get_data(unsigned int row = 0) { return &m_values[row * m_cols]; }

    std::vector<real>::iterator get_data_real(unsigned int row,
                                              unsigned int col) {
        return m_real.begin() + row * m_cols + col;
    }

    std::vector<real>::iterator get_data_imag(unsigned int row,
                                              unsigned int col) {
        return m_imag.begin() + row * m_cols + col;
    }

    real *get_data_real_raw() { return m_real.data(); }

    real *get_data_imag_raw() { return m_imag.data(); }

    const std::vector<real>& get_data_real() const { return m_real; }

    const std::vector<real>& get_data_imag() const { return m_imag; }

    std::vector<complex> get_data_complex() const;

};

}

#endif
