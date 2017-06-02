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

#ifndef MBSOLVE_COPY_LIST_ENTRY_INT_H
#define MBSOLVE_COPY_LIST_ENTRY_INT_H

#include <memory>
#include <mbsolve.hpp>

namespace mbsolve {

/**
 * This internal class stores all data required to copy result data from
 * raw arrays to the \ref result class. Optionally, an extra set of addresses
 * can be stored to separate collecting results and transferring them back
 * to a host device.
 * \ingroup MBSOLVE_LIB
 */
class copy_list_entry
{
private:
    std::shared_ptr<result> m_result;
    std::shared_ptr<const record> m_record;
    real *m_scratch_real;
    real *m_scratch_imag;
    real *m_real;
    real *m_imag;

    unsigned int m_rows;
    unsigned int m_cols;
    unsigned int m_interval;
    unsigned int m_position;

public:
    copy_list_entry(std::shared_ptr<const record> rec,
                    std::shared_ptr<const scenario> scen) :
        m_record(rec), m_scratch_real(NULL), m_scratch_imag(NULL),
        m_real(NULL), m_imag(NULL)
    {
        m_rows = scen->get_endtime()/rec->get_interval();
        m_interval = ceil(1.0 * scen->get_num_timesteps()/m_rows);

        if (rec->get_position() < 0.0) {
            /* copy complete grid */
            m_position = 0;
            m_cols = scen->get_num_gridpoints();
        } else {
            m_position = std::round(rec->get_position()/
                                    scen->get_gridpoint_size());
            m_cols = 1;
        }

	/* create result */
        m_result = std::make_shared<result>(rec->get_name(), m_cols, m_rows);
    }

    bool hasto_record(unsigned int timestep) const {
        return (timestep % m_interval) == 0;
    }

    bool is_complex() const {
        return ((m_scratch_imag) && (m_imag));
    }

    unsigned int get_interval() const { return m_interval; }

    unsigned int get_position() const { return m_position; }

    unsigned int get_cols() const { return m_cols; }

    unsigned int get_rows() const { return m_rows; }

    unsigned int get_size() const { return m_cols * m_rows; }

    std::shared_ptr<const record> get_record() const { return m_record; };

    std::shared_ptr<result> get_result() const { return m_result; }

    std::vector<real>::iterator
    get_result_real(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_result->get_data_real(timestep/m_interval, gridpoint);
    }

    std::vector<real>::iterator
    get_result_imag(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_result->get_data_imag(timestep/m_interval, gridpoint);
    }

    real *
    get_scratch_real(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_scratch_real + (timestep/m_interval) * m_cols + gridpoint;
    }

    real *
    get_scratch_imag(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_scratch_imag + (timestep/m_interval) * m_cols + gridpoint;
    }

    real *
    get_real(unsigned int gridpoint = 0) const {
        return m_real + gridpoint;
    }

    real *
    get_imag(unsigned int gridpoint = 0) const {
        return m_imag + gridpoint;
    }

    void set_real(real *real_data) { m_real = real_data; }

    void set_imag(real *imag_data) { m_imag = imag_data; }

    void set_scratch_real(real *real_data) { m_scratch_real = real_data; }

    void set_scratch_imag(real *imag_data) { m_scratch_imag = imag_data; }
};

}

#endif
