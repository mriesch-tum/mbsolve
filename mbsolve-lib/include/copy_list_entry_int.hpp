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
//#include <mbsolve.hpp>
#include <record.hpp>
#include <result.hpp>
#include <scenario.hpp>
#include <types.hpp>

/* TODO: needs improvement */
#ifdef XEON_PHI_OFFLOAD
#define __mb_on_device __attribute__((target(mic)))
#else
#ifdef __CUDACC__
#define __mb_on_device __host__ __device__
#include <math.h>
#else
#define __mb_on_device
#endif
#endif

namespace mbsolve {

class copy_list_entry_dev
{

    /*TODO make members private -> friend/nested with copy_list_entry? */

    /* TODO make precompiler variable __attribute etc .. ? */

public:
    unsigned int m_cols;
    unsigned int m_interval_idx;
    unsigned int m_position_idx;

    real m_timestep;
    real m_interval;

    unsigned int m_scratch_offset;

    record::type m_type;


    __mb_on_device bool hasto_record(unsigned int iteration) const {
        real t_now = iteration * m_timestep;
        real t_next = t_now + m_timestep;
        real t_sample = floor(1.0 * iteration / m_interval_idx) * m_interval;

        return ((t_now <= t_sample) && (t_next >= t_sample));
    }

    bool is_complex() const {
        //return ((m_scratch_imag) && (m_imag));
    }

    unsigned int get_interval() const { return m_interval; }


    __mb_on_device record::type get_type() const {
        return m_type;
    }

    __mb_on_device unsigned int get_position() const {
        return m_position_idx;
    }

    __mb_on_device unsigned int get_cols() const {
        return m_cols;
    }
    /*
    __attribute__((target(mic))) real *
    get_real(unsigned int gridpoint = 0, unsigned int thread_id = 0) const {
        //  if (m_use_threaded) {
            return m_real_threaded[thread_id] + gridpoint;
            //} else {
            //    return m_real + gridpoint;
            //}
            }*/

    __mb_on_device unsigned int
    get_scratch_real_offset(unsigned int timestep,
                            unsigned int gridpoint = 0) const {
        return m_scratch_offset + (timestep/m_interval_idx) * m_cols
            + gridpoint;
    }

    /*    __attribute__((target(mic))) real *
    get_scratch_real(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_scratch_real + (timestep/interval_idx) * cols + gridpoint;
    }
    /*
    real *
    get_scratch_imag(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_scratch_imag + (timestep/interval_idx) * cols + gridpoint;
    }


    real *
    get_imag(unsigned int gridpoint = 0, unsigned int thread_id = 0) const {
        if (m_use_threaded) {
            return m_imag_threaded[thread_id] + gridpoint;
        } else {
            return m_imag + gridpoint;
        }
    }

    void set_real(real **real_data) {
        m_real_threaded = real_data;
        m_use_threaded = true;
    }

    void set_real(real *real_data) {
        m_real = real_data;
        m_use_threaded = false;
    }

    void set_imag(real **imag_data) {
        m_imag_threaded = imag_data;
        m_use_threaded = true;
    }

    void set_imag(real *imag_data) {
        m_imag = imag_data;
        m_use_threaded = false;
    }

    void set_scratch_real(real *real_data) { m_scratch_real = real_data; }

    void set_scratch_imag(real *imag_data) { m_scratch_imag = imag_data; }

    */
};


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
    real **m_real_threaded;
    real **m_imag_threaded;

    bool m_use_threaded;

    unsigned int m_rows;
    unsigned int m_cols;
    unsigned int m_interval_idx;
    unsigned int m_position_idx;

    real m_timestep;
    real m_interval;



public:
    copy_list_entry(std::shared_ptr<const record> rec,
                    std::shared_ptr<const scenario> scen) :
        m_record(rec), m_scratch_real(NULL), m_scratch_imag(NULL),
        m_real(NULL), m_imag(NULL), m_use_threaded(false)
    {
        m_timestep = scen->get_timestep_size();

        m_rows = ceil(scen->get_endtime()/rec->get_interval()) + 1;
        m_interval = scen->get_endtime()/(m_rows - 1);

        m_interval_idx =
            floor(1.0 * (scen->get_num_timesteps() - 1)/(m_rows - 1));

        if (rec->get_position() < 0.0) {
            /* copy complete grid */
            m_position_idx = 0;
            m_cols = scen->get_num_gridpoints();
        } else {
            m_position_idx = std::round(rec->get_position()/
                                        scen->get_gridpoint_size());
            m_cols = 1;
        }

	/* create result */
        m_result = std::make_shared<result>(rec->get_name(), m_cols, m_rows);
        /*
        m_dev = cle_dev(m_rows, m_cols, m_interval_idx, m_position,
        scen->get_timestep_size(), rec->get_interval());*/


        m_dev.m_interval = rec->get_interval();
        m_dev.m_interval_idx = m_interval_idx;
        m_dev.m_position_idx = m_position_idx;
        m_dev.m_timestep = scen->get_timestep_size();
        m_dev.m_cols = m_cols;

    }

    copy_list_entry_dev m_dev;

    const copy_list_entry_dev& get_dev() const { return m_dev; }

    record::type get_type() const { return m_dev.get_type(); }

    bool hasto_record(unsigned int iteration) const {
        real t_now = iteration * m_timestep;
        real t_next = t_now + m_timestep;
        real t_sample = floor(iteration / m_interval_idx) * m_interval;

        return ((t_now <= t_sample) && (t_next >= t_sample));
    }

    bool is_complex() const {
        return ((m_scratch_imag) && (m_imag));
    }

    unsigned int get_interval() const { return m_interval; }

    unsigned int get_position() const { return m_position_idx; }

    unsigned int get_cols() const { return m_cols; }

    unsigned int get_rows() const { return m_rows; }

    unsigned int get_size() const { return m_cols * m_rows; }

    std::shared_ptr<const record> get_record() const { return m_record; };

    std::shared_ptr<result> get_result() const { return m_result; }

    std::vector<real>::iterator
    get_result_real(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_result->get_data_real(timestep/m_interval_idx, gridpoint);
    }

    std::vector<real>::iterator
    get_result_imag(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_result->get_data_imag(timestep/m_interval_idx, gridpoint);
    }

    unsigned int
    get_scratch_real_offset(unsigned int timestep,
                            unsigned int gridpoint = 0) const {
        return m_dev.get_scratch_real_offset(timestep, gridpoint);
    }

    real *
    get_scratch_real(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_scratch_real + (timestep/m_interval_idx) * m_cols + gridpoint;
    }

    real *
    get_scratch_imag(unsigned int timestep, unsigned int gridpoint = 0) const {
        return m_scratch_imag + (timestep/m_interval_idx) * m_cols + gridpoint;
    }


    real *
    get_real(unsigned int gridpoint = 0, unsigned int thread_id = 0) const {
        if (m_use_threaded) {
            return m_real_threaded[thread_id] + gridpoint;
        } else {
            return m_real + gridpoint;
        }
    }

    real *
    get_imag(unsigned int gridpoint = 0, unsigned int thread_id = 0) const {
        if (m_use_threaded) {
            return m_imag_threaded[thread_id] + gridpoint;
        } else {
            return m_imag + gridpoint;
        }
    }

    void set_real(real **real_data) {
        m_real_threaded = real_data;
        m_use_threaded = true;
    }

    void set_real(real *real_data) {
        m_real = real_data;
        m_use_threaded = false;
        }

    void set_imag(real **imag_data) {
        m_imag_threaded = imag_data;
        m_use_threaded = true;
    }

    void set_imag(real *imag_data) {
        m_imag = imag_data;
        m_use_threaded = false;
    }

    void set_scratch_real(real *real_data, unsigned int offset = 0) {
        m_scratch_real = real_data;
        m_dev.m_scratch_offset = offset;
    }

    void set_scratch_imag(real *imag_data) { m_scratch_imag = imag_data; }
};

}

#endif
