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

#ifndef MBSOLVE_DEVICE_H
#define MBSOLVE_DEVICE_H

#include <set>
#include <string>
#include <memory>
#include <vector>
#include <material.hpp>
#include <types.hpp>

namespace mbsolve {

/**
 * Represents a section in the device or setup, has a certain \ref material.
 * \ingroup MBSOLVE_LIB
 */
class region
{
private:
    /* region name */
    std::string m_name;

    /* region material */
    std::shared_ptr<material> m_mat;

    /* dimensions */
    real m_x_start;
    real m_x_end;

public:
    region(const std::string& name,
           std::shared_ptr<material> mat,
           real x_start,
           real x_end) :
        m_name(name),
        m_mat(mat),
        m_x_start(x_start),
        m_x_end(x_end)
    {
    }

    /**
     * Get region name.
     */
    const std::string& get_name() const
    {
        return m_name;
    }

    /**
     * Get region length.
     */
    real get_length() const
    {
        /* assert negative length */
        /* TODO: implement asserts */

        return (m_x_end - m_x_start);
    }

    /**
     * Get region start position.
     */
    real get_start() const { return m_x_start; }

    /**
     * Get region end position.
     */
    real get_end() const { return m_x_end; }

    /**
     * Get material.
     */
    std::shared_ptr<material> get_material() const { return m_mat; }

};

/**
 * Represents a certain device or setup, consists of one or more \ref region.
 * \ingroup MBSOLVE_LIB
 */
class device
{
private:
    std::string m_name;

    std::vector<std::shared_ptr<region> > m_regions;

    std::set<std::string> m_used_materials;

    /* TODO: boundary conditions for fields */
    /* choices: periodic (ring cavity), PML, PMC (Fabry-Perot cavity) ... */

    /* TODO: boundary conditions for density matrix ? */
    /* choices: injection current. rather put it in scenario? as source? */
    /* choices: periodic */

public:

    device(const std::string& name);

    device(const std::string& file, const std::vector<material *>& materials);

    ~device();

    /**
     * Add new region to device.
     */
    void add_region(std::shared_ptr<region> reg);

    /**
     * Get all regions of device.
     */
    const std::vector<std::shared_ptr<region> >& get_regions() const;

    /**
     * Get IDs of used materials.
     */
    const std::set<std::string>& get_used_materials() const;

    /**
     * Get device name.
     */
    const std::string& get_name() const;

    /**
     * Get device length.
     */
    real get_length() const;

    /**
     * Get the minimum relative permittivity value.
     */
    real get_minimum_permittivity() const;

};

}

#endif
