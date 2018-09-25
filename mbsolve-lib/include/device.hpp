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
    /**
     * Constructs region with a certain material in 1D.
     *
     * \param [in] name    Name of the region.
     * \param [in] mat     Material of the region.
     * \param [in] x_start Start position in x-direction.
     * \param [in] x_end   End position in x-direction.
     */
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
     * Gets region name.
     */
    const std::string& get_name() const
    {
        return m_name;
    }

    /**
     * Gets region length.
     */
    real get_length() const
    {
        /* assert negative length */
        /* TODO: implement asserts */

        return (m_x_end - m_x_start);
    }

    /**
     * Gets region start position.
     */
    real get_x_start() const { return m_x_start; }

    /**
     * Gets region end position.
     */
    real get_x_end() const { return m_x_end; }

    /**
     * Gets material.
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
    /**
     * Constructs device.
     *
     * \param [in] name Name of the device.
     */
    device(const std::string& name);

    device(const std::string& file, const std::vector<material *>& materials);

    ~device();

    /**
     * Adds new region to device.
     */
    void add_region(std::shared_ptr<region> reg);

    /**
     * Gets all regions of device.
     */
    const std::vector<std::shared_ptr<region> >& get_regions() const;

    /**
     * Gets IDs of used materials.
     */
    const std::set<std::string>& get_used_materials() const;

    /**
     * Gets device name.
     */
    const std::string& get_name() const;

    /**
     * Gets device length.
     */
    real get_length() const;

    /**
     * Gets the minimum relative permittivity value.
     */
    real get_minimum_permittivity() const;

};

}

#endif
