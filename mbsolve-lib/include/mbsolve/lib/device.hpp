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

#ifndef MBSOLVE_LIB_DEVICE_H
#define MBSOLVE_LIB_DEVICE_H

#include <memory>
#include <set>
#include <string>
#include <vector>
#include <mbsolve/lib/material.hpp>
#include <mbsolve/lib/types.hpp>

namespace mbsolve {

/**
 * Abstract base class that represents the boundary conditions of the
 * electromagnetic field.
 * \ingroup MBSOLVE_LIB
 */
class bc_field
{
public:
    virtual ~bc_field() {}
};

/**
 * Represents boundary conditions for the electromagnetic field that yield
 * the specified reflectivity values at the ends of the device.
 * \ingroup MBSOLVE_LIB
 */
class bc_field_reflectivity : public bc_field
{
private:
    real m_refl_left;
    real m_refl_right;

public:
    /**
     * Constructs boundary conditions that yield a given reflectivity
     * values, which default to perfect reflection.
     *
     * \param [in] refl_left  Reflectivity at left end of the device.
     * \param [in] refl_right Reflectivity at right end of the device.
     */
    explicit bc_field_reflectivity(
        real refl_left = 1.0,
        real refl_right = 1.0)
      : m_refl_left(refl_left), m_refl_right(refl_right)
    {}

    /**
     * Gets right reflectivity value.
     */
    real get_refl_right() { return m_refl_right; }

    /**
     * Gets left reflectivity value.
     */
    real get_refl_left() { return m_refl_left; }
};

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
    region(
        const std::string& name,
        std::shared_ptr<material> mat,
        real x_start,
        real x_end)
      : m_name(name), m_mat(mat), m_x_start(x_start), m_x_end(x_end)
    {}

    /**
     * Gets region name.
     */
    const std::string& get_name() const { return m_name; }

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

    /* boundary conditions for electromagnetic field */
    std::shared_ptr<bc_field> m_bc_field;

public:
    /**
     * Constructs device.
     *
     * \param [in] name              Name of the device.
     * \param [in] field_boundary    Boundary conditions for the EM field.
     */
    device(
        const std::string& name,
        std::shared_ptr<bc_field> field_boundary =
            std::make_shared<bc_field_reflectivity>(1.0, 1.0));

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

    /**
     * Sets boundary condition for the field.
     */
    std::shared_ptr<bc_field> get_bc_field() const;

    /**
     * Gets boundary condition for the field.
     */
    void set_bc_field(std::shared_ptr<bc_field> boundary_field);
};
}

#endif
