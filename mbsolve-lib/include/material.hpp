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

#ifndef MBSOLVE_MATERIAL_H
#define MBSOLVE_MATERIAL_H

#include <map>
#include <memory>
#include <string>
#include <qm_description.hpp>
#include <types.hpp>

namespace mbsolve {

/**
 * The material contains electromagnetic properties and -- if specified -- a
 * quantum mechanical description.
 * \ingroup MBSOLVE_LIB
 */
class material
{
private:
    /* material id */
    std::string m_id;

    /* electromagnetic properties */
    real m_rel_permittivity;
    real m_rel_permeability;
    real m_overlap_factor;
    real m_losses;

    /* quantum mechanical description of active material */
    std::shared_ptr<qm_description> m_qm;

    /* material library */
    static std::map<std::string, std::shared_ptr<material> > m_materials;

public:

    explicit material(const std::string& id,
                      std::shared_ptr<qm_description> qm = nullptr,
                      real rel_permittivity = 1.0,
                      real overlap_factor = 1.0,
                      real losses = 0.0,
                      real rel_permeability = 1.0) :
        m_id(id),
        m_qm(qm),
        m_rel_permittivity(rel_permittivity),
        m_rel_permeability(rel_permeability),
        m_losses(losses),
        m_overlap_factor(overlap_factor)
    {
    }

    ~material()
    {
    }

    /* TODO copy constructor */

    /* TODO assignment operator */

    /**
     * Get the material ID.
     */
    const std::string& get_id() const { return m_id; }

    /**
     * Get pointer to quantum mechanical description.
     */
    std::shared_ptr<qm_description> get_qm() const { return m_qm; }

    /**
     * Get relative permittivity &epsilon;<sub>r</sub>
     */
    real get_rel_permittivity() const { return m_rel_permittivity; }

    /**
     * Get relative permeability &mu;<sub>r</sub>
     */
    real get_rel_permeability() const { return m_rel_permeability; }

    /**
     * Get losses &alpha;.
     */
    real get_losses() const { return m_losses; }

    /**
     * Get overlap factor &Gamma;.
     */
    real get_overlap_factor() const { return m_overlap_factor; }

    //    void add_to_library() const;

    static void add_to_library(std::shared_ptr<material> mat);

    static void add_to_library(const material& mat);

    static std::shared_ptr<material> get_from_library(const std::string& id);
};

}

#endif
