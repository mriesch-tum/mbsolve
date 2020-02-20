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

#ifndef MBSOLVE_LIB_MATERIAL_H
#define MBSOLVE_LIB_MATERIAL_H

#include <map>
#include <memory>
#include <string>
#include <mbsolve/lib/qm_description.hpp>
#include <mbsolve/lib/types.hpp>

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
    /**
     * Constructs material from given parameters and quantum mechanical
     * description (if provided). Defaults to vacuum.
     *
     * \param [in] id               Identifier of material, must be unique.
     * \param [in] qm               Quantum mechanical description.
     * \param [in] rel_permittivity Relative permittivity.
     * \param [in] overlap_factor   Overlap of active region an electric field
                                    between 0.0 and 1.0.
     * \param [in] losses           Optical losses.
     * \param [in] rel_permeabiliy  Relative Permeability.
     */
    explicit material(
        const std::string& id,
        std::shared_ptr<qm_description> qm = nullptr,
        real rel_permittivity = 1.0,
        real overlap_factor = 1.0,
        real losses = 0.0,
        real rel_permeability = 1.0)
      : m_id(id), m_qm(qm), m_rel_permittivity(rel_permittivity),
        m_rel_permeability(rel_permeability), m_losses(losses),
        m_overlap_factor(overlap_factor)
    {}

    ~material() {}

    /* TODO copy constructor */

    /* TODO assignment operator */

    /**
     * Gets the material ID.
     */
    const std::string& get_id() const { return m_id; }

    /**
     * Gets pointer to quantum mechanical description.
     */
    std::shared_ptr<qm_description> get_qm() const { return m_qm; }

    /**
     * Gets relative permittivity &epsilon;<sub>r</sub>
     */
    real get_rel_permittivity() const { return m_rel_permittivity; }

    /**
     * Gets relative permeability &mu;<sub>r</sub>
     */
    real get_rel_permeability() const { return m_rel_permeability; }

    /**
     * Gets losses &alpha;.
     */
    real get_losses() const { return m_losses; }

    /**
     * Gets overlap factor &Gamma;.
     */
    real get_overlap_factor() const { return m_overlap_factor; }

    /**
     * Adds material to library.
     *
     * \param [in] mat Material to be added.
     */
    static void add_to_library(std::shared_ptr<material> mat);

    /**
     * Gets material from library.
     *
     * \param [in] id Material identifier.
     */
    static std::shared_ptr<material> get_from_library(const std::string& id);
};
}

#endif
