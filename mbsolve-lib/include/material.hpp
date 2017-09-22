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
#include <Eigen/Dense>
#include <types.hpp>

namespace mbsolve {

/**
 * Provides the quantum mechanical description of an active \ref region.
 * \ingroup MBSOLVE_LIB
 */
class qm_description
{
private:

    /* density of charge carriers */
    real m_carrier_density;

    /*
      real m_period_length;
    */

public:
    explicit qm_description(real carrier_density) :
        m_carrier_density(carrier_density)
    {
    }

    virtual ~qm_description()
    {
    }

    /*
      use some sparse structures that provide unique access to a given
      element. transition frequencies are dense, but coupling and anticrossing
      are likely to be sparse.
     */

    real get_carrier_density() const
    {
        return m_carrier_density;
    }

};

/**
 * Quantum mechanical description of a 2-level system.
 * \ingroup MBSOLVE_LIB
 */
class qm_desc_2lvl : public qm_description
{
private:
    /* transition frequency */
    real m_trans_freq;

    /* dipole moment */
    real m_dipole_mom;

    /* scattering rate from upper laser level to lower laser level */
    real m_scattering;

    /* dephasing rate */
    real m_dephasing;

    /* equilibrium population inversion */
    real m_equi_inv;

public:
    explicit qm_desc_2lvl(real carrier_density = 0.0,
                          real transition_freq = 0.0,
                          real dipole_moment = 0.0,
                          real scattering_rate = 0.0,
                          real dephasing_rate = 0.0,
                          real equilibrium_inversion = -1.0) :
        qm_description(carrier_density),
        m_trans_freq(transition_freq),
        m_dipole_mom(dipole_moment),
        m_scattering(scattering_rate),
        m_dephasing(dephasing_rate),
        m_equi_inv(equilibrium_inversion)
    {
    }

    ~qm_desc_2lvl()
    {
    }

    /**
     * Get transition frequency between upper and lower laser level.
     */
    real get_transition_freq() const { return m_trans_freq; }

    /**
     * Get dipole moment between upper and lower laser level.
     */
    real get_dipole_moment() const { return m_dipole_mom; }

    /**
     * Get scattering rate between upper and lower laser level.
     */
    real get_scattering_rate() const { return m_scattering; }

    /**
     * Get dephasing rate.
     */
    real get_dephasing_rate() const { return m_dephasing; }

    /**
     * Get equilibrium population inversion.
     */
    real get_equilibrium_inversion() const { return m_equi_inv; }
};

/**
 * Quantum mechanical description of a 3-level system.
 * \ingroup MBSOLVE_LIB
 */
class qm_desc_3lvl : public qm_description
{
private:
    /* N x N Hamiltonian */
    Eigen::Matrix<real, 3, 3> m_h;

    /* N x N dipole moment operator */
    Eigen::Matrix<real, 3, 3> m_u;

    /* N^2 x N^2 Lindblad superoperator */
    Eigen::Matrix<real, 9, 9> m_g;

public:

    explicit qm_desc_3lvl(real carrier_density,
                          const Eigen::Matrix<real, 3, 3>& hamiltonian,
                          const Eigen::Matrix<real, 3, 3>& dipole_op,
                          const Eigen::Matrix<real, 9, 9>& lindblad_op) :
        qm_description(carrier_density),
        m_h(hamiltonian), m_u(dipole_op), m_g(lindblad_op) { }


    const Eigen::Matrix<real, 3, 3>& get_hamiltonian() const { return m_h; }

    const Eigen::Matrix<real, 3, 3>& get_dipole_op() const { return m_u; }

    const Eigen::Matrix<real, 9, 9>& get_lindblad_op() const { return m_g; }
};

/**
 * Quantum mechanical description of a N-level system.
 * \ingroup MBSOLVE_LIB
 */
class qm_desc_nlvl : public qm_description
{
private:
    /* N x N Hamiltonian */
    /* TODO: matrix */

    /* N x N dipole moment operator */
    /* TODO: matrix */

    /* N^2 x N^2 Lindblad superoperator */
    /* TODO: matrix */


public:
    explicit qm_desc_nlvl(real carrier_density) :
        qm_description(carrier_density) { }

    /* TODO: either this way, or make this the base class */
    /* TODO: base class for all descriptions or only the _Nlvl ones? */
    /* TODO: level count as template argument? */
    //qm_desc_nlvl(const qm_desc_2lvl& desc) { }



};

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
