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

#ifndef MBSOLVE_QM_DESCRIPTION_H
#define MBSOLVE_QM_DESCRIPTION_H

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
 * Quantum mechanical description of a n-level system, where n must be known
 * at compile time.
 * \ingroup MBSOLVE_LIB
 */
template<unsigned int n_lvl>
class qm_desc_clvl : public qm_description
{
    typedef Eigen::Matrix<complex, n_lvl, n_lvl> matrix_t;
    typedef matrix_t (*callback_t)(const matrix_t&);
private:
    /* Hamiltonian */
    matrix_t m_h;

    /* dipole moment operator */
    matrix_t m_u;

    /* Lindblad superoperator */
    callback_t m_g;

    /* initial density matrix */
    matrix_t m_d_init;

public:

    explicit qm_desc_clvl(real carrier_density,
                          const matrix_t& hamiltonian,
                          const matrix_t& dipole_op,
                          const callback_t& lindblad_op,
                          const matrix_t& d_init) :
        qm_description(carrier_density),
        m_h(hamiltonian), m_u(dipole_op), m_g(lindblad_op), m_d_init(d_init)
    {
    }

    const matrix_t& get_hamiltonian() const { return m_h; }

    const matrix_t& get_dipole_op() const { return m_u; }

    const callback_t& get_lindblad_op() const { return m_g; }

    const matrix_t& get_d_init() const { return m_d_init; }
};

typedef qm_desc_clvl<3> qm_desc_3lvl;

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

    /* Lindblad superoperator */
    /* TODO: functor */


public:
    explicit qm_desc_nlvl(real carrier_density) :
        qm_description(carrier_density) { }

    /* TODO: either this way, or make this the base class */
    /* TODO: base class for all descriptions or only the _Nlvl ones? */
    /* TODO: level count as template argument? */
    //qm_desc_nlvl(const qm_desc_2lvl& desc) { }



};

}

#endif
