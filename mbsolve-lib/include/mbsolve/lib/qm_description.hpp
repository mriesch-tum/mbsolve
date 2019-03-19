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

#ifndef MBSOLVE_LIB_QM_DESCRIPTION_H
#define MBSOLVE_LIB_QM_DESCRIPTION_H

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <mbsolve/lib/types.hpp>

namespace mbsolve {

/**
 * Quantum operator in matrix form, e.g. for the Hamiltonian, dipole moment
 * operator, ...
 * \ingroup MBSOLVE_LIB
 */
class qm_operator
{
private:
    /* number of levels */
    unsigned int m_num_levels;

    /* main diagonal entries */
    std::vector<real> m_main_diag;

    /* off-diagonal entries */
    std::vector<std::complex<real> > m_off_diag;

public:
    /**
     * Constructs quantum mechanical operator using two vectors.
     */
    explicit qm_operator(
        const std::vector<real>& main_diagonal,
        /**< [in] Main diagonal entries (real) */
        const std::vector<std::complex<real> >& off_diagonal =
            std::vector<std::complex<real> >()
        /**< [in] Top half of the off-diagonal entries
         *   (complex), arranged in column-major ordering */
        )
      : m_num_levels(main_diagonal.size()), m_main_diag(main_diagonal),
        m_off_diag(off_diagonal)
    {
        /* TODO: assert or exception: vector sizes do not match */
    }

    virtual ~qm_operator() {}

    /**
     * Gets number of levels.
     */
    unsigned int get_num_levels() const { return m_num_levels; }

    /**
     * Gets main diagonal entries
     */
    const std::vector<real>& get_main_diagonal() const { return m_main_diag; }

    /**
     * Gets off-diagonal entries (top half, column major ordering)
     */
    const std::vector<std::complex<real> >& get_off_diagonal() const
    {
        return m_off_diag;
    }
};

/**
 * Quantum superoperator, e.g. for the relaxation superoperator.
 * \ingroup MBSOLVE_LIB
 */
class qm_superop
{
private:
protected:
    /* number of levels */
    unsigned int m_num_levels;

public:
    /**
     * Constructs quantum mechanical superoperator.
     *
     * \param [in] num_levels Number of energy levels.
     */
    explicit qm_superop(unsigned int num_levels) : m_num_levels(num_levels) {}

    virtual ~qm_superop() {}

    /**
     * Executes the superoperator on the operator \p arg. Alias to \ref
     * action.
     *
     * \param [in] arg Operator argument to the superoperator.
     */
    qm_operator operator()(const qm_operator& arg) const
    {
        return this->action(arg);
    }

    /**
     * Executes the superoperator on the operator \p arg. Must be
     * implemented by subclass.
     *
     * \param [in] arg Operator argument to the superoperator.
     */
    virtual qm_operator action(const qm_operator& arg) const
    {
        /* return zero by default */
        std::vector<real> zeros(m_num_levels, 0.0);
        return qm_operator(zeros);
    }

    /**
     * Get number of levels.
     */
    unsigned int get_num_levels() const { return m_num_levels; }
};

/**
 * Lindblad form relaxation superoperator.
 * \ingroup MBSOLVE_LIB
 */
class qm_lindblad_relaxation : public qm_superop
{
private:
    std::vector<std::vector<real> > m_scattering;
    std::vector<real> m_dephasing;

public:
    /**
     * Constructs the Lindblad relaxation superoperator using the matrix
     * \p rates.
     *
     * \param [in] rates Matrix with non-negative entries, whose off-diagonal
     * terms with index (m, n) describe the scattering rates from energy level
     * m to energy level n. The main-diagonal entries (m, m) of the matrix
     * represent pure dephasing.
     */
    explicit qm_lindblad_relaxation(
        const std::vector<std::vector<real> >& rates)
      : qm_superop(rates.size()),
        m_scattering(rates.size(), std::vector<real>(rates.size(), 0)),
        m_dephasing((rates.size() * (rates.size() - 1)) / 2, 0)
    {
        /* TODO exception or assert: rates not quadratic? */
        /* TODO exception or assert: all physical constraints fulfilled? */
        /* TODO in particular pure dephasing issue?! */
        /* TODO last main diagonal element must be zero? */

        /** Implementation info: see https://doi.org/10.1364/OE.19.016784
         * Eqs. (9ab)
         */

        /* determine relaxation rates */
        for (int m = 0; m < m_num_levels; m++) {
            real relaxation_rate = 0;

            for (int j = 0; j < m_num_levels; j++) {
                if (j != m) {
                    relaxation_rate += rates[j][m];
                    m_scattering[m][j] = rates[m][j];
                }
            }
            m_scattering[m][m] = -relaxation_rate;
        }

        /* determine dephasing rate for transition mn */
        unsigned int idx_dephasing = 0;
        for (int n = 0; n < m_num_levels; n++) {
            for (int m = 0; m < n; m++) {
                real dephasing = 0;

                /* dephasing due to relaxation */
                dephasing -= 0.5 *
                    (std::abs(m_scattering[m][m]) +
                     std::abs(m_scattering[n][n]));
                /* pure dephasing */
                dephasing -= 0.5 * (rates[m][m] + rates[n][n]);

                m_dephasing[idx_dephasing] = dephasing;
                idx_dephasing++;
            }
        }
    }

    /**
     * Executes the superoperator on the operator \p arg.
     *
     * \param [in] arg Operator argument to the superoperator.
     */
    qm_operator action(const qm_operator& arg) const
    {
        /* TODO exception or assert: operator size does not match */
        /* TODO in particular off-diagonal elements/rates? */

        /** Implementation info: see https://doi.org/10.1364/OE.19.016784
         * Eqs. (9ab)
         */

        std::vector<real> populations(m_num_levels);

        for (int m = 0; m < m_num_levels; m++) {
            real pop = 0;
            for (int n = 0; n < m_num_levels; n++) {
                pop += m_scattering[m][n] * arg.get_main_diagonal()[n];
            }
            populations[m] = pop;
        }

        std::vector<std::complex<real> > coherences = arg.get_off_diagonal();
        for (int i = 0; i < coherences.size(); i++) {
            coherences[i] *= m_dephasing[i];
        }

        return qm_operator(populations, coherences);
    }

    /**
     * Gets the scattering rates.
     */
    const std::vector<std::vector<real> >& get_scattering_rates() const
    {
        return m_scattering;
    }

    /**
     * Gets the dephasing rates.
     */
    const std::vector<real>& get_dephasing_rates() const
    {
        return m_dephasing;
    }
};

/**
 * Provides the quantum mechanical description of an active \ref region.
 * \ingroup MBSOLVE_LIB
 */
class qm_description
{
private:
protected:
    /* number of levels */
    unsigned int m_num_levels;

    /* density of charge carriers */
    real m_carrier_density;

    /*
      real m_period_length;
    */

    /* hamiltonian */
    qm_operator m_hamiltonian;

    /* dipole moment operator */
    qm_operator m_dipole_op;

    /* relaxation superoperator */
    std::shared_ptr<qm_superop> m_relax_superop;

public:
    /**
     * Constructs quantum mechanical description.
     *
     * \param [in] carrier_density    Density of particles in the system.
     * \param [in] hamiltonian        Hamiltonian operator.
     * \param [in] dipole_operator    Dipole moment operator.
     * \param [in] relaxation_superop Relaxation superoperator.
     */
    explicit qm_description(
        real carrier_density,
        const qm_operator& hamiltonian,
        const qm_operator& dipole_operator,
        std::shared_ptr<qm_superop> relaxation_superop)
      : m_num_levels(hamiltonian.get_num_levels()),
        m_carrier_density(carrier_density), m_hamiltonian(hamiltonian),
        m_dipole_op(dipole_operator), m_relax_superop(relaxation_superop)
    {
        /* TODO: assert or exception: level count different? */
    }

    virtual ~qm_description() {}

    /**
     * Gets number of levels.
     */
    unsigned int get_num_levels() const { return m_num_levels; }

    /**
     * Gets carrier density.
     */
    real get_carrier_density() const { return m_carrier_density; }

    /**
     * Gets Hamiltonian.
     */
    const qm_operator& get_hamiltonian() const { return m_hamiltonian; }

    /**
     * Gets dipole moment operator.
     */
    const qm_operator& get_dipole_operator() const { return m_dipole_op; }

    /**
     * Gets relaxation superoperator.
     */
    std::shared_ptr<qm_superop> get_relaxation_superop() const
    {
        return m_relax_superop;
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
    /**
     * Constructs a two-level system description with the Hamiltonian
     * \f[ \hat H = \frac{\hbar}{2} \begin{bmatrix} -\omega_0 & 0 \\ 0 &
     \omega_0 \end{bmatrix} \f]
     * and the dipole moment operator \f[ \hat \mu =
     \begin{bmatrix} 0 & e d \\ e d & 0 \end{bmatrix}, \f]
    * where \f$ e \f$ is the elementary charge.

     * \param [in] carrier_density         Density of particles in the system.
     * \param [in] transition_freq         = \f$ \omega_0 \f$
     * \param [in] dipole_moment           = \f$ d \f$
     * \param [in] scattering_rate         Describes the decay of the
     * population inversion towards \p equilibrium_inversion.
     * \param [in] dephasing_rate          Describes the decay of the
     * coherence term \f$ \rho_{12} = \rho_{21}^* \f$.
     * \param [in] equilibrium_inversion   Equilibrium value of the population
     * inversion. By default, all particles are in the lower energy level.
     */
    explicit qm_desc_2lvl(
        real carrier_density,
        real transition_freq,
        real dipole_moment,
        real scattering_rate,
        real dephasing_rate,
        real equilibrium_inversion = -1.0)
      : qm_description(
            carrier_density,
            qm_operator(
                { -HBAR * transition_freq / 2, +HBAR * transition_freq / 2 }),
            qm_operator({ 0, 0 }, { E0 * dipole_moment }),
            std::make_shared<qm_lindblad_relaxation>(
                std::vector<std::vector<real> >{
                    { 2 * dephasing_rate - scattering_rate,
                      scattering_rate * 0.5 * (1 - equilibrium_inversion) },
                    { scattering_rate * 0.5 * (1 + equilibrium_inversion),
                      0 } })),
        m_trans_freq(transition_freq), m_dipole_mom(dipole_moment),
        m_scattering(scattering_rate), m_dephasing(dephasing_rate),
        m_equi_inv(equilibrium_inversion)
    {
        /* TODO: initial density matrix argument? default? */
        /*
           const qm_operator& rho_init
        */

        /* TODO exception if 2 * dephasing_rate - scattering_rate < 0 */
    }

    ~qm_desc_2lvl() {}

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
}

#endif
