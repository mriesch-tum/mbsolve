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

#ifndef MBSOLVE_COHERENCE_VECTOR_REPRESENTATION_H
#define MBSOLVE_COHERENCE_VECTOR_REPRESENTATION_H

#include <Eigen/Dense>
#include <qm_description.hpp>

namespace mbsolve {

/**
 * Provides helper functions for the coherence vector representation.
 * \ingroup MBSOLVE_LIB_INT
 */
class cv_representation
{
private:
    const unsigned int m_num_levels;
    const unsigned int m_num_adj;

    typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> real_matrix_t;
    typedef Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic>
    complex_matrix_t;

    real_matrix_t m_hamiltonian;
    real_matrix_t m_dipole_op;
    real_matrix_t m_dipole_op_vec;
    real_matrix_t m_relax_superop;
    real_matrix_t m_equi;
    real_matrix_t m_d_init;

    std::vector<complex_matrix_t> m_generators;

    void
    setup_generators()
    {
        m_generators.clear();

        /* set up SU(N) generators u -- real part off-diagonal elements */
        for (int k = 0; k < m_num_levels; k++) {
            for (int j = 0; j < k; j++) {
                complex_matrix_t g =
                    complex_matrix_t::Zero(m_num_levels, m_num_levels);
                g(j, k) = 1;
                g(k, j) = 1;
                m_generators.push_back(g);
            }
        }

        /* set up SU(N) generators v -- imag part off-diagonal elements */
        for (int k = 0; k < m_num_levels; k++) {
            for (int j = 0; j < k; j++) {
                complex_matrix_t g =
                    complex_matrix_t::Zero(m_num_levels, m_num_levels);
                g(j, k) = complex(0, -1);
                g(k, j) = complex(0, +1);
                m_generators.push_back(g);
            }
        }

        /* set up SU(N) generators w -- main-diagonal elements */
        for (int l = 0; l < m_num_levels - 1; l++) {
            complex_matrix_t g =
                complex_matrix_t::Zero(m_num_levels, m_num_levels);
            int j = 0;

            for (j = 0; j <= l; j++) {
                g(j, j) = -1.0;
            }
            g(j, j) = (l + 1);

            int k = ((l + 1) * (l + 2))/2;
            g /= sqrt(k);
            m_generators.push_back(g);
        }
    }

    /* TODO: move to qm_operator? */
    /**
     * Transforms qm_operator to complex matrix.
     */
    complex_matrix_t
    convert_qm_operator(const qm_operator& op)
    {
        complex_matrix_t op_mat =
            complex_matrix_t::Zero(m_num_levels, m_num_levels);

        for (int i = 0; i < m_num_levels; i++) {
            /* main diagonal elements */
            op_mat(i, i) = op.get_main_diagonal()[i];
        }

        for (int i = 0, row = 0, col = 1; i < op.get_off_diagonal().size();
             i++) {

            /* off-diagonal, top half */
            op_mat(row, col) = op.get_off_diagonal()[i];

            /* off-diagonal, bottom half */
            op_mat(col, row) = std::conj(op.get_off_diagonal()[i]);

            if (row == col - 1) {
                row = 0;
                col++;
            } else {
                row++;
            }
        }

        return op_mat;
    }

    qm_operator
    convert_qm_operator(complex_matrix_t mat)
    {
        /* main diagonal elements */
        std::vector<real> main_diag;
        for (int i = 0; i < m_num_levels; i++) {
            main_diag.push_back(mat(i, i).real());
        }

        /* TODO: hermiticity check? */
        /* off-diagonal elements */
        std::vector<complex> off_diag;
        for (int i = 1; i < m_num_levels; i++) {
            for (int j = 0; j < i; j++) {
                off_diag.push_back(mat(j, i));
            }
        }

        return qm_operator(main_diag, off_diag);
    }

    real_matrix_t
    calc_liouvillian(const qm_operator& op)
    {
        complex_matrix_t op_mat = convert_qm_operator(op);

        /* determine Liouvillian with respect to operator op */
        real_matrix_t ret(m_num_adj, m_num_adj);
        for (int i = 0; i < m_num_adj; i++) {
            for (int j = 0; j < m_num_adj; j++) {
                complex_matrix_t result = op_mat *
                    (m_generators[i] * m_generators[j] -
                     m_generators[j] * m_generators[i]);
                complex c = complex(0, 1) * 0.5 * result.trace();
                ret(i, j) = c.real() * 1.0/HBAR;
            }
        }
        return ret;
    }

    real_matrix_t
    calc_op_vec(const qm_operator& op)
    {
        complex_matrix_t op_mat = convert_qm_operator(op);
        real_matrix_t ret(m_num_adj, 1);

        for (int i = 0; i < m_num_adj; i++) {
            complex_matrix_t m = op_mat * m_generators[i];
            ret(i) = m.real().trace();
        }

        return ret;
    }

    real_matrix_t
    calc_relaxation_superop(std::shared_ptr<qm_superop> op)
    {
        real_matrix_t ret = real_matrix_t::Zero(m_num_adj, m_num_adj);

        for (int i = 0; i < m_num_adj; i++) {
            for (int j = 0; j < m_num_adj; j++) {
                qm_operator a = convert_qm_operator(m_generators[j]);
                qm_operator b = (*op)(a);
                complex_matrix_t c = convert_qm_operator(b);
                complex_matrix_t d = c * m_generators[i];
                complex e = 0.5 * d.trace();
                ret(i, j) = e.real();
            }
        }
        return ret;
    }

    real_matrix_t
    calc_equi_vec(std::shared_ptr<qm_superop> op)
    {
        real_matrix_t ret = real_matrix_t::Zero(m_num_adj, 1);
        complex_matrix_t a =
            complex_matrix_t::Identity(m_num_levels, m_num_levels);
        qm_operator b = convert_qm_operator(a);
        qm_operator c = (*op)(b);
        complex_matrix_t d = convert_qm_operator(c);

        for (int i = 0; i < m_num_adj; i++) {
            complex_matrix_t m = d * m_generators[i];
            ret(i) = m.real().trace() * 1.0/m_num_levels;
        }
        return ret;
    }

public:
    explicit cv_representation(std::shared_ptr<qm_description> qm_desc) :
        m_num_levels(qm_desc->get_num_levels()),
        m_num_adj(qm_desc->get_num_levels() * qm_desc->get_num_levels() - 1),
        m_hamiltonian(m_num_adj, m_num_adj),
        m_dipole_op(m_num_adj, m_num_adj),
        m_dipole_op_vec(m_num_adj, 1),
        m_relax_superop(m_num_adj, m_num_adj),
        m_equi(m_num_adj, 1),
        m_d_init(m_num_adj, 1)
    {
        /* set up generators of Lie algebra su(num_levels) */
        setup_generators();

        /* calculate Liouvillian of Hamiltonian */
        m_hamiltonian = calc_liouvillian(qm_desc->get_hamiltonian());

        /* calculate Liouvillian of dipole moment operator */
        m_dipole_op = calc_liouvillian(qm_desc->get_dipole_operator());

        /* calculate dipole moment operator in vector form */
        m_dipole_op_vec = calc_op_vec(qm_desc->get_dipole_operator());

        /* calculate superoperator in coherence vector representation */
        m_relax_superop =
            calc_relaxation_superop(qm_desc->get_relaxation_superop());

        /* calculate equilibrium vector */
        m_equi = calc_equi_vec(qm_desc->get_relaxation_superop());

        /* calculate inital vector */
        m_d_init = calc_op_vec(qm_desc->get_rho_init());
    }

    /**
     * Gets the contribution of the Hamiltonian to the Liouvillian
     * \f$ \mathcal{L} = -\mathrm{i}\hbar^{-1} \left[ H, \cdot \right] \f$.
     */
    const real_matrix_t&
    get_hamiltonian() const
    {
        return m_hamiltonian;
    }

    /**
     * Gets the contribution of the dipole moment operator to the
     * Liouvillian
     * \f$ \mathcal{L} = -\mathrm{i}\hbar^{-1} \left[ \mu, \cdot \right] \f$.
     */
    const real_matrix_t&
    get_dipole_operator() const
    {
        return m_dipole_op;
    }

    /**
     * Gets the dipole moment operator in vector form.
     */
    const real_matrix_t&
    get_dipole_operator_vec() const
    {
        return m_dipole_op_vec;
    }

    /**
     * Gets relaxation superoperator in coherence vector representation.
     */
    const real_matrix_t&
    get_relaxation_superop() const
    {
        return m_relax_superop;
    }

    /**
     * Gets equilibrium vector.
     */
    const real_matrix_t&
    get_equilibrium_vec() const
    {
        return m_equi;
    }

    /**
     * Gets initial coherence vector.
     */
    const real_matrix_t&
    get_initial_vec() const
    {
        return m_d_init;
    }

    /**
     * Calculates population \f$ \rho_{mm} \f$ for a given coherence vector
     * \param d and index \param m.
     */
    template<unsigned int num_lvl, unsigned int num_adj>
    static real
    calc_population(const Eigen::Matrix<real, num_adj, 1>& d, unsigned int m)
    {
        /* TODO d.size() instead of template param? */

        /* TODO test and make general */

        if (num_lvl == 2) {
            if (m == 0) {
                return 0.5 - 0.5 * d(2);
            } else {
                return 0.5 + 0.5 * d(2);
            }
        } else if (num_lvl == 3) {
            switch (m) {
            case 0:
                return 1.0/3 + 0.5 * (-1.0 * d(6) - 1.0/sqrt(3) * d(7));
            case 1:
                return 1.0/3 + 0.5 * (+1.0 * d(6) - 1.0/sqrt(3) * d(7));
            case 2:
                return 1.0/3 + 0.5 * (+2.0/sqrt(3) * d(7));
            default:
                return 0.0;
            }
        } else {
            real sum = 0.0;

            for (int i = 0; i < num_lvl - 1; i++) {
                real factor;

                if (m < i + 1) {
                    factor = -1.0;
                } else if (m == i + 1) {
                    factor = m;
                } else {
                    factor = 0.0;
                }

                int k = ((i + 1) * (i + 2))/2;

                sum += factor/sqrt(k) * d(num_lvl * (num_lvl - 1) + i);
            }

            return 1.0/num_lvl + 0.5 * sum;
        }
    }



};

}

#endif
