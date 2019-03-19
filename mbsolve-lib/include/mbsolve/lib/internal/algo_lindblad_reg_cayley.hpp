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

#ifndef MBSOLVE_ALGO_LINDBLAD_REG_CAYLEY_H
#define MBSOLVE_ALGO_LINDBLAD_REG_CAYLEY_H

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <unsupported/Eigen/MatrixFunctions>
#include <mbsolve/lib/qm_description.hpp>

namespace mbsolve {

/**
 * Solves the Lindblad equation using the symmetric Strang operator splitting,
 * the Cayley transform to approximate the matrix exponential.
 *
 * For details see literature:
 * Bidégaray, B. et al., Introducing physical relaxation terms in Bloch
 * equations, J. Comp. Phys., Vol. 170, Issue 2, 2001
 * https://doi.org/10.1006/jcph.2001.6752.
 *
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB_INT
 */
template<unsigned int num_lvl = 0>
class lindblad_reg_cayley
{
public:
    static std::string name() { return "reg-cayley"; }

    typedef Eigen::Matrix<real, num_lvl, num_lvl> real_matrix_t;
    typedef Eigen::Matrix<std::complex<real>, num_lvl, num_lvl>
        complex_matrix_t;

    /* TODO: is there a Eigen Hermitian Matrix type? */
    typedef complex_matrix_t hermitian_matrix_t;

    typedef hermitian_matrix_t density;

    class sim_constants
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /* time step size */
        real dt;

        /* time-independent Hamiltonian operator */
        hermitian_matrix_t H0;

        /* dipole moment operator */
        hermitian_matrix_t mu;

        /* dissipation rhs -- population part */
        real_matrix_t d1;

        /* dissipation rhs -- dephasing part */
        real_matrix_t d2;

        /* dissipation propagator -- population part */
        real_matrix_t D1;

        /* dissipation propagator -- dephasing part */
        real_matrix_t D2;

        /* carrier density */
        real carrier_dens;
    };

    typedef Eigen::aligned_allocator<sim_constants> allocator;

    static inline hermitian_matrix_t convert_qm_operator(
        const qm_operator& op)
    {
        hermitian_matrix_t op_mat = hermitian_matrix_t::Zero();

        /* main diagonal elements */
        for (int i = 0; i < num_lvl; i++) {
            op_mat(i, i) = op.get_main_diagonal()[i];
        }

        /* off-diagonal elements */
        for (int i = 0, row = 0, col = 1; (row < num_lvl) &&
             (col < num_lvl) && (i < op.get_off_diagonal().size());
             i++) {

            /* top half */
            op_mat(row, col) = op.get_off_diagonal()[i];

            /* bottom half */
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

    static inline sim_constants get_qm_constants(
        std::shared_ptr<const qm_description> qm,
        real time_step)
    {
        sim_constants sc;

        /* time step size */
        sc.dt = time_step;

        if (qm) {
            /* check whether number of levels matches solver */
            if (qm->get_num_levels() != num_lvl) {
                throw std::invalid_argument(
                    "Number of energy levels does not "
                    "match selected solver!");
            }

            /* Hamiltonian */
            sc.H0 = convert_qm_operator(qm->get_hamiltonian());

            /* dipole moment operator */
            sc.mu = convert_qm_operator(qm->get_dipole_operator());

            /* TODO & Warning
             * Here, we assume the Lindblad form.
             * - Is this acceptable for other cases as well?
             */
            std::shared_ptr<qm_lindblad_relaxation> sop =
                std::dynamic_pointer_cast<qm_lindblad_relaxation>(
                    qm->get_relaxation_superop());
            if (!sop) {
                throw std::invalid_argument("Unsupported superoperator!");
            }

            /* dissipation propagator for populations is the matrix
             * exponential of the population scattering rate matrix
             */
            sc.d1 = real_matrix_t::Zero();
            std::vector<std::vector<real> > smat =
                sop->get_scattering_rates();

            for (int row = 0; row < num_lvl; row++) {
                for (int col = 0; col < num_lvl; col++) {
                    sc.d1(row, col) = smat[row][col];
                }
            }
            sc.D1 = (sc.d1 * time_step / 2).exp();

            /* dissipation propagator for dephasing terms is the elementwise
             * exponential of the dephasing rates
             */
            sc.d2 = real_matrix_t::Zero();
            std::vector<real> deph = sop->get_dephasing_rates();
            for (int i = 0, row = 0, col = 1;
                 (row < num_lvl) && (col < num_lvl) && (i < deph.size());
                 i++) {
                sc.d2(row, col) = deph[i];
                sc.d2(col, row) = deph[i];
                if (row == col - 1) {
                    row = 0;
                    col++;
                } else {
                    row++;
                }
            }
            sc.D2 = (sc.d2.array() * time_step / 2).exp();

            /* carrier density */
            sc.carrier_dens = qm->get_carrier_density();
        } else {
            sc.H0 = hermitian_matrix_t::Zero();
            sc.mu = hermitian_matrix_t::Zero();
            sc.d1 = real_matrix_t::Zero();
            sc.D1 = real_matrix_t::Zero();
            sc.D2 = real_matrix_t::Zero();
            sc.d2 = real_matrix_t::Zero();
            sc.carrier_dens = 0.0;
        }

        return sc;
    }

private:
    static inline density propagate_dissipation(
        const sim_constants& sc,
        const density& d)
    {
        /* update coherence terms */
        density d1 = sc.D2.array() * d.array();

        /* update populations */
        d1.diagonal() = sc.D1 * d.diagonal();

        return d1;
    }

    static inline complex_matrix_t matrix_exponential(
        const complex_matrix_t& m)
    {
        /* Cayley approximation from Bidégaray et al. */
        return (complex_matrix_t::Identity() - 0.5 * m).inverse() *
            (complex_matrix_t::Identity() + 0.5 * m);
    }

public:
    static inline void
    update(const sim_constants& sc, density& d, real e, real* p_t)
    {
        /* determine complete Hamiltonian */
        hermitian_matrix_t H = sc.H0 - e * sc.mu;

        /* prepare Hamiltonian propagator */
        complex_matrix_t U2 = matrix_exponential(
            -std::complex<real>(0, 1) * 1.0 / HBAR * sc.dt * H);

        /* dissipation propagator half step */
        density d1 = propagate_dissipation(sc, d);

        /* Hamiltonian propagator full step */
        density d2 = U2 * d1 * U2.adjoint();

        /* dissipation propagator half step */
        density d3 = propagate_dissipation(sc, d2);

        /* return result */
        d = d3;

        /* update polarization term */
        complex_matrix_t rhs_d = sc.d2.array() * d.array();
        rhs_d.diagonal() = sc.d1 * d.diagonal();
        complex_matrix_t rhs =
            -std::complex<real>(0, 1) * 1.0 / HBAR * (sc.H0 * d - d * sc.H0) +
            rhs_d;
        *p_t = sc.carrier_dens * (rhs * sc.mu).trace().real();
    }

    static inline real calc_inversion(const density& d)
    {
        return (d(1, 1) - d(0, 0)).real();
    }

    static inline real calc_population(const density& d, unsigned int idx)
    {
        return d(idx, idx).real();
    }

    static inline density get_density(const qm_operator& op)
    {
        return convert_qm_operator(op);
    }

    static inline density get_density() { return density::Zero(); }
};

}

#endif
