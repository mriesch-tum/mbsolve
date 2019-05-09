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

#ifndef MBSOLVE_ALGO_LINDBLAD_CVR_RODR_H
#define MBSOLVE_ALGO_LINDBLAD_CVR_RODR_H

#include <numeric>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <unsupported/Eigen/MatrixFunctions>
#include <internal/coherence_vector_representation.hpp>

namespace mbsolve {

/**
 * Solves the Lindblad equation in coherence vector representation using the
 * generalized Rodrigues formula.
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB_INT
 */
template<unsigned int num_lvl>
class lindblad_cvr_rodr
{
private:
    static const unsigned int vec_len = num_lvl * num_lvl - 1;

    typedef Eigen::Matrix<real, vec_len, vec_len> real_matrix_t;

public:
    static std::string name() { return "cvr-rodr"; }

    /* TODO outsource to cvr class */
    typedef Eigen::Matrix<real, vec_len, 1> density;

    class sim_constants {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /* does material even have a quantum mechanical description? */
        bool has_qm;

        /* polarization constant */
        real M_CP;

        /* time step size */
        real d_t;

        /* time-independent part of the Lindblad equation rhs */
        real_matrix_t M;

        /* dipole moment Liouvillian in cvr form (matrix) + squared */
        real_matrix_t U;

        /* rodrigues formula precalc (two-level case) */
        real theta_1;
        real_matrix_t U2;

        /* rodrigues formula precalc (generalized case) */
        real_matrix_t coeff_1[vec_len/2];
        real_matrix_t coeff_2[vec_len/2];
        real theta[vec_len/2];

        /* time-independent propagator ( contains time-independent part of
         * Hamiltonian plus dissipation superoperator contribution)
         * applied on a half step.
         */
        real_matrix_t A1;

        /* dipole operator in cvr form (vector) */
        density v;

        /* equilibrium vector terms */
        density d_eq;
        density d_in;
    };

    typedef Eigen::aligned_allocator<sim_constants> allocator;

private:
    static inline void
    fill_coefficients(sim_constants& sc)
    {
        /* get eigenvalues and eigenvectors in real form */
        Eigen::RealSchur<real_matrix_t> schur(sc.U);
        real_matrix_t eigenval = schur.matrixT();
        real_matrix_t eigenvec = schur.matrixU();

        /* creating sorting order (descending eigenvalues) */
        std::vector<size_t> perm_idx(vec_len - 1);
        std::iota(perm_idx.begin(), perm_idx.end(), 1);
        std::sort(perm_idx.begin(), perm_idx.end(),
                  [&eigenval](size_t i1, size_t i2) {
                      return std::abs(eigenval(i1, i1 - 1)) >
                          std::abs(eigenval(i2, i2 - 1));
                  });

        /* TODO optimize
         * ignore eigenvalues = 0
         * group eigenvalues with multiplicity >= 2
         */
        for (int i = 0; i < vec_len/2; i++) {
            unsigned int i1 = perm_idx[i];

            real_matrix_t b = real_matrix_t::Zero();
            b(i1, i1 - 1) = +1.0;
            b(i1 - 1, i1) = -1.0;

            sc.coeff_1[i] = eigenvec * b * eigenvec.transpose();
            sc.coeff_2[i] = eigenvec * b * b * eigenvec.transpose();
            sc.theta[i] = eigenval(i1, i1 - 1);
            std::cout << "theta: "<< std::endl << sc.theta[i] << std::endl;
        }
        std::cout << "eigenval: " << eigenval << std::endl;
    }

public:
    /*
     * Determines the simulation constants of the quantum mechanical
     * description for the generalized Rodrigues formula.
     */
    static sim_constants
    get_qm_constants(std::shared_ptr<const qm_description> qm,
                     real time_step) {
        sim_constants sc;

        if (qm) {
            /* quantum mechanical description present */
            sc.has_qm = true;

            /* check whether number of levels matches solver */
            if (qm->get_num_levels() != num_lvl) {
                throw std::invalid_argument("Number of energy levels does not "
                                            "match selected solver!");
            }

            /* create coherence vector representation */
            cv_representation cvr(qm);

            /* polarization constant */
            sc.M_CP = 0.5 * qm->get_carrier_density();

            /* time step size */
            sc.d_t = time_step;

            /* dissipation term in cvr */
            real_matrix_t G = cvr.get_relaxation_superop();

            /* time-independent Hamiltonian contribution in cvr */
            real_matrix_t M_0 = cvr.get_hamiltonian();

            /* time-independent part of the Lindblad equation rhs */
            sc.M = M_0 + G;

            /* dipole moment Liouvillian in cvr form (matrix) + squared */
            sc.U = -cvr.get_dipole_operator();
            sc.U2 = sc.U * sc.U;
            sc.theta_1 = sqrt(pow(sc.U(0, 1), 2) + pow(sc.U(0, 2), 2)
                              + pow(sc.U(1, 2), 2));

            /* time-independent propagator, half step */
            sc.A1 = (sc.M * time_step/2).exp();

            /* prepare time-dependent propagator */
            fill_coefficients(sc);

            /* dipole moment operator in cvr form (vector) */
            sc.v = cvr.get_dipole_operator_vec();

            /* equilibrium terms */
            sc.d_eq = cvr.get_equilibrium_vec();
            real eps = std::numeric_limits<real>::epsilon();
            if (sc.d_eq.isZero(eps)) {
                sc.d_in = density::Zero();
            } else {
                /* solve equation system M * d_in = d_eq */
                sc.d_in = sc.M.fullPivLu().solve(sc.d_eq);
                real err = (sc.M * sc.d_in - sc.d_eq).norm()/sc.d_eq.norm();
                if (err > 1e-3) {
                    throw std::invalid_argument("Time-indepent matrix not "
                                                "invertible!");
                }
            }
        } else {
            /* no quantum mechanical description */
            sc.has_qm = false;
            sc.M_CP = 0.0;
            sc.d_t = 0.0;
            sc.M = real_matrix_t::Zero();
            sc.U = real_matrix_t::Zero();
            sc.U2 = real_matrix_t::Zero();
            sc.theta_1 = 0.0;
            sc.A1 = real_matrix_t::Zero();
            sc.v = density::Zero();
            sc.d_eq = density::Zero();
            sc.d_in = density::Zero();

            /* TODO initialize properly?:
               real_matrix_t coeff_1[num_lvl/2];
               real_matrix_t coeff_2[num_lvl/2];
               real theta[num_lvl/2];
            */
        }

        return sc;
    }

    /**
     * Updates the density matrix and the derivative of the polarization.
     */
    static inline void
    update(const sim_constants &sc, density& d, real e, real *p_t) {
        if (sc.has_qm) {
            /* time-indepedent half step */
            density d1 = sc.A1 * (d + sc.d_in) - sc.d_in;

            /* determine time-dependent propagator */
            real_matrix_t A2;

            if (num_lvl == 2) {
                /* original Rodrigues' formula */
                A2 = sin(sc.theta_1 * e * sc.d_t)/sc.theta_1 * sc.U
                    + (1 - cos(sc.theta_1 * e * sc.d_t))/
                    (sc.theta_1 * sc.theta_1) * sc.U2
                    + real_matrix_t::Identity();
            } else {
                A2 = real_matrix_t::Identity();
                for (int i = 0; i < vec_len/2; i++) {
                    /* TODO nolias()? */
                    A2 += sin(sc.theta[i] * e * sc.d_t) * sc.coeff_1[i]
                        + (1 - cos(sc.theta[i] * e * sc.d_t)) *
                        sc.coeff_2[i];
                }
            }

            /* time-dependent full step */
            density d2 = A2 * d1;

            /* time-indepedent half step */
            d = sc.A1 * (d2 + sc.d_in) - sc.d_in;

            /* update polarization */
            *p_t = sc.M_CP * sc.v.transpose() * (sc.M * d + sc.d_eq);
        }
    }

    /*
     * Returns the population inversion (which is the difference
     * \f$ \rho_{22} - \rho_{11} \f$).
     */
    static inline real
    calc_inversion(const density& d) {
        return d(num_lvl * (num_lvl - 1));
    }

    /*
     * Returns the population specified by \param idx (zero-based).
     */
    static inline real
    calc_population(const density& d, unsigned int idx) {
        return cv_representation::calc_population<num_lvl, vec_len>(d, idx);
    }

    /*
     * Converts the density matrix \param op into the coherence vector
     * representation.
     */
    static inline density
    get_density(const qm_operator& op) {
        return cv_representation::coherence_vector<num_lvl, vec_len>(op);
    }

    /*
     * Returns a zero coherence vector (for initialization purposes only).
     */
    static inline density
    get_density() {
        return density::Zero();
    }
};

}

#endif
