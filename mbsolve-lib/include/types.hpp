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

#ifndef MBSOLVE_TYPES_H
#define MBSOLVE_TYPES_H

#include <complex>
#include <string>

namespace mbsolve {

/* TODO: define type switch single/double */

/**
 * May be set to single or double precision.
 * \ingroup MBSOLVE_LIB
 */
typedef double real;
//typedef float real;

/**
 * May be set to single or double precision.
 * \ingroup MBSOLVE_LIB
 */
typedef std::complex<real> complex;

/**
 * Reduced Planck's constant.
 * \ingroup MBSOLVE_LIB
 */
static const real HBAR = 1.05457266e-34;

/**
 * Vacuum permeability.
 * \ingroup MBSOLVE_LIB
 */
static const real MU0 = M_PI * 4e-7;

/**
 * Vacuum permittivity.
 * \ingroup MBSOLVE_LIB
 */
static const real EPS0 = 8.854187817e-12;

/**
 * Elementary charge.
 * \ingroup MBSOLVE_LIB
 */
static const real E0 = 1.60217733e-19;

}

#endif
