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

#ifndef MBSOLVE_LIB_H
#define MBSOLVE_LIB_H

/**
 * \mainpage
 *
 * The mbsolve project provides an open-source solver tool for the
 * Maxwell-Bloch equations.
 */

/**
 * \defgroup MBSOLVE_LIB mbsolve-lib
 * Provides data types and base classes.
 */

#include <mbsolve/lib/device.hpp>
#include <mbsolve/lib/material.hpp>
#include <mbsolve/lib/qm_description.hpp>
#include <mbsolve/lib/record.hpp>
#include <mbsolve/lib/result.hpp>
#include <mbsolve/lib/scenario.hpp>
#include <mbsolve/lib/solver.hpp>
#include <mbsolve/lib/source.hpp>
#include <mbsolve/lib/types.hpp>
#include <mbsolve/lib/writer.hpp>

/**
 * \defgroup MBSOLVE_LIB_INT mbsolve-lib-internal
 * Provides helper functions and data structures in order to factor out
 * common tasks of solvers (and writers).
 */

#endif
