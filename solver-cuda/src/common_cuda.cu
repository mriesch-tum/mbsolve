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

#include <common_cuda.hpp>

namespace mbsolve {

/* material properties in constant GPU memory */
/*
__device__ __constant__
sim_constants_2lvl l_sim_consts[MB_CUDA_MAX_MATERIALS];
*/

/* source properties in constant GPU memory */
__device__ __constant__ sim_source l_sim_sources[MB_CUDA_MAX_SOURCES];

/* copy list in constant GPU memory */
__device__ __constant__ copy_list_entry_dev l_copy_list[MB_CUDA_MAX_CLE];

}
