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

#ifndef MBSOLVE_ALGO_LINDBLAD_NOOP_H
#define MBSOLVE_ALGO_LINDBLAD_NOOP_H

#include <mbsolve/lib/qm_description.hpp>

namespace mbsolve {

/**
 * Does not solve the Lindblad equation. Using this algorithm, mbsolve is
 * effictively a 1D Maxwell's equations solver.
 * Use only internally to implement solvers.
 * \ingroup MBSOLVE_LIB_INT
 */
template<unsigned int num_lvl = 0>
class lindblad_noop
{
private:
    class dummy
    {
    public:
        void* operator new[](std::size_t count) { return new char; }

        void operator delete[](void* ptr) { delete (ptr); }
    };

public:
    static std::string name() { return "noop"; }

    typedef char density;

    typedef char sim_constants;

    typedef std::allocator<sim_constants> allocator;

    static inline sim_constants get_qm_constants(
        std::shared_ptr<const qm_description> qm,
        real time_step)
    {
        return 0;
    }

    static inline void
    update(const sim_constants& sc, density& d, real e, real* p_t)
    {
        /* no-op */
    }

    static inline real calc_inversion(const density& d) { return 0; }

    static inline real calc_population(const density& d, unsigned int idx)
    {
        return 0;
    }

    static inline density get_density()
    {
        /* return empty density */
        return 0; // density();
    }

    static inline density get_density(const qm_operator& /* op -- unused*/)
    {
        return get_density();
    }
};
}

#endif
