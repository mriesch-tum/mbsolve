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

#ifndef DEVICE_H
#define DEVICE_H

#include <string>
#include <vector>
#include <Quantity.hpp>
#include <types.hpp>

namespace mbsolve {

class Region
{
public:
    Region() : RelPermeability(1.0),
	       RelPermittivity(1.0),
	       Overlap(1.0) { }

    std::string Name;

    /* dimensions */
    Quantity X0;
    Quantity XDim;

    /* electromagnetic properties */
    Quantity RelPermeability; /* default 1.0 */
    Quantity RelPermittivity; /* default 1.0 */

    Quantity Overlap;
    Quantity Losses;

    /* carrier properties */
    Quantity DopingDensity;
    Quantity PeriodLength;

    /*
      use some sparse structures that provide unique access to a given
      element. transition frequencies are dense, but coupling and anticrossing
      are likely to be sparse. use operator()(int lvl1, int lvl2)
     */

    /* transistion frequencies */
    /* (n-1)*n/2, dense */
    std::vector<Quantity> TransitionFrequencies;

    /* dipole moments */
    /* max. (n-1)*n/2, probably sparse */
    std::vector<Quantity> DipoleMoments;

    /* anticrossing energies */
    /* max. n-2, probably sparse */
    std::vector<Quantity> AnticrossingEnergies;

    /* scattering matrix */
    /* n^2, probably dense */
    std::vector<Quantity> ScatteringRates;

    /* dephasing rates */
    /* (n-1)*n/2, dense */
    std::vector<Quantity> DephasingRates;


};

class Device
{
public:
    std::string Name;
    std::vector<Region> Regions;

    /* TODO: what happens if no Regions inserted */

    Quantity XDim() const {
	Quantity total;
	std::vector<Region>::const_iterator it;
	for (it = Regions.begin(); it != Regions.end(); it++) {
	    total += it->XDim;
	}
	return total;
    }

    Quantity MinRelPermittivity() const {
	Quantity min(1e12);
	std::vector<Region>::const_iterator it;
	for (it = Regions.begin(); it != Regions.end(); it++) {
	    if (it->RelPermittivity < min) {
		min = it->RelPermittivity;
	    }
	}
	return min;
    }

};

}

#endif
