#ifndef DEVICE_H
#define DEVICE_H

#include <string>
#include <vector>
#include <Quantity.hpp>
#include <Types.hpp>

namespace mbsolve {

class Region
{
public:
    std::string Name;

    /* dimensions */
    Quantity X0;
    Quantity XDim;

    /* electromagnetic properties */
    Quantity Permeability;
    Quantity Permittivity;
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



};

}

#endif
