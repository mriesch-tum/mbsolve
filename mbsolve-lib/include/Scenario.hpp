#ifndef SCENARIO_H
#define SCENARIO_H

#include <string>
#include <vector>
#include <Record.hpp>
#include <Source.hpp>

namespace mbsolve {

class Scenario
{
public:
    std::string Name;

    unsigned int NumTimeSteps;

    unsigned int NumGridPoints;

    real TimeStepSize;

    real GridPointSize;

    real SimEndTime;

    std::vector<IRecord *> Records;

    /* TODO: add sources vector */
    std::vector<ISource *> Sources;


};

}

#endif
