#ifndef SCENARIO_H
#define SCENARIO_H

#include <string>

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

};

}

#endif
