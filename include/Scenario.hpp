#ifndef SCENARIO_H
#define SCENARIO_H

#include <string>

namespace mbsolve {

class Scenario
{
public:
    std::string Name;

    unsigned int N_t;

    unsigned int N_x;

    real d_t;

    real d_x;

    real t_e;

    real mod_a;

    real mod_f;

};

}

#endif
