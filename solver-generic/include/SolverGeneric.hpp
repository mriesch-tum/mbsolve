#ifndef SOLVERGENERIC_H
#define SOLVERGENERIC_H

#include <Solver.hpp>

namespace mbsolve {

class SolverGeneric : public ISolver
{
public:
    SolverGeneric(const Device& device, const Scenario& scenario);

    ~SolverGeneric();

    std::string getName() const;

    void run() const;
};

}

#endif
