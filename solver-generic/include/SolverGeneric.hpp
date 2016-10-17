#ifndef SOLVERGENERIC_H
#define SOLVERGENERIC_H

#include <Solver.hpp>

namespace mbsolve {

class SolverGeneric : public Solver
{
public:
    SolverGeneric();

    ~SolverGeneric();

    void do_setup(const Device& device, const Scenario& scenario);

    void do_cleanup();

    void do_run(std::vector<Result *>& results);
};

}

#endif
