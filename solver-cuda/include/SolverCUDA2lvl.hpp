#ifndef SOLVERCUDA2LVL_H
#define SOLVERCUDA2LVL_H

#include <Solver.hpp>

namespace mbsolve {

class SolverCUDA2lvl : public Solver
{
public:
    SolverCUDA2lvl();

    ~SolverCUDA2lvl();

    void do_setup(const Device& device, const Scenario& scenario);

    void do_cleanup();

    void do_run(std::vector<Result *>& results);
};

}

#endif
