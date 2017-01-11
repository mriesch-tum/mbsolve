#ifndef OPENMP_2LVL_PC_H
#define OPENMP_2LVL_PC_H

#include <Solver.hpp>

namespace mbsolve {

struct sim_constants
{
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    real w12;
    real d12;
    real tau1;
    real gamma12;

    unsigned int idx_start;
    unsigned int idx_end;

    real d_x;
    real d_t;

    real dm11_init;
    real dm22_init;
};

class SolverOMP_2lvl_pc : public ISolver
{
public:
    SolverOMP_2lvl_pc(const Device& device, const Scenario& scenario);

    ~SolverOMP_2lvl_pc();

    std::string getName() const;

    void run() const;

private:

    //    CUDADensityMatrix *m_dm;

    real *m_h;
    real *m_e;
    real *m_e_est;

    //    std::vector<CopyListEntry *> m_copyListBlack;
    //    std::vector<CopyListEntry *> m_copyListRed;
};

}

#endif
