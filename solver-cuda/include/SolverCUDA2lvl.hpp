#ifndef SOLVERCUDA2LVL_H
#define SOLVERCUDA2LVL_H

#include <Solver.hpp>

namespace mbsolve {

static const unsigned int NumEntries = 4;
static const unsigned int NumMultistep = 5;
static const unsigned int MaxRegions = 8;

struct sim_constants
{
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    real w12;
    real d12;
    real gamma1;
    real gamma2;

    unsigned int idx_start;
    unsigned int idx_end;

    real d_x;
    real d_t;
};

class DensityMatrix
{
public:
    DensityMatrix();

    ~DensityMatrix();

    void initialize(unsigned int numGridPoints);

    void next();

    __device__ __inline__ real *OldDM(unsigned int entry) const;

    __device__ __inline__ real *NewDM(unsigned int entry) const;

    __device__ __inline__ real *RHS(unsigned int entry, unsigned int row) const;

private:
    real *dm_a[NumEntries];
    real *dm_b[NumEntries];
    real *rhs[NumEntries][NumMultistep];

    bool a_is_old;
    unsigned int head;
};


class SolverCUDA2lvl : public Solver
{
public:
    SolverCUDA2lvl();

    ~SolverCUDA2lvl();

    void do_setup(const Device& device, Scenario& scenario);

    void do_cleanup();

    void do_run(std::vector<Result *>& results);

private:
    cudaStream_t comp_maxwell;
    cudaStream_t comp_density;
    cudaStream_t copy;

    DensityMatrix dm;

    real *h;
    real *e;
};

}

#endif
