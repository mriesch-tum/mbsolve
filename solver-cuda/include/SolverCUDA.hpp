#include <Solver.hpp>

namespace mbsolve {

class SolverCUDA : public Solver
{
public:
    SolverCUDA();

    ~SolverCUDA();

    void setup();

    void cleanup();

    void run();
};

}
