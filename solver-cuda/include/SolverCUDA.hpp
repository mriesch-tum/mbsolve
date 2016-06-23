#include <Solver.hpp>

namespace mbsolve {

class SolverCUDA : public Solver
{
public:
    SolverCUDA();

    ~SolverCUDA();

    void do_setup(const Device& device, const Scenario& scenario);

    void do_cleanup();

    void do_run();
};

}
