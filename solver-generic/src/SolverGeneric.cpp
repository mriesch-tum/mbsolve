#include <SolverGeneric.hpp>

namespace mbsolve {

static SolverFactory<SolverGeneric> factory("generic");

SolverGeneric::SolverGeneric(const Device& device, const Scenario& scenario) :
    ISolver(device, scenario)
{
}

SolverGeneric::~SolverGeneric()
{
}

std::string
SolverGeneric::getName() const
{
    return std::string("Generic dummy solver");
}

void
SolverGeneric::run(const std::vector<Result *>& results) const
{
}


}
