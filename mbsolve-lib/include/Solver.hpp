#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <map>
#include <vector>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Types.hpp>

namespace mbsolve {

class ISolverFactory;

class Solver
{
private:
    bool m_initialized;
    static std::map<std::string, ISolverFactory *> m_solvers;

protected:
    std::string m_name;

    virtual void do_setup(const Device& device, const Scenario& scenario) { }

    virtual void do_run(std::vector<Result *>& results) { }

    virtual void do_cleanup() { }

public:
    explicit Solver(std::string name);

    ~Solver();

    const std::string& name();

    void setup(const Device& device, const Scenario& scenario);

    void run(std::vector<Result *>& results);

    static void register_new(const std::string& name, ISolverFactory *factory);

    static Solver *create(const std::string& name);

};

class ISolverFactory
{
public:
    virtual Solver* create() = 0;
};

template<typename T>
class SolverFactory : ISolverFactory
{
public:
    explicit SolverFactory(const std::string& name) {
	Solver::register_new(name, this);
    }

    Solver* create() { return new T; }
};

}

#endif
