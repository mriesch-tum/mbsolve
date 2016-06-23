#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <vector>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Types.hpp>

namespace mbsolve {

class Solver
{
private:
    bool m_initialized;

protected:
    std::string m_name;

    virtual void do_setup(const Device& device, const Scenario& scenario) { }

    virtual void do_run(std::vector<Result *>& results) { }

    virtual void do_cleanup() { }

public:
    explicit Solver(std::string name) : m_initialized(false), m_name(name) { }

    ~Solver() {
	if (m_initialized) {
	    do_cleanup();
	}
    }

    const std::string& name() { return m_name; }

    void setup(const Device& device, const Scenario& scenario) {
	do_setup(device, scenario);
	m_initialized = true;
    }

    void run(std::vector<Result *>& results) {
	return do_run(results);
    }
};

}

#endif
