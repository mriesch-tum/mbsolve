#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <Device.hpp>
#include <Scenario.hpp>

namespace mbsolve {

class Solver
{
private:
    bool m_initialized;

protected:
    std::string m_name;

public:
    explicit Solver(std::string name, bool setup = false) :
	m_name(name), m_initialized(false) {
	if (setup) {
	    this->setup();
	    this->m_initialized = true;
	}
    }

    ~Solver() {
	if (m_initialized) {
	    this->cleanup();
	}
    }

    const std::string& name() { return m_name; }

    virtual void setup() { }

    virtual void cleanup() { }

    virtual void run() { }
};

}

#endif
