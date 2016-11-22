#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Types.hpp>

namespace mbsolve {

class ISolver
{
protected:
    Scenario m_scenario;
    Device m_device;

    std::vector<Result *> m_results;
public:
    ISolver(const Device& device, const Scenario& scenario) :
	m_device(device), m_scenario(scenario) { }

    virtual ~ISolver();

    const Scenario& getScenario() const { return m_scenario; }

    const Device& getDevice() const { return m_device; }

    virtual std::string getName() const = 0;

    virtual void run() const = 0;

    const std::vector<Result *>& getResults() const
    {
	return m_results;
    }
};

class ISolverFactory;

class Solver
{
private:
    static std::map<std::string, ISolverFactory *> m_factories;
    ISolver *m_solver;

public:
    Solver(const std::string& name, const Device& device,
	   const Scenario& scenario);

    ~Solver();

    std::string getName() const;

    void run() const;

    const std::vector<Result *>& getResults() const;

    static void registerFactory(const std::string& name,
				ISolverFactory *factory);
};

class ISolverFactory
{
public:
    ISolverFactory() { }
    virtual ~ISolverFactory() { }
    virtual ISolver* createInstance(const Device& device,
				    const Scenario& scenario) const = 0;
};

template<typename T>
class SolverFactory : ISolverFactory
{
private:
    std::string m_name;
public:
    explicit SolverFactory(const std::string& name) : m_name(name) {
	Solver::registerFactory(name, this);
    }

    ISolver* createInstance(const Device& device,
			    const Scenario& scenario) const {
	return new T(device, scenario);
    }

    const std::string& getName() const { return m_name; }
};

}

#endif
