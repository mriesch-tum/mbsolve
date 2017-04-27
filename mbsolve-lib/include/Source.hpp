#ifndef SOURCE_H
#define SOURCE_H

#include <map>
#include <string>
#include <vector>
#include <Quantity.hpp>
#include <types.hpp>

namespace mbsolve {

class ISource
{
public:
    virtual Quantity operator()(const Quantity& x) const = 0;
    /* TODO: add position. how?
	   virtual const Quantity& position
    */
    /* TODO: add source type: hard, soft, thevenin */
    /* TODO: for thevenin: add internal resistance */
};

class SineSource : ISource
{
private:
    Quantity m_ampl;
    Quantity m_freq;
    Quantity m_phase;
public:
    SineSource(const Quantity& ampl, const Quantity& freq,
	       const Quantity& phase) :
	m_ampl(ampl), m_freq(freq), m_phase(phase)
    {
    }

    Quantity operator()(const Quantity& x) const
    {
	/* TODO fill in */
	return m_ampl;
    }

};

class GaussianPulse : ISource
{

public:


};

}

#endif
