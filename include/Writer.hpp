#ifndef WRITER_H
#define WRITER_H

#include <string>
#include <vector>
#include <Device.hpp>
#include <Scenario.hpp>
#include <Types.hpp>

namespace mbsolve {

class Writer
{
public:
    virtual void write(const std::string& file,
		       const std::vector<Result *>& results,
		       const Device& device,
		       const Scenario& scenario)
    {
    }

};

}

#endif
