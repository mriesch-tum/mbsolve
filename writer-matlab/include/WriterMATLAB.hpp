#ifndef WRITERMATLAB_H
#define WRITERMATLAB_H

#include <Writer.hpp>

namespace mbsolve {

class WriterMATLAB : public Writer
{
public:
    void write(const std::string& file, const std::vector<Result *>& results,
	       const Device& device, const Scenario& scenario);

};

}

#endif
